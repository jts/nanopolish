//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variants -- tools for calling variants
//
#include <algorithm>
#include <map>
#include <iterator>
#include <iomanip>
#include "nanopolish_profile_hmm.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_model_names.h"
#include "nanopolish_variant_db.h"

#define DEBUG_HAPLOTYPE_SELECTION 1

// return a new copy of the string with gap symbols removed
std::string remove_gaps(const std::string& str)
{
    std::string ret = str;
    ret.erase( std::remove(ret.begin(), ret.end(), '-'), ret.end());
    return ret;
}

void filter_variants_by_count(std::vector<Variant>& variants, int min_count)
{
    std::map<std::string, std::pair<Variant, int>> map;

    for(size_t i = 0; i < variants.size(); ++i) {
        std::string key = variants[i].key();
        auto iter = map.find(key);
        if(iter == map.end()) {
            map.insert(std::make_pair(key, std::make_pair(variants[i], 1)));
        } else {
            iter->second.second += 1;
        }
    }

    variants.clear();

    for(auto iter = map.begin(); iter != map.end(); ++iter) {
        Variant& v = iter->second.first;
        if(iter->second.second >= min_count) {
            v.add_info("BaseCalledReadsWithVariant", iter->second.second);
            variants.push_back(v);
        }
    }
}

void filter_out_non_snp_variants(std::vector<Variant>& variants)
{
    std::vector<Variant> tmp;
    for(size_t i = 0; i < variants.size(); ++i) {
        bool is_snp = variants[i].ref_seq.length() == 
                      variants[i].alt_seq.length();
        if(is_snp) {
            tmp.push_back(variants[i]);
        }
    }
    variants.swap(tmp);
}

// this doesn't handle triallele ("1/2") genotypes, yet
std::string make_genotype(size_t alt_alleles, size_t ploidy)
{
    assert(alt_alleles <= ploidy);
    std::string out(ploidy + ploidy - 1, '/');
    int num_ref = ploidy - alt_alleles;
    for(size_t pos = 0; pos < out.size(); pos += 2) {
        out[pos] = (pos / 2) < num_ref ? '0' : '1';
    }
    return out;
}

std::vector<Variant> call_variants(const VariantGroup& variant_group,
                                   Haplotype base_haplotype, 
                                   const std::vector<HMMInputData>& input,
                                   const int max_haplotypes,
                                   const int ploidy,
                                   const bool genotype_all_input_variants,
                                   const uint32_t alignment_flags)
{
    size_t num_variants = variant_group.size();

    // Determine the maximum number of variants we can jointly test
    // without exceeding the maximum number of haplotypes
    size_t sum_num_haplotypes = 0;
    size_t max_r = 1;
    
    while(max_r <= num_variants) {
        size_t num_haplotypes_r = nChoosek(num_variants, max_r);
        if(num_haplotypes_r + sum_num_haplotypes < max_haplotypes) {
            sum_num_haplotypes += num_haplotypes_r;
        } else {
            break;
        }
        //printf("n: %zu r: %zu nCr: %zu sum: %zu\n", num_variants, max_r, num_haplotypes_r, sum_num_haplotypes);
        max_r += 1;
    }
    max_r -= 1;

    // Construct haplotypes (including the base haplotype with no variants)
    std::vector<Haplotype> haplotypes;
    std::vector< VariantCombination > haplotype_variant_combinations;

    for(size_t r = 0; r <= max_r; ++r) {

        Combinations combinations(num_variants, r);
        while(!combinations.done()) {
            VariantCombination vc(variant_group.getID(), combinations.get());

            // Apply variants to haplotype
            Haplotype current_haplotype = base_haplotype;
            bool good_haplotype = current_haplotype.apply_variants(variant_group.get_variants(vc));
            if(good_haplotype) {
                haplotypes.push_back(current_haplotype);
                haplotype_variant_combinations.push_back(vc);
            }
            combinations.next();
        }
    }

    // Score each haplotype
    DoubleMatrix read_haplotype_scores;
    allocate_matrix(read_haplotype_scores, input.size(), haplotypes.size());

    // Score all reads against all haplotypes
    std::vector<double> read_sum(input.size(), -INFINITY);
    
    #pragma omp parallel for
    for(size_t ri = 0; ri < input.size(); ++ri) {
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            double score = profile_hmm_score(haplotypes[hi].get_sequence(), input[ri], alignment_flags);
            
            #pragma omp critical
            {
                set(read_haplotype_scores, ri, hi, score);
                read_sum[ri] = add_logs(read_sum[ri], score);
            }
        }
    }

    // Select the haplotype with the highest score relative to the base haplotype (no variants)
    double base_score = -INFINITY;
    double best_score = -INFINITY;
    double best_supporting_reads = 0;
    std::vector<size_t> best_haplotype_set;

    // Select the combination of haplotypes that maximizes the log-likelihood

    // Uses the bitvector approach from https://www.mathsisfun.com/combinatorics/combinations-permutations.html
    // (among many other places) to select haplotypes.
    std::vector<bool> haplotype_selector(haplotypes.size() + ploidy - 1);
    std::fill(haplotype_selector.begin(), haplotype_selector.begin() + ploidy, true);
   
#ifdef DEBUG_HAPLOTYPE_SELECTION 
    fprintf(stderr, "selecting haplotypes\n");
#endif

    do {
        // Convert haplotype mask into an array of haplotype indices
        std::vector<size_t> current_haplotypes;

        size_t curr_hi = 0;
        for(const auto& bit : haplotype_selector) {
            if(bit) {
                assert(curr_hi < haplotypes.size());
                current_haplotypes.push_back(curr_hi);
            } else {
                curr_hi++;
            }
        }
        
        bool is_hom = true;
        bool is_base_set = true;
        for(size_t i = 0; i < current_haplotypes.size(); ++i) {
            is_base_set = is_base_set && (current_haplotypes[i] == 0);
            is_hom = is_hom && (i == 0 || current_haplotypes[i] == current_haplotypes[i-1]);
        }

        double haplotype_set_score = 0.0f;
        std::vector<double> read_support(current_haplotypes.size(), 0.0f);

        for(size_t ri = 0; ri < input.size(); ++ri) {
            double set_sum = -INFINITY;
            for(size_t j = 0; j < current_haplotypes.size(); ++j) {
                size_t hi = current_haplotypes[j];
                double rhs = get(read_haplotype_scores, ri, hi);
                set_sum = add_logs(set_sum, rhs);
                read_support[j] += exp(rhs - read_sum[ri]);
            }
            haplotype_set_score += set_sum;
        }
        
        /*
        // HACK: penalize hets
        haplotype_set_score += is_hom ? 0.0 : -5.0;
        */

        if(is_base_set) {
            base_score = haplotype_set_score;
        }
        
#ifdef DEBUG_HAPLOTYPE_SELECTION 
        fprintf(stderr, "Haplotype set score: %.5lf\t", haplotype_set_score);
        for(size_t i = 0; i < current_haplotypes.size(); ++i) {
            fprintf(stderr, "\t%zu:%.2lf", current_haplotypes[i], read_support[i]);
        }
        fprintf(stderr, "\n");
#endif

        if(haplotype_set_score > best_score) {
            best_score = haplotype_set_score;
            best_haplotype_set = current_haplotypes;
        }
    } while(std::prev_permutation(haplotype_selector.begin(), haplotype_selector.end()));

    // TODO: set an appropriate threshold
    if(best_score - base_score < 5) {
        best_haplotype_set.assign(ploidy, 0); // set the called haplotypes to be the base (0)
    }

    std::vector<Variant> output_variants;
    for(size_t vi = 0; vi < num_variants; vi++) {

        size_t var_count = 0;
        for(size_t j = 0; j < best_haplotype_set.size(); ++j) {
            const VariantCombination& curr_vc = haplotype_variant_combinations[best_haplotype_set[j]];

            for(size_t k = 0; k < curr_vc.get_num_variants(); ++k) {
                var_count += curr_vc.get_variant_id(k) == vi;
            }
        }

        if( !(genotype_all_input_variants || var_count > 0)) {
            continue;
        }

        Variant v = variant_group.get(vi);
        if(var_count > 0) {
            v.quality = best_score - base_score;
        } else {
            v.quality = 0.0;
        }
        v.add_info("TotalReads", input.size());
        v.add_info("AlleleCount", var_count);
        v.genotype = make_genotype(var_count, ploidy);

        //v.add_info("SupportFraction", best_supporting_reads / input.size());

        output_variants.push_back(v);
    }

    free_matrix(read_haplotype_scores);
    return output_variants;
}

std::vector<Variant> select_positive_scoring_variants(std::vector<Variant>& candidate_variants,
                                                      Haplotype base_haplotype, 
                                                      const std::vector<HMMInputData>& input,
                                                      const uint32_t alignment_flags)
{
    std::vector<Variant> selected_variants;
    double base_score = 0.0f;
    #pragma omp parallel for
    for(size_t j = 0; j < input.size(); ++j) {

        double score = profile_hmm_score(base_haplotype.get_sequence(), input[j], alignment_flags);

        #pragma omp atomic
        base_score += score;
    }

    for(size_t vi = 0; vi < candidate_variants.size(); ++vi) {

        Haplotype current_haplotype = base_haplotype;
        current_haplotype.apply_variant(candidate_variants[vi]);
        
        double haplotype_score = 0.0f;
        #pragma omp parallel for
        for(size_t j = 0; j < input.size(); ++j) {
            double score = profile_hmm_score(current_haplotype.get_sequence(), input[j], alignment_flags);

            #pragma omp atomic
            haplotype_score += score;
        }

        if(haplotype_score > base_score) {
            candidate_variants[vi].quality = haplotype_score - base_score;
            selected_variants.push_back(candidate_variants[vi]);
        }
    }

    return selected_variants;
}

Variant score_variant(const Variant& input_variant,
                      Haplotype base_haplotype, 
                      const std::vector<HMMInputData>& input,
                      const uint32_t alignment_flags)
{
    Variant out_variant = input_variant;

    double base_score = 0.0f;
    #pragma omp parallel for
    for(size_t j = 0; j < input.size(); ++j) {

        double score = profile_hmm_score(base_haplotype.get_sequence(), input[j], alignment_flags);

        #pragma omp atomic
        base_score += score;
    }

    base_haplotype.apply_variant(input_variant);
        
    double haplotype_score = 0.0f;
#pragma omp parallel for
    for(size_t j = 0; j < input.size(); ++j) {
        double score = profile_hmm_score(base_haplotype.get_sequence(), input[j], alignment_flags);

#pragma omp atomic
        haplotype_score += score;
    }

    out_variant.quality = haplotype_score - base_score;
    return out_variant;
}

