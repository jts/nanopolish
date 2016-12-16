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

void score_variant_group(VariantGroup& variant_group,
                         Haplotype base_haplotype, 
                         const std::vector<HMMInputData>& input,
                         const int max_haplotypes,
                         const int ploidy,
                         const bool genotype_all_input_variants,
                         const uint32_t alignment_flags)
{
    size_t num_variants = variant_group.get_num_variants();

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
    // Track the variant combination ID within the group
    std::vector<std::pair<Haplotype, size_t>> haplotypes;

    for(size_t r = 0; r <= max_r; ++r) {

        Combinations combinations(num_variants, r);
        while(!combinations.done()) {
            VariantCombination vc(combinations.get());

            // Apply variants to haplotype
            Haplotype current_haplotype = base_haplotype;
            bool good_haplotype = current_haplotype.apply_variants(variant_group.get_variants(vc));
            if(good_haplotype) {
                size_t vc_idx = variant_group.add_combination(vc);
                haplotypes.push_back(std::make_pair(current_haplotype, vc_idx));
            }
            combinations.next();
        }
    }

    std::vector<std::string> read_ids;
    for(size_t i = 0; i < input.size(); ++i) {
        std::stringstream ss;
        ss << input[i].read->read_name << ":" << input[i].strand;
        read_ids.push_back(ss.str());
    }

/*
    // Score each haplotype
    DoubleMatrix read_haplotype_scores;
    allocate_matrix(read_haplotype_scores, input.size(), haplotypes.size());

    // Score all reads against all haplotypes
    std::vector<double> read_sum(input.size(), -INFINITY);
*/
  
    #pragma omp parallel for
    for(size_t ri = 0; ri < input.size(); ++ri) {
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            const auto& current = haplotypes[hi];
            double score = profile_hmm_score(current.first.get_sequence(), input[ri], alignment_flags);
            
            #pragma omp critical
            {
                variant_group.set_combination_read_score(current.second, read_ids[ri], score);
//                set(read_haplotype_scores, ri, hi, score);
//                read_sum[ri] = add_logs(read_sum[ri], score);
            }
        }
    }

#if 0
#endif
}

std::vector<Variant> simple_call(VariantGroup& variant_group,
                                 const int ploidy,
                                 const bool genotype_all_input_variants)
{
    // Select the best set of haplotypes that maximizes the probability of the data
    double base_score = -INFINITY;
    double best_score = -INFINITY;
    std::vector<size_t> best_set;
    std::vector<size_t> base_set;

#ifdef DEBUG_HAPLOTYPE_SELECTION 
    fprintf(stderr, "Selecting haplotypes\n");
#endif
    
    // 
    const std::vector<std::string> group_reads = variant_group.get_read_ids();

    Combinations vc_sets(variant_group.get_num_combinations(), ploidy, CO_WITH_REPLACEMENT);
    while(!vc_sets.done()) {

        // The current combination is represented as a vector of haplotype IDs
        std::vector<size_t> current_set = vc_sets.get();
        
        // Check if the current set consists of entirely of haplotypes without variants
        bool is_base_set = true;
        for(size_t i = 0; i < current_set.size(); ++i) {
            const VariantCombination& vc = variant_group.get_combination(current_set[i]);
            is_base_set = is_base_set && (variant_group.get_variants(vc).size() == 0);
        }

        double set_score = 0.0f;
        //std::vector<double> read_support(current_haplotypes.size(), 0.0f);

        for(const auto& read_id : group_reads) {
            double set_sum = -INFINITY;
            for(size_t j = 0; j < current_set.size(); ++j) {
                size_t vc_id = current_set[j];
                double rhs = variant_group.get_combination_read_score(vc_id, read_id);
                set_sum = add_logs(set_sum, rhs);
                //read_support[j] += exp(rhs - read_sum[ri]);
            }
            set_score += set_sum;
        }
        
        if(is_base_set) {
            base_score = set_score;
            base_set = current_set;
        }
        
#ifdef DEBUG_HAPLOTYPE_SELECTION 
        fprintf(stderr, "Current set score: %.5lf\t", set_score);
        for(size_t i = 0; i < current_set.size(); ++i) {
            fprintf(stderr, "\t%zu:%.2lf", current_set[i], /*read_support[i]*/0.0f);
        }
        fprintf(stderr, "\n");
#endif

        if(set_score > best_score) {
            best_score = set_score;
            best_set = current_set;
        }
        vc_sets.next();
    }

    // TODO: set an appropriate threshold
    if(best_score - base_score < 5) {
        best_set = base_set;
    }

    std::vector<Variant> output_variants;
    for(size_t vi = 0; vi < variant_group.get_num_variants(); vi++) {

        size_t var_count = 0;
        for(size_t j = 0; j < best_set.size(); ++j) {
            const VariantCombination& curr_vc = variant_group.get_combination(best_set[j]);

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
        v.add_info("TotalReads", group_reads.size());
        v.add_info("AlleleCount", var_count);
        v.genotype = make_genotype(var_count, ploidy);
        output_variants.push_back(v);
    }

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

