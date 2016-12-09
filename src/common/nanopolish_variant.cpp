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

//#define DEBUG_HAPLOTYPE_SELECTION 1

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

std::vector<Variant> select_variants(const std::vector<Variant>& candidate_variants,
                                     Haplotype base_haplotype,
                                     const std::vector<HMMInputData>& input)
{
    // make a copy of the variant set to modify
    std::vector<Variant> all_variants = candidate_variants;

    // Calculate baseline probablilty
    std::vector<Variant> selected_variants;

    double base_lp = profile_hmm_score(base_haplotype.get_sequence(), input);

    while(!all_variants.empty()) {
 
        double best_variant_lp = -INFINITY;
        size_t best_variant_idx = 0;
        size_t best_supporting_reads = 0;

        std::vector<double> base_lp_by_read; 
        for(size_t j = 0; j < input.size(); ++j) {
            double tmp = profile_hmm_score(base_haplotype.get_sequence(), input[j]);
            base_lp_by_read.push_back(tmp);
        }

        for(size_t i = 0; i < all_variants.size(); ++i) {
        
            // apply the variant to get a new haplotype
            Variant& v = all_variants[i];
            Haplotype derived = base_haplotype;
            derived.apply_variant(v);

            // score the haplotype
            double variant_lp = 0.0f;
            size_t supporting_reads = 0;
            #pragma omp parallel for
            for(size_t j = 0; j < input.size(); ++j) {
                double tmp = profile_hmm_score(derived.get_sequence(), input[j]);
                #pragma omp critical
                {
                    variant_lp += tmp;
                    supporting_reads += tmp > base_lp_by_read[j];
                }
            }
            
            if(variant_lp > best_variant_lp) {
                best_variant_lp = variant_lp;
                best_variant_idx = i;
                best_supporting_reads = supporting_reads;
            }
        }

        if(best_variant_lp - base_lp > 0.1) {
            // move the best variant from the all list to the selected list
            Variant& best_variant = all_variants[best_variant_idx];
            best_variant.add_info("TotalReads", input.size());
            best_variant.add_info("SupportingReads", best_supporting_reads);
            best_variant.add_info("SupportFraction", (double)best_supporting_reads / input.size());
            best_variant.quality = best_variant_lp - base_lp;
            selected_variants.push_back(best_variant);

            //printf("SELECTED %zu from %zu: \n\t", best_variant_idx, all_variants.size() + 1);

            // apply the variant to the base haplotype
            base_haplotype.apply_variant(selected_variants.back());
            base_lp = best_variant_lp;
        } else {
            // no variant improved upon the base haplotype, stop
            break;
        }
    }

    return selected_variants;
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

std::vector<Variant> call_variants(const std::vector<Variant>& candidate_variants,
                                   Haplotype base_haplotype, 
                                   const std::vector<HMMInputData>& input,
                                   const int max_haplotypes,
                                   const int ploidy,
                                   const bool genotype_all_input_variants,
                                   const uint32_t alignment_flags)
{
    size_t num_variants = candidate_variants.size();

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
    std::vector< std::vector<bool> > haplotype_variant_mask;
    std::vector<bool> is_haplotype_valid;

    for(size_t r = 0; r <= max_r; ++r) {
        // From: http://stackoverflow.com/questions/9430568/generating-combinations-in-c
        std::vector<bool> variant_selector(num_variants);
        std::fill(variant_selector.begin(), variant_selector.begin() + r, true);

        do {
            Haplotype current_haplotype = base_haplotype;
            bool good_haplotype = true;
            for(size_t vi = 0; vi < num_variants; vi++) {
                if(!variant_selector[vi]) {
                    continue;
                }

                good_haplotype = good_haplotype && current_haplotype.apply_variant(candidate_variants[vi]);
            }
            
            // skip the haplotype if all the variants couldnt be added to it
            haplotypes.push_back(current_haplotype);
            is_haplotype_valid.push_back(good_haplotype);
            haplotype_variant_mask.push_back(variant_selector);

            if(!good_haplotype) {
                continue;
            }
        } while(std::prev_permutation(variant_selector.begin(), variant_selector.end()));
    }

    DoubleMatrix read_haplotype_scores;
    allocate_matrix(read_haplotype_scores, input.size(), haplotypes.size());

    // Score all reads against all haplotypes
    #pragma omp parallel for
    for(size_t ri = 0; ri < input.size(); ++ri) {
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            double score = 0.0f;
            if(is_haplotype_valid[hi]) {
                score = profile_hmm_score(haplotypes[hi].get_sequence(), input[ri], alignment_flags);
            }
            
            #pragma omp critical
            set(read_haplotype_scores, ri, hi, score);
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
    
    printf("Selecting haplotypes (%zu of %zu): \n", ploidy, haplotypes.size());
    do {

        // Convert haplotype mask into an array of haplotype indices
        std::vector<size_t> current_haplotypes;
        bool all_haplotypes_valid = true;

        size_t curr_hi = 0;
        for(const auto& bit : haplotype_selector) {
            if(bit) {
                assert(curr_hi < haplotypes.size());
                all_haplotypes_valid = all_haplotypes_valid && is_haplotype_valid[curr_hi];
                current_haplotypes.push_back(curr_hi);
            } else {
                curr_hi++;
            }
        }

        bool is_base_set = true;
/*
        printf("\thaplotypes: ");
        for(size_t i = 0; i < current_haplotypes.size(); ++i) {
            printf("%zu ", current_haplotypes[i]);
        }
        printf("\n");
*/
        for(size_t i = 0; i < current_haplotypes.size(); ++i) {
            is_base_set = is_base_set && (current_haplotypes[i] == 0);
        }

        if(!all_haplotypes_valid) {
            continue;
        }

        double haplotype_set_score = 0.0f;

        for(size_t ri = 0; ri < input.size(); ++ri) {
            double set_sum = -INFINITY;
            for(size_t j = 0; j < current_haplotypes.size(); ++j) {
                size_t hi = current_haplotypes[j];
                set_sum = add_logs(set_sum, get(read_haplotype_scores, ri, hi));
            }
            haplotype_set_score += set_sum;
        }

        if(is_base_set) {
            base_score = haplotype_set_score;
        }

        if(haplotype_set_score > best_score) {
            best_score = haplotype_set_score;
            best_haplotype_set = current_haplotypes;
//            best_supporting_reads = supporting_reads;
        }
    } while(std::prev_permutation(haplotype_selector.begin(), haplotype_selector.end()));

    // TODO: set an appropriate threshold
    if(best_score - base_score < 1) {
        best_haplotype_set.assign(ploidy, 0); // set the called haplotypes to be the base (0)
    }

    printf("Selected haplotypes: ");
    for(size_t i = 0; i < best_haplotype_set.size(); ++i) {
        printf("%zu ", best_haplotype_set[i]);
    }
    printf("%.2lf\n", best_score - base_score);

    /*
    // Calculate the set of variants to call genotypes on and return
    // -in variant discovery mode this is the union of all variants on called haplotypes
    // -in genotype mode this is the complete set of input variants
    std::vector<bool> output_variant_mask(candidate_variants.size(), false);
    if(genotype_all_input_variants) {
        output_variant_mask.assign(candidate_variants.size(), true);
    } else {
        for(size_t vi = 0; vi < num_variants; vi++) {
            for(size_t j = 0; j < best_haplotype_set.size(); ++j) {
                output_variant_mask[vi] = output_variant_mask[vi] || haplotype_variant_mask[best_haplotype_set[j]][vi];
            }
        }
    }
    */
    std::vector<Variant> output_variants;
    for(size_t vi = 0; vi < num_variants; vi++) {

        /*
        if(!output_variant_mask[vi]) {
            continue;
        }
        */

        size_t var_count = 0;
        for(size_t j = 0; j < best_haplotype_set.size(); ++j) {
            var_count += haplotype_variant_mask[best_haplotype_set[j]][vi];
        }

        if( !(genotype_all_input_variants || var_count > 0)) {
            continue;
        }

        Variant v = candidate_variants[vi];
        v.quality = best_score - base_score;
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

