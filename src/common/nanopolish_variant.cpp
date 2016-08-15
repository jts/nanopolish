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

std::vector<Variant> select_variant_set(const std::vector<Variant>& candidate_variants,
                                        Haplotype base_haplotype, 
                                        const std::vector<HMMInputData>& input,
                                        const int max_haplotypes,
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

    // Calculate the likelihood of the haplotype with no additional variants added
    // also do some bookkeeping about per-read/per-model likelihoods
    double base_lp = 0.0f;
    double base_lp_by_model_strand[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    int read_counts[6] = { 0, 0, 0, 0, 0, 0 };

    std::vector<double> base_lp_by_read(input.size()); 
    
    #pragma omp parallel for
    for(size_t j = 0; j < input.size(); ++j) {

        double tmp = profile_hmm_score(base_haplotype.get_sequence(), input[j], alignment_flags);

        #pragma omp critical
        {
            base_lp_by_read[j] = tmp;
            base_lp += tmp;

            int mid = input[j].read->pore_model[input[j].strand].metadata.model_idx;
            int cid = 2 * mid + input[j].rc;
            base_lp_by_model_strand[cid] += tmp;
            read_counts[cid] += 1;
        }
    }

    double best_lp = -INFINITY;
    std::vector<Variant> best_variant_set;

    // Score haplotypes by adding 1, 2, ..., max_r variant sets to it
    for(size_t r = 1; r <= max_r; ++r) {
        // From: http://stackoverflow.com/questions/9430568/generating-combinations-in-c
        std::vector<bool> variant_selector(num_variants);
        std::fill(variant_selector.begin(), variant_selector.begin() + r, true);

        do {
            Haplotype current_haplotype = base_haplotype;
            std::vector<Variant> current_variant_set;
            bool good_haplotype = true;

            for(size_t vi = 0; vi < num_variants; vi++) {
                if(!variant_selector[vi]) {
                    continue;
                }

                current_variant_set.push_back(candidate_variants[vi]);
                good_haplotype = good_haplotype && current_haplotype.apply_variant(current_variant_set.back());
            }
            
            // skip the haplotype if all the variants couldnt be added to it
            if(!good_haplotype) {
                continue;
            }

            // score the haplotype
            double current_lp = 0.0f;
            double current_lp_by_model_strand[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
            size_t supporting_reads = 0;
            std::vector<double> relative_lp_by_read(input.size(), 0.0f);

            #pragma omp parallel for
            for(size_t j = 0; j < input.size(); ++j) {
                double tmp = profile_hmm_score(current_haplotype.get_sequence(), input[j], alignment_flags);

                #pragma omp critical
                {
                    current_lp += tmp;
                    supporting_reads += tmp > base_lp_by_read[j];
                    int mid = input[j].read->pore_model[input[j].strand].metadata.model_idx;
                    int cid = 2 * mid + input[j].rc;
                    current_lp_by_model_strand[cid] += tmp;
                    relative_lp_by_read[j] = tmp - base_lp_by_read[j];
                }
            }

            if(current_lp > best_lp && current_lp - base_lp > 0.1) {
                best_lp = current_lp;
                best_variant_set = current_variant_set;

                // Annotate variants
                for(size_t vi = 0; vi < best_variant_set.size(); ++vi) {
                    Variant& v = best_variant_set[vi];
                    v.add_info("TotalReads", input.size());
                    v.add_info("SupportingReads", supporting_reads);
                    v.add_info("SupportFraction", (double)supporting_reads / input.size());

                    // Annotate variants with qualities from the three possible models
                    std::string names[3] = { "Template", "Comp.P1", "Comp.P2" };

                    for(int mid = 0; mid < 3; mid++) {
                        int cid = 2 * mid;
                        double s0 = current_lp_by_model_strand[cid] - base_lp_by_model_strand[cid];
                        int c0 = read_counts[cid];

                        double s1 = current_lp_by_model_strand[cid + 1] - base_lp_by_model_strand[cid + 1];
                        int c1 = read_counts[cid + 1];

                        std::stringstream ss;
                        ss << std::setprecision(4) << s0 / c0 << "," << s1 / c1;
                        v.add_info(names[mid], ss.str());
                    }

                    /*
                    v.add_info("TemplateQuality", current_lp_by_strand[0] - base_lp_by_strand[0]);
                    v.add_info("ComplementQuality", current_lp_by_strand[1] - base_lp_by_strand[1]);
                    v.add_info("ForwardQuality", current_lp_by_rc[0] - base_lp_by_rc[0]);
                    v.add_info("ReverseQuality", current_lp_by_rc[1] - base_lp_by_rc[1]);
                    v.add_info("TAvgQuality", (current_lp_by_model[0] - base_lp_by_model[0]) / model_count[0]);
                    v.add_info("C1AvgQuality", (current_lp_by_model[1] - base_lp_by_model[1]) / model_count[1]);
                    v.add_info("C2AvgQuality", (current_lp_by_model[2] - base_lp_by_model[2]) / model_count[2]);
                    */

                    std::stringstream counts;
                    std::ostream_iterator<int> rc_out(counts, ",");
                    std::copy(std::begin(read_counts), std::end(read_counts), rc_out);
                    std::string rc_str = counts.str();
                    v.add_info("ReadCounts", rc_str.substr(0, rc_str.size() - 1));

                    std::stringstream scores;
                    std::ostream_iterator<float> scores_out(scores, ",");
                    std::copy(std::begin(relative_lp_by_read), std::end(relative_lp_by_read), scores_out);
                    std::string scores_str = scores.str();
                    v.add_info("Scores", scores_str.substr(0, scores_str.size() - 1));
                    
                    v.quality = best_lp - base_lp;
                }
            }
#ifdef DEBUG_HAPLOTYPE_SELECTION
            std::stringstream ss;
            for(size_t vi = 0; vi < current_variant_set.size(); ++vi) {
                const Variant& v = current_variant_set[vi];
                ss << (vi > 0 ? "," : "") << v.key();
            }
            fprintf(stderr, "haplotype: %zu variants: %s relative score: %.2lf\n", hi, ss.str().c_str(), current_lp - base_lp);
#endif
        } while(std::prev_permutation(variant_selector.begin(), variant_selector.end()));
    }
    return best_variant_set;
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

