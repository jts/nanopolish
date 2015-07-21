//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variants -- tools for calling variants
//
#include <algorithm>
#include <map>
#include "nanopolish_profile_hmm.h"
#include "nanopolish_variants.h"
#include "nanopolish_haplotype.h"

// return a new copy of the string with gap symbols removed
std::string remove_gaps(const std::string& str)
{
    std::string ret = str;
    ret.erase( std::remove(ret.begin(), ret.end(), '-'), ret.end());
    return ret;
}

// extract differences between the pair of strings
std::vector<Variant> extract_variants(const std::string& reference, 
                                      const std::string& haplotype)
{
    AlnParam par = aln_param_nt2nt;
    par.band_width = std::max(20, abs(reference.size() - haplotype.size()) * 2);
    AlnAln* aln = aln_stdaln(reference.c_str(), haplotype.c_str(), &par, 1, 1);
    
    // Make aligned strings where gaps are padded with '-'
    std::string pad_ref(aln->out1);
    std::string pad_hap(aln->out2);

    assert(pad_ref.size() == pad_hap.size());
    
    //std::cout << "PR: " << pad_ref << "\n";
    //std::cout << "PH: " << pad_hap << "\n";

    // parse variants from the alignment
    std::vector<Variant> variants;

    // generate a map from padded bases to positions in the original reference sequence
    std::vector<size_t> ref_positions(pad_ref.size(), 0);
    size_t pos = 0;
    for(size_t i = 0; i < pad_ref.size(); ++i) {
        ref_positions[i] = pad_ref[i] != '-' ? pos : std::string::npos;
        pos += pad_ref[i] != '-';
    }

    // diff_start iterates over the places where these sequences are different
    size_t diff_start = 0;
    while(1) {
        
        // find the start point of the next difference between the strings
        while(diff_start < pad_ref.size() && pad_ref[diff_start] == pad_hap[diff_start]) {
            diff_start++;
        }
 
        // check for end of alignment
        if(diff_start == pad_ref.size())
            break;

        // find the end point of the difference
        bool is_indel = false;
        size_t diff_end = diff_start;
        while(diff_end < pad_ref.size() && pad_ref[diff_end] != pad_hap[diff_end]) {
            is_indel = is_indel || pad_ref[diff_end] == '-' || pad_hap[diff_end] == '-';
            diff_end++;
        }

        // If the difference is an indel, we include the previous matching reference base
        diff_start -= is_indel;
    
        Variant v;
        v.ref_name = "noctg";

        assert(ref_positions[diff_start] != std::string::npos);
        v.ref_position = ref_positions[diff_start];
        v.ref_seq = remove_gaps(pad_ref.substr(diff_start, diff_end - diff_start).c_str());
        v.alt_seq = remove_gaps(pad_hap.substr(diff_start, diff_end - diff_start).c_str());
        
        variants.push_back(v);
        diff_start = diff_end;
    }

    aln_free_AlnAln(aln);
    return variants;
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
            v.add_info("ReadDepth", iter->second.second);
            variants.push_back(v);
        }
    }
}

void filter_out_non_snp_variants(std::vector<Variant>& variants)
{
    std::vector<Variant> tmp;
    for(size_t i = 0; i < variants.size(); ++i) {
        bool is_snp = variants[i].ref_seq.length() == 1 &&
                      variants[i].alt_seq.length() == 1;
        if(is_snp) {
            tmp.push_back(variants[i]);
        }
    }
    variants.swap(tmp);
}

// Parse variants from the called haplotype and calculate
// quality scores for them
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

        if(best_variant_lp - base_lp > 1.0) {
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

