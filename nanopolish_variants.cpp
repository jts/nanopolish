//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variants -- tools for calling variants
//
#include <algorithm>
#include "nanopolish_profile_hmm.h"
#include "nanopolish_variants.h"
#include "nanopolish_haplotype.h"
#include "overlapper.h"

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
    // Align the reference sequence and haplotype
    OverlapperParams params = affine_default_params;
    params.type = ALT_GLOBAL;
    SequenceOverlap overlap = Overlapper::computeAlignmentAffine(reference, haplotype, params);
    
    // Make aligned strings where gaps are padded with '-'
    std::string pad_ref;
    std::string pad_hap;
    overlap.makePaddedMatches(reference, haplotype, &pad_ref, &pad_hap);
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
            is_indel = pad_ref[diff_end] == '-' || pad_hap[diff_end] == '-';
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
    return variants;
}

// Apply a variant to the reference sequence to generate a new sequence
std::string apply_variant(const std::string& reference,
                          Variant v)
{
    std::string tmp = reference;
    assert(reference.substr(v.ref_position, v.ref_seq.length()) == v.ref_seq);
    tmp.replace(v.ref_position, v.ref_seq.length(), v.alt_seq);
    return tmp;
}

// Parse variants from the called haplotype and calculate
// quality scores for them
std::vector<Variant> evaluate_variants(const std::string& reference, 
                                       const std::string& haplotype, 
                                       const std::vector<HMMInputData>& input)
{
    // Calculate baseline probablilty

    std::vector<Variant> all_variants = extract_variants(reference, haplotype);
    std::vector<Variant> selected_variants;

    //
    // Test each variant, greedily selecting the best one at each step
    //
    
    Haplotype base(reference);
    double base_lp = profile_hmm_score(base.get_sequence(), input);

    while(!all_variants.empty()) {
 
        double best_variant_lp = -INFINITY;
        size_t best_variant_idx = 0;

        for(size_t i = 0; i < all_variants.size(); ++i) {
        
            // apply the variant to get a new haplotype
            Variant& v = all_variants[i];
            Haplotype derived = base;
            derived.apply_variant(v);

            // score the haplotype
            double variant_lp = profile_hmm_score(derived.get_sequence(), input);
            
            if(variant_lp > best_variant_lp) {
                best_variant_lp = variant_lp;
                best_variant_idx = i;
            }
        }

        if(best_variant_lp > base_lp) {
            // move the best variant from the all list to the selected list
            selected_variants.push_back(all_variants[best_variant_idx]);

            all_variants.erase(all_variants.begin() + best_variant_idx);
            
            // calculate a quality score for the variant
            selected_variants.back().quality = best_variant_lp - base_lp;

            //printf("SELECTED %zu from %zu: \n\t", best_variant_idx, all_variants.size() + 1);
            //selected_variants.back().write_vcf(stdout);

            // apply the variant to the base haplotype
            base.apply_variant(selected_variants.back());
            base_lp = best_variant_lp;
        } else {
            // no variant improved upon the base haplotype, stop
            break;
        }
    }

    /*
    // Calculate qualities for the variants in isolation
    // TODO: make this haplotype based
    for(size_t i = 0; i < variants.size(); ++i) {
        Variant& v = variants[i];
        std::string derived = apply_variant(reference, v);
        double lp_derived = profile_hmm_score(derived, input);
        v.quality = lp_derived - lp_base;
        v.quality = std::max(0.0, v.quality);
    }
    */
    return selected_variants;
}

