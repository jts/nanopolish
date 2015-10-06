//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variant -- tools for calling variants
//
#ifndef NANOPOLISH_VARIANT_H
#define NANOPOLISH_VARIANT_H

#include <sstream>
#include "stdaln.h"
#include "nanopolish_common.h"

// forward declare
class Haplotype;

struct Variant
{
    static void write_vcf_header(FILE* fp) {
        fprintf(fp, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
    }

    Variant() { }
    Variant(const std::string& line) { read_vcf(line); }

    // generate a unique identifier for this variant
    std::string key() const
    {
        std::stringstream out;
        out << ref_name << ':' << ref_position << ':' << ref_seq << ':' << alt_seq;
        return out.str();
    }

    void write_vcf(FILE* fp)
    {
        fprintf(fp, "%s\t%zu\t%s\t", ref_name.c_str(), ref_position + 1, ".");
        fprintf(fp, "%s\t%s\t%.1lf\t", ref_seq.c_str(), alt_seq.c_str(), quality);
        fprintf(fp, "%s\t%s\n", "PASS", info.c_str());
    }

    void read_vcf(const std::string& line)
    {
        std::stringstream ss(line);
        std::string dummy;
        ss >> ref_name;
        ss >> ref_position;
        ss >> dummy; // ID, not used
        ss >> ref_seq;
        ss >> alt_seq;
        ss >> quality;
        ss >> dummy; // FILTER, not used
        ss >> info;

        // VCF is 1-based but we internally represent a variant as 0-based
        ref_position -= 1;

        assert(!ref_name.empty());
        assert(!ref_seq.empty());
        assert(!alt_seq.empty());
        assert(ref_position >= 0);
        assert(quality >= 0.0f);
    }

    template<typename T>
    void add_info(const std::string& key, T value)
    {
        std::stringstream ss;
        ss << key << "=" << value;
        if(info.empty()) {
            info = ss.str();
        } else {
            info.append(1, ';');
            info.append(ss.str());
        }
    }

    std::string ref_name;
    size_t ref_position;
    std::string ref_seq;
    std::string alt_seq;
    double quality;
    std::string info;
};

inline bool sortByPosition(const Variant& a, const Variant& b) 
{ 
    return a.ref_name == b.ref_name ? 
        a.ref_position < b.ref_position : 
        a.ref_name < b.ref_name; 
}

// Determine potential variants between the reference and haplotype string
std::vector<Variant> extract_variants(const std::string& reference, 
                                      const std::string& haplotype);

// Remove variants that are in the vector fewer than min_count times
void filter_variants_by_count(std::vector<Variant>& variants, int min_count);

// Remove snps or indels 
void filter_out_non_snp_variants(std::vector<Variant>& variants);

// Select variants to add to the base haplotype one-by-one
std::vector<Variant> select_variants(const std::vector<Variant>& candidate_variants,
                                     Haplotype base_haplotype, 
                                     const std::vector<HMMInputData>& input);

// Select groups of variants to add to the base haplotype
std::vector<Variant> select_variant_set(const std::vector<Variant>& candidate_variants,
                                        Haplotype base_haplotype, 
                                        const std::vector<HMMInputData>& input,
                                        const uint32_t alignment_flags);


#endif
