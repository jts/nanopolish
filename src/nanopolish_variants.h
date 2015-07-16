//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variants -- tools for calling variants
//
#ifndef NANOPOLISH_VARIANTS_H
#define NANOPOLISH_VARIANTS_H

#include <sstream>
#include "stdaln.h"
#include "nanopolish_common.h"

struct Variant
{
    static void write_vcf_header(FILE* fp) {
        fprintf(fp, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
    }

    // generate a unique identifier for this variant
    std::string key() const
    {
        std::stringstream out;
        out << ref_name << ':' << ref_position << ':' << ref_seq << ':' << alt_seq;
        return out.str();
    }

    void write_vcf(FILE* fp)
    {
        fprintf(fp, "%s\t%zu\t%s\t", ref_name.c_str(), ref_position, ".");
        fprintf(fp, "%s\t%s\t%.1lf\t", ref_seq.c_str(), alt_seq.c_str(), quality);
        fprintf(fp, "%s\t%s\n", "PASS", info.c_str());
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

// Determine potential variants between the reference and haplotype string
std::vector<Variant> extract_variants(const std::string& reference, 
                                      const std::string& haplotype);

// Remove redundant variants from the vector
// Optionally remove variants that have been seen fewer than min_count times
void deduplicate_variants(std::vector<Variant>& variants);

// Parse variants from the called haplotype and calculate
// quality scores for them
std::vector<Variant> evaluate_variants(const std::string& reference, 
                                       const std::string& haplotype, 
                                       const std::vector<HMMInputData>& input);

#endif
