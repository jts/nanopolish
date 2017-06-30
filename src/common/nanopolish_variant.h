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
class VariantGroup;

struct Variant
{
    static void write_vcf_header(FILE* fp, 
                                 const std::vector<std::string>& tag_lines = std::vector<std::string>());

    static std::string make_vcf_tag_string(const std::string& tag,
                                           const std::string& id,
                                           int count,
                                           const std::string& type,
                                           const std::string& description);

    Variant() { }
    Variant(const std::string& line) { read_vcf(line); }

    // generate a unique identifier for this variant
    std::string key() const
    {
        std::stringstream out;
        out << ref_name << ':' << ref_position << ':' << ref_seq << ':' << alt_seq;
        return out.str();
    }

    void write_vcf(FILE* fp) const
    {
        assert(fp != NULL);
        const char* gt_def = genotype.empty() ? NULL : "GT";
        const char* gt_str = genotype.empty() ? NULL : genotype.c_str();

        fprintf(fp, "%s\t%zu\t%s\t", ref_name.c_str(), ref_position + 1, ".");
        fprintf(fp, "%s\t%s\t%.1lf\t", ref_seq.c_str(), alt_seq.c_str(), quality);
        fprintf(fp, "%s\t%s\t%s\t%s\n", "PASS", info.c_str(), gt_def, gt_str);
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
        ss >> dummy; // GT tag
        ss >> genotype;

        // VCF is 1-based but we internally represent a variant as 0-based
        ref_position -= 1;

        assert(!ref_name.empty());
        assert(!ref_seq.empty());
        assert(!alt_seq.empty());
        //assert(ref_position >= 0);
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

    bool is_snp() const { return ref_seq.length() == 1 && alt_seq.length() == 1; }

    std::string ref_name;
    size_t ref_position;
    std::string ref_seq;
    std::string alt_seq;
    double quality;
    std::string info;
    std::string genotype;
};

inline bool sortByPosition(const Variant& a, const Variant& b) 
{ 
    return a.ref_name == b.ref_name ? 
        a.ref_position < b.ref_position : 
        a.ref_name < b.ref_name; 
}

class VariantKeyComp
{
    public: 
        inline bool operator()(const Variant& a, const Variant& b)
        {
            return a.key() < b.key();
        }
};

// Read a collection of variants from a VCF file
std::vector<Variant> read_variants_from_file(const std::string& filename);
std::vector<Variant> read_variants_for_region(const std::string& filename,
                                              const std::string& contig,
                                              int region_start,
                                              int region_end);

// Remove variants that are in the vector fewer than min_count times
void filter_variants_by_count(std::vector<Variant>& variants, int min_count);

// Remove snps or indels 
void filter_out_non_snp_variants(std::vector<Variant>& variants);

// Score the variants contained within the input group using the nanopolish HMM
void score_variant_group(VariantGroup& variant_group,
                         Haplotype base_haplotype, 
                         const std::vector<HMMInputData>& input,
                         const int max_haplotypes,
                         const int ploidy,
                         const bool genotype_all_input_variants,
                         const uint32_t alignment_flags);

// Call genotypes for the variants in this group using a simple model
std::vector<Variant> simple_call(VariantGroup& variant_group,
                                 const int ploidy,
                                 const bool genotype_all_input_variants);

// Call genotypes for the variants in this group using support from nearby variants
std::vector<Variant> multi_call(VariantGroup& variant_group,
                                std::vector<const VariantGroup*> neighbor_groups,
                                const int ploidy,
                                const bool genotype_all_input_variants);

// Select variants that have a positive score wrt the base haplotype
std::vector<Variant> select_positive_scoring_variants(std::vector<Variant>& candidate_variants,
                                                      Haplotype base_haplotype, 
                                                      const std::vector<HMMInputData>& input,
                                                      const uint32_t alignment_flags);

// Score a single variant, stopping when the absolute value of the score relative
// to the reference meets a threshold
Variant score_variant_thresholded(const Variant& input_variant,
                                  Haplotype base_haplotype, 
                                  const std::vector<HMMInputData>& input,
                                  const uint32_t alignment_flags,
                                  const uint32_t score_threshold);

#endif
