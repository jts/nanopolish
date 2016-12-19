//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variant_db -- data structure for working with
// a collection of variants
//
#ifndef NANOPOLISH_VARIANT_DB_H
#define NANOPOLISH_VARIANT_DB_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include "nanopolish_variant.h"

// typedefs
typedef size_t VariantGroupID;

enum CombinationOption
{
    CO_WITHOUT_REPLACEMENT,
    CO_WITH_REPLACEMENT
};

typedef std::vector<std::vector<size_t>> SizeTVecVec;

SizeTVecVec cartesian_product(const SizeTVecVec& input);

// Algorithm from: http://stackoverflow.com/questions/9430568/generating-combinations-in-c
class Combinations
{
    public:
        Combinations(size_t N, size_t k, CombinationOption option = CO_WITHOUT_REPLACEMENT);

        // Have all combinations been processed?
        bool done();

        // Get the current combination, as a vector of indices
        std::vector<size_t> get() const;
        std::string get_as_string() const;
        
        // Get the rank of the current combination (the rank of the N-th combination generated is N)
        size_t get_rank() const { return m_rank; }
        // Go to the next combination
        void next();

    private:

        std::vector<size_t> _get_without_replacement() const;
        std::vector<size_t> _get_with_replacement() const;
        
        size_t m_rank;
        bool m_done;
        CombinationOption m_option;
        std::vector<bool> setmask;
        std::vector<bool> initial_setmask;
};

// A variant combination is a subset of n variants
// from a VariantGroup. The Variants are represented
// by a vector of IDs.
class VariantCombination
{
    public:
        VariantCombination(const std::vector<size_t>& variant_ids);
        
        size_t get_num_variants() const { return m_variant_ids.size(); }

        // Get the ID of the variant at idx in the combination
        size_t get_variant_id(size_t idx) const { return m_variant_ids[idx]; }

    private:
        std::vector<size_t> m_variant_ids;
};

// A collection of variants from a region
// of a genome that are close enough together to
// require joint calling
class VariantGroup
{

    public:
        VariantGroup(VariantGroupID id, const std::vector<Variant>& v) : m_group_id(id), m_variants(v) {}

        // Return the ID of this group
        VariantGroupID getID() const { return m_group_id; }
        
        //
        // Variant access
        //

        // get a single variant
        const Variant& get(size_t i) const { return m_variants[i]; }
        
        // get a subset of variants according to the input combination
        const std::vector<Variant> get_variants(const VariantCombination& vc) const;
        
        // Return the number of variants in the group
        size_t get_num_variants() const { return m_variants.size(); }
    
        // Add a new variant combination to the collection. Returns its index.
        size_t add_combination(const VariantCombination& vc);
        const VariantCombination& get_combination(size_t idx) const { return m_combinations[idx]; }
        size_t get_num_combinations() const { return m_combinations.size(); }
        std::string get_vc_allele_string(size_t idx) const;

        // Set the score computed by the HMM for a variant combination for a single read
        void set_combination_read_score(size_t combination_idx, const std::string& read_id, double score);
        double get_combination_read_score(size_t combination_idx, const std::string& read_id) const;

        // Return the IDs and sum scores of all reads used in this group
        std::vector<std::pair<std::string, double>> get_read_sum_scores() const;

    private:

        VariantGroupID m_group_id;
        std::vector<Variant> m_variants;

        typedef std::map<std::string, double> ReadScoreMap;
        std::vector<VariantCombination> m_combinations;
        std::vector<ReadScoreMap> m_scores;
        std::map<std::string, double> m_read_score_sum;
};

class VariantDB
{
    public:
        VariantDB() {}

        // Add a new variant group into the collection
        // Returns the (numeric) ID of this group
        size_t add_new_group(const std::vector<Variant>& variants);

        // 
        size_t get_num_groups() const { return m_variant_groups.size(); }

        //
        VariantGroup& get_group(size_t i) { return m_variant_groups[i]; }

    private:
        std::vector<VariantGroup> m_variant_groups;
};

#endif
