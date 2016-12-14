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
#include <assert.h>
#include "nanopolish_variant.h"

// typedefs
typedef size_t VariantGroupID;

enum CombinationOption
{
    CO_WITHOUT_REPLACEMENT,
    CO_WITH_REPLACEMENT
};

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

        // Go to the next combination
        void next();

    private:

       std::vector<size_t> _get_without_replacement() const;
       std::vector<size_t> _get_with_replacement() const;

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
        VariantCombination(VariantGroupID group_id, const std::vector<size_t>& variant_ids);
        
        //
        size_t get_num_variants() const { return m_variant_ids.size(); }

        size_t get_variant_id(size_t idx) const { return m_variant_ids[idx]; }


    private:
        VariantGroupID m_group_id;
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

        const Variant& get(size_t i) const { return m_variants[i]; }
        
        const std::vector<Variant> get_variants(const VariantCombination& vc) const;

        // Return the number of variants in the group
        size_t size() const { return m_variants.size(); }

        VariantCombination get_next_combination();

    private:

        VariantGroupID m_group_id;
        std::vector<Variant> m_variants;
};

class VariantDB
{
    public:
        VariantDB();

        // Add a new variant group into the collection
        // Returns the (numeric) ID of this group
        size_t add_new_variant_group(const std::vector<Variant>& variants);

        // 
        size_t get_num_variant_groups() const { return m_variant_groups.size(); }

    private:
        std::vector<VariantGroup> m_variant_groups;
};

#endif
