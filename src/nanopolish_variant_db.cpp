//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variant_db -- data structures for working with
// collections of variants
//
#include <algorithm>
#include <sstream>
#include "nanopolish_variant_db.h"

Combinations::Combinations(size_t N, size_t k) 
{ 
    assert(k <= N);
    initial_setmask.resize(N, false);
    std::fill(initial_setmask.begin(), initial_setmask.begin() + k, true);
    setmask = initial_setmask;
    m_done = false;
}

bool Combinations::done()
{
    return m_done;
}

std::vector<size_t> Combinations::get() const
{
    assert(!m_done);
    std::vector<size_t> out;
    for(size_t i = 0; i < setmask.size(); ++i) {
        if(setmask[i]) {
            out.push_back(i);
        }
    }
    return out;
}

std::string Combinations::get_as_string() const
{
    std::vector<size_t> comb = get();

    std::stringstream ss;
    for(size_t i = 0; i < comb.size(); ++i) {

        ss << comb[i];
        if(i < comb.size() - 1) {
            ss << ",";
        }
    }
    return ss.str();
}

void Combinations::next()
{
    std::prev_permutation(setmask.begin(), setmask.end());
    m_done = (setmask == initial_setmask);
}

//
// VariantCombinations
//
VariantCombination::VariantCombination(VariantGroupID group_id, const std::vector<size_t>& variant_ids) :
                                                                                   m_group_id(group_id),
                                                                                   m_variant_ids(variant_ids)
{

}


//
// VariantGroup
//
const std::vector<Variant> VariantGroup::get_variants(const VariantCombination& vc) const
{
    std::vector<Variant> out;
    for(size_t i = 0; i < vc.get_num_variants(); ++i) {
        out.push_back(m_variants[vc.get_variant_id(i)]);
    }
    return out;
}


//
// VariantDB
//
size_t VariantDB::add_new_variant_group(const std::vector<Variant>& variants)
{
    m_variant_groups.push_back(VariantGroup(m_variant_groups.size(), variants));
    return m_variant_groups.back().getID();
}


