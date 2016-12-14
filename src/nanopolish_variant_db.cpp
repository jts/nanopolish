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


Combinations::Combinations(size_t N, size_t k, CombinationOption option) 
{ 
    assert(k <= N);
    m_option = option;

    if(m_option == CO_WITHOUT_REPLACEMENT) {
        initial_setmask.resize(N, false);
        std::fill(initial_setmask.begin(), initial_setmask.begin() + k, true);
    } else {
        // Uses the bitvector approach from https://www.mathsisfun.com/combinatorics/combinations-permutations.html
        initial_setmask.resize(N + k - 1, false);
        std::fill(initial_setmask.begin(), initial_setmask.begin() + k, true);
    }
    setmask = initial_setmask;
    m_done = false;
}

bool Combinations::done()
{
    return m_done;
}

std::vector<size_t> Combinations::get() const
{
    if(m_option == CO_WITHOUT_REPLACEMENT) {
        return _get_without_replacement();
    } else {
        return _get_with_replacement();
    }
}

std::vector<size_t> Combinations::_get_without_replacement() const
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

std::vector<size_t> Combinations::_get_with_replacement() const
{
    std::vector<size_t> out;
    size_t curr_id = 0;
    for(size_t i = 0; i < setmask.size(); ++i) {

        if(setmask[i]) {
            out.push_back(curr_id);
        } else {
            curr_id++;
        }
    }
    return out;
}

std::string Combinations::get_as_string() const
{
    std::vector<size_t> input = get();
    std::stringstream ss;
    for(size_t i = 0; i < input.size(); ++i) {
        ss << input[i];
        if(i < input.size() - 1) {
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


