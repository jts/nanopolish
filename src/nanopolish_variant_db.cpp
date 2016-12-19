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

// http://stackoverflow.com/a/17050528
SizeTVecVec cartesian_product(const SizeTVecVec& input)
{
    SizeTVecVec s = {{}};
    for (const auto& u : input) {
        SizeTVecVec r;
        for (const auto& x : s) {
            for (const auto y : u) {
                r.push_back(x);
                r.back().push_back(y);
            }
        }
        s = move(r);
    }
    return s;
}

Combinations::Combinations(size_t N, size_t k, CombinationOption option) 
{ 
    assert(k <= N);
    m_rank = 0;
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
    m_rank++;
    m_done = (setmask == initial_setmask);
}

//
// VariantCombinations
//
VariantCombination::VariantCombination(const std::vector<size_t>& variant_ids) : m_variant_ids(variant_ids)
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

size_t VariantGroup::add_combination(const VariantCombination& vc)
{
    m_combinations.push_back(vc);
    m_scores.resize(m_combinations.size());
    return m_combinations.size() - 1;
}

std::string VariantGroup::get_vc_allele_string(size_t idx) const
{
    const VariantCombination& vc = m_combinations[idx];
    std::vector<bool> varmask(m_variants.size(), false);
    for(size_t i = 0; i < vc.get_num_variants(); ++i) {
        varmask[vc.get_variant_id(i)] = true;
    }

    std::stringstream ss;
    for(size_t i = 0; i < m_variants.size(); ++i) {
        ss << (varmask[i] ? m_variants[i].ref_seq : m_variants[i].alt_seq) << ",";
    }
    std::string out = ss.str();
    return out.substr(0, out.size() - 1);
}

void VariantGroup::set_combination_read_score(size_t combination_idx, const std::string& read_id, double score)
{
    assert(combination_idx < m_scores.size());
    m_scores[combination_idx][read_id] = score;

    auto itr = m_read_score_sum.find(read_id);

    if(itr == m_read_score_sum.end()) {
        m_read_score_sum[read_id] = score;
    } else {
        itr->second = add_logs(itr->second, score);
    }
}

double VariantGroup::get_combination_read_score(size_t combination_idx, const std::string& read_id) const
{
    assert(combination_idx < m_scores.size());
    const auto& itr = m_scores[combination_idx].find(read_id);
    assert(itr != m_scores[combination_idx].end());
    return itr->second;
}

std::vector<std::pair<std::string, double>> VariantGroup::get_read_sum_scores() const
{
    std::vector<std::pair<std::string, double>> out;
    for(const auto& itr : m_read_score_sum) {
        out.push_back(std::make_pair(itr.first, itr.second));
    }
    return out;
}

//
// VariantDB
//
size_t VariantDB::add_new_group(const std::vector<Variant>& variants)
{
    m_variant_groups.push_back(VariantGroup(m_variant_groups.size(), variants));
    return m_variant_groups.back().getID();
}
