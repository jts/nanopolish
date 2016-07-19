//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_haplotype - a haplotype derived from 
// a reference sequence and a set of variants
//
#include "nanopolish_haplotype.h"

// Definitions
const size_t Haplotype::INSERTED_POSITION = std::string::npos;

Haplotype::Haplotype(const std::string& ref_name,
                     const size_t ref_position,
                     const std::string& ref_sequence) : 
                        m_ref_name(ref_name),
                        m_ref_position(ref_position),
                        m_reference(ref_sequence)
{
    m_sequence = m_reference;
    m_coordinate_map.resize(m_reference.size());
    for(size_t i = 0; i < m_coordinate_map.size(); ++i) {
        m_coordinate_map[i] = m_ref_position + i;
    }
}
 
//       
bool Haplotype::apply_variant(const Variant& v)
{
    // Search the coordinate map for the reference position
    size_t derived_idx = _find_derived_index_by_ref_lower_bound(v.ref_position);

    // if we could not find the reference position in the map
    // this variant is incompatable with the haplotype, do nothing
    if(derived_idx == m_coordinate_map.size() || 
       m_coordinate_map[derived_idx] != v.ref_position) 
    {
        return false;
    }

    // Check that the string matches
    size_t rl = v.ref_seq.length();
    size_t al = v.alt_seq.length();

    // no match, variant conflicts with haplotype sequence
    if(m_sequence.substr(derived_idx, rl) != v.ref_seq) {
        return false;
    }

    // update sequence
    m_sequence.replace(derived_idx, rl, v.alt_seq);

    // update coordinate map
    
    // make a pair of iterators that bound the changed sequence
    std::vector<size_t>::iterator fi = m_coordinate_map.begin() + derived_idx;
    std::vector<size_t>::iterator li = fi + rl;
    
    // erase the positions of the changed bases
    std::vector<size_t>::iterator ii = m_coordinate_map.erase(fi, li);

    // insert new positions for the alt bases with invalid indices
    m_coordinate_map.insert(ii, al, INSERTED_POSITION);
    
    // sanity check
    assert(m_coordinate_map.size() == m_sequence.size());

    m_variants.push_back(v);
    return true;
}

// return a new haplotype subsetted by reference coordinates
Haplotype Haplotype::substr_by_reference(size_t start, size_t end)
{
    assert(start >= m_ref_position);
    assert(start <= m_ref_position + m_reference.length());
    
    assert(end >= m_ref_position);
    assert(end <= m_ref_position + m_reference.length());

    size_t derived_base_start = _find_derived_index_by_ref_lower_bound(start);
    size_t derived_base_end = _find_derived_index_by_ref_lower_bound(end);
    
    // Bump out the reference coordinate to encompass the complete range (start, end)
    while(m_coordinate_map[derived_base_start] > start ||
          m_coordinate_map[derived_base_start] == INSERTED_POSITION)
    { 
        derived_base_start -= 1;
    }

    assert(derived_base_start != m_coordinate_map.size());
    assert(derived_base_end != m_coordinate_map.size());
    assert(m_coordinate_map[derived_base_start] <= start);
    assert(m_coordinate_map[derived_base_end] >= end);

    start = m_coordinate_map[derived_base_start];
    end = m_coordinate_map[derived_base_end];
    
    Haplotype ret(m_ref_name,
                  start,
                  m_reference.substr(start - m_ref_position, end - start + 1));

    ret.m_sequence = m_sequence.substr(derived_base_start, derived_base_end - derived_base_start + 1);
    ret.m_coordinate_map = std::vector<size_t>(m_coordinate_map.begin() + derived_base_start,
                                               m_coordinate_map.begin() + derived_base_end + 1);

    assert(ret.m_coordinate_map.front() == start);
    assert(ret.m_coordinate_map.back() == end);
    assert(ret.m_coordinate_map.size() == ret.m_sequence.size());

    return ret;
}

size_t Haplotype::_find_derived_index_by_ref_lower_bound(size_t ref_index)
{
    for(size_t i = 0; i < m_coordinate_map.size(); ++i) {
        if(m_coordinate_map[i] != INSERTED_POSITION && m_coordinate_map[i] >= ref_index) {
            return i;
        }
    }
    return m_coordinate_map.size();
}

