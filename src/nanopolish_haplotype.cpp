//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_haplotype - a haplotype derived from 
// a reference sequence and a set of variants
//
#include "nanopolish_haplotype.h"

Haplotype::Haplotype(const std::string& reference) : m_reference(reference)
{
    m_sequence = reference;
    m_coordinate_map.resize(m_reference.size());
    for(size_t i = 0; i < m_coordinate_map.size(); ++i) {
        m_coordinate_map[i] = i;
    }
}
 
//       
void Haplotype::apply_variant(const Variant& v)
{
    // Search the coordinate map for the reference position
    size_t derived_idx = 0;
    while(derived_idx < m_coordinate_map.size() && 
          m_coordinate_map[derived_idx] != v.ref_position) {
        derived_idx++;
    }

    // if we could not find the reference position in the map
    // this variant is incompatable with the haplotype, do nothing
    if(derived_idx == m_coordinate_map.size()) {
        return;
    }

    // Check that the string matches
    size_t rl = v.ref_seq.length();
    size_t al = v.alt_seq.length();

    // no match, variant conflicts with haplotype sequence
    if(m_sequence.substr(derived_idx, rl) != v.ref_seq) {
        return;
    }

    // apply variant
    m_sequence.replace(derived_idx, rl, v.alt_seq);

    // update coordinate map
    
    // make a pair of iterators that bound the changed sequence
    std::vector<size_t>::iterator fi = m_coordinate_map.begin() + derived_idx;
    std::vector<size_t>::iterator li = fi + rl;
    
    // erase the positions of the changed bases
    std::vector<size_t>::iterator ii = m_coordinate_map.erase(fi, li);

    // insert new positions for the alt bases with invalid indices
    m_coordinate_map.insert(ii, al, std::string::npos);
    
    // sanity check
    assert(m_coordinate_map.size() == m_sequence.size());

    m_variants.push_back(v);
}
