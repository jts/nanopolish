//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_haplotype - a haplotype derived from 
// a reference sequence and a set of variants
//
#ifndef NANOPOLISH_HAPLOTYPE_H
#define NANOPOLISH_HAPLOTYPE_H

#include "nanopolish_variants.h"

class Haplotype
{
    public:
        Haplotype(const std::string& reference);

        const std::string& get_sequence() const { return m_sequence; } 
    
        void apply_variant(const Variant& v);

    private:

        // the original sequence this haplotype
        // is based on
        std::string m_reference;
        
        // the sequence of the haplotype
        std::string m_sequence;

        // the set of variants this haplotype contains
        std::vector<Variant> m_variants;

        // a mapping from bases of the derived sequence
        // to their original reference position
        std::vector<size_t> m_coordinate_map;
};

#endif
