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

        // constructor
        Haplotype(const std::string& ref_name,
                  const size_t ref_position,
                  const std::string& ref_sequence);
        
        // get the sequence of the haplotype
        const std::string& get_sequence() const { return m_sequence; } 
        
        // get the sequence of the reference
        const std::string& get_reference() const { return m_reference; } 
    
        // add a variant into the haplotype
        void apply_variant(const Variant& v);

        // return a new haplotype subsetted by reference coordinates
        Haplotype substr_by_reference(size_t start, size_t end);

    private:
        
        // functions
        Haplotype(); // not allowed
        
        // Find the first derived index that has a corresponding
        // reference position which is not less than ref_index.
        // This mimics std::lower_bound
        size_t _find_derived_index_by_ref_lower_bound(size_t ref_index);

        //
        // data
        //

        // the name of the reference contig/chromosome this haplotype is from
        std::string m_ref_name;

        // the start position of the reference sequence on the ref contig/chromosome
        size_t m_ref_position;

        // the original sequence this haplotype is based on
        std::string m_reference;
        
        // the sequence of the haplotype
        std::string m_sequence;

        // the set of variants this haplotype contains
        std::vector<Variant> m_variants;

        // a mapping from bases of the derived sequence
        // to their original reference position
        std::vector<size_t> m_coordinate_map;

        // a constant value indicating inserted sequence in the coordinate map
        static const size_t INSERTED_POSITION;
};

#endif
