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

#include "nanopolish_variant.h"

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
    
        // get the reference location
        const std::string get_reference_name() const { return m_ref_name; }
        const size_t get_reference_position() const { return m_ref_position; }
        const size_t get_reference_end() const { return m_ref_position + m_reference.length(); }

        // return the reference position corresponding to base i of the haplotype
        // returns std::string::npos if the base was inserted into the haplotype
        // and therefore has no corresponding reference base
        size_t get_reference_position_for_haplotype_base(size_t i) const;

        // add a variant into the haplotype
        // returns true if the variant is successfully added to the haplotype
        bool apply_variant(const Variant& v);

        // return all the variants on this haplotype
        std::vector<Variant> get_variants() const { return m_variants; }
        
        // Set ref_lower and ref_upper to be valid reference (ie non-deleted/inserted) positions that
        // contain the lower/supper positions on the haplotype.
        // If no such range exists then the return value(s) is set to std::string::npos

        // Extend the haplotype range (hap_lower/hap_upper) until both have a cooresponding
        // reference base. Set ref_lower/ref_upper to these values. If a range cannot be found,
        // the out parameters are set to std::string::npos
        void get_enclosing_reference_range_for_haplotype_range(size_t& hap_lower, size_t& hap_upper,
                                                               size_t& ref_lower, size_t& ref_upper) const;

        // return a new haplotype subsetted by reference coordinates
        Haplotype substr_by_reference(size_t start, size_t end) const;

    private:
        
        // functions
        Haplotype(); // not allowed
        
        // Find the first derived index that has a corresponding
        // reference position which is not less than ref_index.
        // This mimics std::lower_bound
        size_t _find_derived_index_by_ref_lower_bound(size_t ref_index) const;

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
