//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bedmethyl -- tools for working with
// bedmethyl files
//
#ifndef NANOPOLISH_BEDMETHYL_H
#define NANOPOLISH_BEDMETHYL_H

#include <sstream>
#include "stdaln.h"
#include "nanopolish_common.h"

struct BedMethylRecord
{
    BedMethylRecord() { }
    BedMethylRecord(const std::string& line) { read_bed(line); }

    void read_bed(const std::string& line)
    {
        std::stringstream ss(line);
        std::string dummy;
        ss >> ref_name;
        ss >> ref_start_position;
        ss >> ref_end_position;
        ss >> dummy; // name, not used
        ss >> dummy; // score, not used
        ss >> strand;
        ss >> dummy; // display item, not used
        ss >> dummy; // display item, not used
        ss >> dummy; // display item, not used
        ss >> coverage;
        ss >> percent_methylated;
    }

    // generate a unique identifier for this variant
    std::string key() const
    {
        std::stringstream out;
        out << ref_name << ':' << ref_start_position << ':' << strand;
        return out.str();
    }

    std::string ref_name;
    size_t ref_start_position;
    size_t ref_end_position;
    std::string strand;
    int coverage;
    float percent_methylated;
};

inline bool sortBedMethylByPosition(const BedMethylRecord& a, const BedMethylRecord& b) 
{ 
    return a.ref_name == b.ref_name ? 
        a.ref_start_position < b.ref_start_position : 
        a.ref_name < b.ref_name; 
}

class BedMethylKeyComp
{
    public: 
        inline bool operator()(const BedMethylRecord& a, const BedMethylRecord& b)
        {
            return a.key() < b.key();
        }
};

class BedMethylKeyEqualityComp
{
    public: 
        inline bool operator()(const BedMethylRecord& a, const BedMethylRecord& b)
        {
            return a.key() == b.key();
        }
};

// Read a collection of variants from a VCF file
std::vector<BedMethylRecord> read_bedmethyl_from_file(const std::string& filename);
std::vector<BedMethylRecord> read_bedmethyl_for_region(const std::string& filename,
                                                       const std::string& contig,
                                                       int region_start,
                                                       int region_end);

#endif
