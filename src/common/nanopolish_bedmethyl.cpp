//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bedmethyl -- tools for working with
// bedmethyl files
//
#include <algorithm>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
#include "nanopolish_bedmethyl.h"

std::vector<BedMethylRecord> read_bedmethyl_from_file(const std::string& filename)
{
    std::vector<BedMethylRecord> out;
    std::ifstream infile(filename);
    std::string line;
    while(getline(infile, line)) {
        BedMethylRecord bmr(line);
        out.push_back(bmr);
    }
    return out;
}

std::vector<BedMethylRecord> read_bedmethyl_for_region(const std::string& filename,
                                                       const std::string& contig,
                                                       int region_start,
                                                       int region_end)
{
    std::vector<BedMethylRecord> out;
    std::ifstream infile(filename);
    std::string line;
    while(getline(infile, line)) {
        BedMethylRecord bmr(line);
        if(bmr.ref_name == contig &&
           bmr.ref_start_position >= region_start &&
           bmr.ref_end_position <= region_end)
        {
            out.push_back(bmr);
        }
    }
    return out;
}
