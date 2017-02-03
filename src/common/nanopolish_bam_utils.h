//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bam_processor -- framework for iterating
// over a bam file and running an arbitrary function
// on each aligned read in parallel
//
#ifndef NANOPOLISH_BAM_UTILS_H
#define NANOPOLISH_BAM_UTILS_H

#include <string>
#include <vector>
#include "htslib/hts.h"
#include "htslib/sam.h"

// Allocate space for the variable-length fields
// in the bam record, and write them. If aux
// tags will be written to the record later, aux_reserve
// can be used to avoid re-allocation.
void write_bam_vardata(bam1_t* record,
                       const std::string& qname,
                       const std::vector<uint32_t> cigar,
                       const std::string& seq,
                       const std::string& qual,
                       size_t aux_reserve = 0);

#endif
