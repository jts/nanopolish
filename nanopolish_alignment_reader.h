//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_reader -- parse read alignments
// from a bam file
//
#ifndef NANOPOLISH_ALIGNMENT_READER_H
#define NANOPOLISH_ALIGNMENT_READER_H

#include <string>
#include "htslib/htslib/sam.h"

void sample_bam(const std::string& filename, int ref_id, int start, int end, int stride);
void sample_read(bam1_t* record, int start, int end, int stride);

#endif
