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

void parse_bam(const std::string& filename);

#endif
