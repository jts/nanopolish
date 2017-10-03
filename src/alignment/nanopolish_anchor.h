//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_anchor - a collection of data types
// for representing a set of event-to-sequence
// mappings.
#ifndef NANOPOLISH_ANCHOR_H
#define NANOPOLISH_ANCHOR_H

#include <memory>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

struct AlignedPair
{
    int ref_pos;
    int read_pos;
};

struct AlignedPairRefLBComp
{
    bool operator()(const AlignedPair& o, int v) { return o.ref_pos < v; }
};

struct AlignedPairRefUBComp
{
    bool operator()(int v, const AlignedPair& o) { return v < o.ref_pos; }
};

// typedefs
typedef std::vector<AlignedPair>::iterator AlignedPairIter;
typedef std::vector<AlignedPair>::const_iterator AlignedPairConstIter;

// Return a vector specifying read-to-reference alignment segments
// Each segment contains a vector pairs of bases that have been aligned to each other.
// Each segment is defined a sequence of cigar match/ins/deletion options, followed by
// a reference skip (N) operation.
// This function can handle an "event cigar" bam record, which requires the ability
// for event indices to be in ascending or descending order. In the latter case
// read_stride should be -1
typedef std::vector<AlignedPair> AlignedSegment;
std::vector<AlignedSegment> get_aligned_segments(const bam1_t* record, int read_stride = 1);

#endif
