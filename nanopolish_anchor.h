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

#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "nanopolish_fast5_map.h"

struct AlignedPair
{
    int ref_pos;
    int read_pos;
};

// An event index and orientation that gives us a handle
// into the event sequence for some SquiggleRead
struct HMMStrandAnchor
{
    //
    HMMStrandAnchor() : event_idx(-1), rc(false) {}
    HMMStrandAnchor(int ei, bool f) : event_idx(ei), rc(f) {}

    //
    int event_idx;
    bool rc; // with respect to consensus
};

// A pair of vectors containing all of the anchors
// for both strands of a SquiggleRead
struct HMMReadAnchorSet
{
    std::vector<HMMStrandAnchor> strand_anchors[NUM_STRANDS];
};

// This data structure represents a column of a 
// multiple alignment where the base sequence
// is a subsequence of a contig that we
// we have mapped events to. It also holds alternative
// sequences sampled from the reads at starting at this column.
struct HMMAnchoredColumn
{
    std::string base_sequence;
    std::vector<HMMStrandAnchor> anchors;
    std::vector<std::string> alt_sequences;
};

//
struct HMMRealignmentInput
{
    std::vector<SquiggleRead> reads;
    std::vector<HMMAnchoredColumn> anchored_columns;
};

// functions
HMMRealignmentInput build_input_for_region(const std::string& bam_filename, 
                                           const std::string& ref_filename, 
                                           const Fast5Map& read_name_map, 
                                           const std::string& contig_name,
                                           int start, 
                                           int end, 
                                           int stride);



// Return a vector specifying pairs of bases that have been aligned to each other
std::vector<AlignedPair> get_aligned_pairs(const bam1_t* record);

std::vector<int> uniformally_sample_read_positions(const std::vector<AlignedPair>& aligned_pairs,
                                                   int ref_start,
                                                   int ref_end,
                                                   int ref_stride);

#endif
