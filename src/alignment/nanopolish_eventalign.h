//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_eventalign.cpp -- align squiggle events
// to kmers of a sequence
//
#ifndef NANOPOLISH_EVENTALIGN_H
#define NANOPOLISH_EVENTALIGN_H

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "nanopolish_common.h"

//
// Structs
//
struct EventAlignmentParameters
{
    EventAlignmentParameters()
    {
        sr = NULL;
        fai = NULL;
        hdr = NULL;
        record = NULL;
        strand_idx = NUM_STRANDS;
        
        read_idx = -1;
        region_start = -1;
        region_end = -1;
    }

    // Mandatory
    SquiggleRead* sr;
    const faidx_t* fai;
    const bam_hdr_t* hdr;
    const bam1_t* record;
    size_t strand_idx;
    
    // optional
    int read_idx;
    int region_start;
    int region_end;
};

struct EventAlignment
{
    // ref data
    std::string ref_name;
    std::string ref_kmer;
    int ref_position;

    // event data
    size_t read_idx;
    int strand_idx;
    int event_idx;
    bool rc;

    // hmm data
    char hmm_state;
};

// Entry point from nanopolish.cpp
int eventalign_main(int argc, char** argv);

// The main function to realign a read
std::vector<EventAlignment> align_read_to_ref(const EventAlignmentParameters& params);

#endif
