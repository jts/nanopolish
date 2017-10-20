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
#include "nanopolish_alphabet.h"
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
        
        alphabet = "";
        read_idx = -1;
        region_start = -1;
        region_end = -1;
    }

    // returns the pore model that should be used, based on the alphabet
    const PoreModel* get_model() const;

    // Mandatory
    SquiggleRead* sr;
    const faidx_t* fai;
    const bam_hdr_t* hdr;
    const bam1_t* record;
    size_t strand_idx;
    
    // optional
    std::string alphabet;
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
    std::string model_kmer;
    char hmm_state;
};

// Entry point from nanopolish.cpp
int eventalign_main(int argc, char** argv);

// print the alignment as a tab-separated table
void emit_event_alignment_tsv(FILE* fp,
                              const SquiggleRead& sr,
                              uint32_t strand_idx,
                              const EventAlignmentParameters& params,
                              const std::vector<EventAlignment>& alignments);

// The main function to realign a read
std::vector<EventAlignment> align_read_to_ref(const EventAlignmentParameters& params);

// get the specified reference region, threadsafe
std::string get_reference_region_ts(const faidx_t* fai, const char* ref_name, int start, int end, int* fetched_len);

#endif
