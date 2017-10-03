//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_anchor - a collection of data types
// for representing a set of event-to-sequence
// mappings.
#include <vector>
#include <string>
#include <stdio.h>
#include <assert.h>
#include "htslib/faidx.h"
#include "nanopolish_common.h"
#include "nanopolish_anchor.h"
#include "nanopolish_scorereads.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_squiggle_read.h"

std::vector<AlignedSegment> get_aligned_segments(const bam1_t* record, int read_stride)
{
    std::vector<AlignedSegment> out;
    // Initialize first segment
    out.push_back(AlignedSegment());

    // This code is derived from bam_fillmd1_core
    //uint8_t *ref = NULL;
    //uint8_t *seq = bam_get_seq(record);
    uint32_t *cigar = bam_get_cigar(record);
    const bam1_core_t *c = &record->core;

    // read pos is an index into the original sequence that is present in the FASTQ
    // on the strand matching the reference
    int read_pos = 0;

    // query pos is an index in the query string that is recorded in the bam
    // we record this as a sanity check
    //int query_pos = 0;
    
    int ref_pos = c->pos;

    for (int ci = 0; ci < c->n_cigar; ++ci) {
        
        int cigar_len = cigar[ci] >> 4;
        int cigar_op = cigar[ci] & 0xf;

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;
 
        // Process match between the read and the reference
        bool is_aligned = false;
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = true;
            read_inc = read_stride;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL) {
            ref_inc = 1;   
        } else if(cigar_op == BAM_CREF_SKIP) {
            // end the current segment and start a new one
            out.push_back(AlignedSegment());
            ref_inc = 1;
        } else if(cigar_op == BAM_CINS) {
            read_inc = read_stride;
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1; // special case, do not use read_stride
        } else if(cigar_op == BAM_CHARD_CLIP) {
            read_inc = 0;
        } else {
            printf("Cigar: %d\n", cigar_op);
            assert(false && "Unhandled cigar operation");
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(is_aligned) {
                out.back().push_back({ref_pos, read_pos});
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
    return out;
}

std::vector<int> uniformally_sample_read_positions(const std::vector<AlignedPair>& aligned_pairs,
                                                   int ref_start,
                                                   int ref_end,
                                                   int ref_stride)
{
    uint32_t num_anchors = ((ref_end - ref_start) / ref_stride) + 1;
    std::vector<int> out(num_anchors, -1);

    for(size_t pair_idx = 0; pair_idx < aligned_pairs.size(); ++pair_idx) {

        // We use a loop here in case there is no read base
        // aligned to one of the anchor positions we are interested in.
        // In this situation the loop will catch it and emit the last seen
        // read position
        int ref_pos = aligned_pairs[pair_idx].ref_pos;
        int end_pos = pair_idx + 1 != aligned_pairs.size() ? aligned_pairs[pair_idx + 1].ref_pos : ref_pos + 1;

        for(; ref_pos < end_pos; ++ref_pos) {
            if(ref_pos >= ref_start && ref_pos <= ref_end && ref_pos % ref_stride == 0) {
                uint32_t anchor_id = (ref_pos - ref_start) / ref_stride;
                assert(anchor_id < num_anchors);
                out[anchor_id] = aligned_pairs[pair_idx].read_pos;
            }
        }
    }
    return out;
}
