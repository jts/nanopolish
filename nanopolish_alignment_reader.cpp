//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_reader -- parse read alignments
// from a bam file
//
#include <stdio.h>
#include <assert.h>
#include "nanopolish_alignment_reader.h"
#include "htslib/htslib/sam.h"

void sample_bam(const std::string& filename, int ref_id, int start, int end, int stride)
{
    // load bam file
    htsFile* bam_fh = sam_open(filename.c_str(), "r");
    assert(bam_fh != NULL);

    // load index file
    std::string index_filename = filename + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // Initialize iteration    
    bam1_t* record = bam_init1();
    hts_itr_t* itr = sam_itr_queryi(bam_idx, ref_id, start, end);
    
    // Iterate over reads aligned here
    printf("Iter: %d %d %d\n", itr->tid, itr->beg, itr->end);
    
    int result;
    while((result = sam_itr_next(bam_fh, itr, record)) >= 0) {
        sample_read(record, start, end, stride);
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_destroy1(record);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
}

void sample_read(bam1_t* record, int start, int end, int stride)
{
    int endpos = bam_endpos(record);
    printf("Record start: %d end: %d name: %s\n", record->core.pos, endpos, (char*)record->data);

    // This code is derived from bam_fillmd1_core
    uint8_t *ref = NULL;
    uint8_t *seq = bam_get_seq(record);
    uint32_t *cigar = bam_get_cigar(record);
    bam1_core_t *c = &record->core;
    kstring_t *str;

    int read_pos = 0;
    int ref_pos = c->pos;

    for (int ci = 0; ci < c->n_cigar && ref_pos <= end; ++ci) {
        
        int cigar_len = cigar[ci] >> 4;
        int cigar_op = cigar[ci] & 0xf;

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;
        
        // Process match between the read and the reference
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            read_inc = 1;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL || cigar_op == BAM_CREF_SKIP) {
            ref_inc = 1;   
        } else if (cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else {
            assert(false && "Unhandled cigar operation");
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(ref_pos >= start && ref_pos <= end && ref_pos % stride == 0) {
                printf("Match %d %d\n", ref_pos, read_pos);
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }    
}

