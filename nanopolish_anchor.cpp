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
#include "htslib/htslib/faidx.h"
#include "nanopolish_common.h"
#include "nanopolish_anchor.h"

void build_input_for_region(const std::string& bam_filename, 
                            const std::string& ref_filename, 
                            const Fast5Map& read_name_map, 
                            const std::string& contig_name,
                            int start, 
                            int end, 
                            int stride)
{
    // load bam file
    htsFile* bam_fh = sam_open(bam_filename.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = bam_filename + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    int contig_id = bam_name2id(hdr, contig_name.c_str());
    
    // load reference fai file
    faidx_t *fai = fai_load(ref_filename.c_str());

    // load the reference sequence for this region
    int fetched_len = 0;
    char* ref_segment = faidx_fetch_seq(fai, contig_name.c_str(), start, end, &fetched_len);

    // Initialize iteration    
    bam1_t* record = bam_init1();
    hts_itr_t* itr = sam_itr_queryi(bam_idx, contig_id, start, end);
    
    // Iterate over reads aligned here
    printf("Iter: %d %d %d\n", itr->tid, itr->beg, itr->end);
    
    std::vector<SquiggleRead> squiggle_reads;
    std::vector<HMMReadAnchorSet> read_anchors;
    
    // Load the SquiggleReads aligned to this region and the bases
    // that are mapped to our reference anchoring positions
    int result;
    while((result = sam_itr_next(bam_fh, itr, record)) >= 0) {

        // Load a squiggle read for the mapped read
        std::string read_name = bam_get_qname(record);
        std::string fast5_path = read_name_map.get_path(read_name);

        // load read
        squiggle_reads.push_back(SquiggleRead(read_name, fast5_path));

        // TODO: just read this from the fast5
        int read_len = record->core.l_qseq;
        std::string read_seq(read_len, '\0');
        uint8_t* packed_seq = bam_get_seq(record);
        
        // unpack
        for(int i = 0; i < read_len; ++i) {
            read_seq[i] = seq_nt16_str[bam_seqi(packed_seq, i)];
        }

        squiggle_reads.back().twod_sequence = read_seq;

        // parse alignments to reference
        read_anchors.push_back(build_anchors_for_read(record, start, end, stride));
    }

    // The HMMReadAnchorSet contains anchors for each strand of a read
    // laid out in a vector. Transpose this data so we have one anchor
    // for every read column-wise.
    assert(!read_anchors.empty());
    assert(read_anchors.front().strand_anchors[T_IDX].size() == 
           read_anchors.back().strand_anchors[T_IDX].size());

    // same for all anchors, so we use the first
    size_t num_anchors = read_anchors.front().strand_anchors[T_IDX].size(); 
    size_t num_strands = read_anchors.size() * 2;
    
    std::vector<HMMAnchoredColumn> columns(num_anchors);
    for(size_t ai = 0; ai < num_anchors; ++ai) {
        
        HMMAnchoredColumn& column = columns[ai];

        for(size_t rai = 0; rai < read_anchors.size(); ++rai) {
            HMMReadAnchorSet& ras = read_anchors[rai];
            assert(ai < ras.strand_anchors[T_IDX].size());
            assert(ai < ras.strand_anchors[C_IDX].size());

            column.anchors.push_back(ras.strand_anchors[T_IDX][ai]);
            column.anchors.push_back(ras.strand_anchors[C_IDX][ai]);
        }

        // Add sequences except for last anchor
        if(ai != num_anchors - 1) {
            
            // base
            column.base_sequence = std::string(ref_segment + ai * stride, stride);
            printf("Base: %s\n", column.base_sequence.c_str());

            // alts
            for(size_t rai = 0; rai < read_anchors.size(); ++rai) {
                HMMReadAnchorSet& ras = read_anchors[rai];
                int32_t b1 = ras.strand_anchors[T_IDX][ai].base_idx;
                int32_t b2 = ras.strand_anchors[T_IDX][ai + 1].base_idx;
                if(b1 >= 0 && b2 >= 0)  {
                    column.alt_sequences.push_back(squiggle_reads[rai].twod_sequence.substr(b1, b2 - b1 + 1));
                    printf("Alt[%zu]:  %s\n", rai, column.alt_sequences.back().c_str());
                }
            }
        }

    }

    for(size_t rai = 0; rai < read_anchors.size(); ++rai) {
        printf("T\t");
        for(size_t sai = 0; sai < read_anchors[rai].strand_anchors[T_IDX].size(); ++sai) {
            printf("%d\t", read_anchors[rai].strand_anchors[T_IDX][sai].event_idx);
        }
        printf("\n");
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    bam_destroy1(record);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    free(ref_segment);
}

HMMReadAnchorSet build_anchors_for_read(bam1_t* record, int start, int end, int stride)
{
    //
    printf("Record start: %d end: %d name: %s\n", record->core.pos, bam_endpos(record), (char*)record->data);

    // We want an anchor for every stride-th base, even if this read is not aligned there
    // The missing anchors will be flagged as -1
    uint32_t num_anchors = ((end - start) / stride) + 1;

    HMMReadAnchorSet out;
    out.strand_anchors[T_IDX].resize(num_anchors);
    out.strand_anchors[C_IDX].resize(num_anchors);

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
        } else if(cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if(cigar_op == BAM_CHARD_CLIP) {
            // no increment   
        } else {
            printf("Cigar: %d\n", cigar_op);
            assert(false && "Unhandled cigar operation");
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(ref_pos >= start && ref_pos <= end && ref_pos % stride == 0 && ref_inc > 0) {

                //printf("Match %d %d\n", ref_pos, read_pos);
                uint32_t anchor_id = (ref_pos - start) / stride;
                assert(anchor_id < num_anchors);

                // Process both strands of the SquiggleRead here and generate strand anchors
                // mapping from a reference base to a nanopore event
                // FAKE DATA
                out.strand_anchors[T_IDX][anchor_id] = HMMStrandAnchor(read_pos, read_pos, false);
                out.strand_anchors[C_IDX][anchor_id] = HMMStrandAnchor(10000 - read_pos, 10000 - read_pos, true);
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
    return out;
}

