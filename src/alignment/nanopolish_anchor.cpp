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

HMMRealignmentInput build_input_for_region(const std::string& bam_filename,
                                           const std::string& ref_filename,
                                           const Fast5Map& read_name_map,
                                           const std::string& contig_name,
                                           int start,
                                           int end,
                                           int stride)
{
    // Initialize return data
    HMMRealignmentInput ret;

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

    // Adjust end position to make sure we don't go out-of-range
    end = std::min(end, faidx_seq_len(fai, contig_name.c_str()));

    // load the reference sequence for this region
    int fetched_len = 0;
    char* ref_segment = faidx_fetch_seq(fai, contig_name.c_str(), start, end, &fetched_len);
    if(ref_segment == NULL) {
        exit(EXIT_FAILURE);
    }
    ret.original_sequence = ref_segment;

    // Initialize iteration
    bam1_t* record = bam_init1();
    hts_itr_t* itr = sam_itr_queryi(bam_idx, contig_id, start, end);
   
    // Iterate over reads aligned here
    std::vector<HMMReadAnchorSet> read_anchors;
    std::vector<std::vector<std::string>> read_substrings;

    // kmer size of pore model
    uint32_t k = 0;

    // Load the SquiggleReads aligned to this region and the bases
    // that are mapped to our reference anchoring positions
    int result;
    while((result = sam_itr_next(bam_fh, itr, record)) >= 0) {

        // Load a squiggle read for the mapped read
        std::string read_name = bam_get_qname(record);
        std::string fast5_path = read_name_map.get_path(read_name);

        // load read
        ret.reads.push_back(std::unique_ptr<SquiggleRead>(new SquiggleRead(read_name, fast5_path)));
        SquiggleRead& sr = *ret.reads.back();

        // Recalibrate each strand
        for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; strand_idx++) {
            if(!sr.has_events_for_strand(strand_idx)) {
                continue;
            }

            std::vector<EventAlignment> ao = alignment_from_read(sr, strand_idx, -1,
                                                                 "", fai, hdr,
                                                                 record, -1, -1);
            recalibrate_model(sr, strand_idx, ao, &gDNAAlphabet, true, true);
        }

        k = sr.pore_model[T_IDX].k;

        // parse alignments to reference
        std::vector<AlignedPair> aligned_pairs = get_aligned_pairs(record);
        std::vector<int> read_bases_for_anchors = 
            uniformally_sample_read_positions(aligned_pairs, start, end, stride);

        /*
        for(size_t i = 0; i < read_bases_for_anchors.size(); ++i) {
            printf("%d ", read_bases_for_anchors[i]);
        }
        printf("\n");
        */

        // Convert the read base positions into event indices for both strands
        HMMReadAnchorSet event_anchors;
        event_anchors.strand_anchors[T_IDX].resize(read_bases_for_anchors.size());
        event_anchors.strand_anchors[C_IDX].resize(read_bases_for_anchors.size());

        bool do_base_rc = bam_is_rev(record);
        bool template_rc = do_base_rc;
        bool complement_rc = !do_base_rc;

        for(size_t ai = 0; ai < read_bases_for_anchors.size(); ++ai) {

            int read_kidx = read_bases_for_anchors[ai];

            // read not aligned to this reference position
            if(read_kidx == -1) {
                continue;
            }

            if(do_base_rc)
                read_kidx = sr.flip_k_strand(read_kidx);

            // If the aligned base is beyong the start of the last k-mer of the read, skip
            if(read_kidx < 0 || read_kidx >= (int)(sr.read_sequence.size() - k + 1)) {
                continue;
            }

            int template_idx = sr.get_closest_event_to(read_kidx, T_IDX);
            int complement_idx = sr.has_events_for_strand(C_IDX) ? sr.get_closest_event_to(read_kidx, C_IDX) : -1;

            assert(template_idx != -1 && (!sr.has_events_for_strand(C_IDX) || complement_idx != -1));
            assert(template_idx < (int)sr.events[T_IDX].size());
            assert(sr.events[C_IDX].empty() || complement_idx < (int)sr.events[C_IDX].size());

            event_anchors.strand_anchors[T_IDX][ai] = { template_idx, template_rc };
            event_anchors.strand_anchors[C_IDX][ai] = { complement_idx, complement_rc };

            // If this is not the last anchor, extract the sequence of the read
            // from this anchor to the next anchor as an alternative assembly
            if(ai < read_bases_for_anchors.size() - 1) {
                int start_kidx = read_bases_for_anchors[ai];
                int end_kidx = read_bases_for_anchors[ai + 1];
                int max_kidx = sr.read_sequence.size() - k;

                // flip
                if(do_base_rc) {
                    start_kidx = sr.flip_k_strand(start_kidx);
                    end_kidx = sr.flip_k_strand(end_kidx);

                    // swap
                    int tmp = end_kidx;
                    end_kidx = start_kidx;
                    start_kidx = tmp;
                }

                // clamp values within range
                start_kidx = start_kidx >= 0 ? start_kidx : 0;
                end_kidx = end_kidx <= max_kidx ? end_kidx : max_kidx;
                
                std::string s = sr.read_sequence.substr(start_kidx, end_kidx - start_kidx + k);

                if(do_base_rc) {
                    s = gDNAAlphabet.reverse_complement(s);
                }

                if(ai >= read_substrings.size())
                    read_substrings.resize(ai + 1);

                read_substrings[ai].push_back(s);
            }
        }

        read_anchors.push_back(event_anchors);
    }

    // If there are no reads aligned to this segment return empty data
    if(read_anchors.empty()) {
        return ret;
    }

    // The HMMReadAnchorSet contains anchors for each strand of a read
    // laid out in a vector. Transpose this data so we have one anchor
    // for every read column-wise.
    assert(read_anchors.front().strand_anchors[T_IDX].size() == 
           read_anchors.back().strand_anchors[T_IDX].size());

    // same for all anchors, so we use the first
    size_t num_anchors = read_anchors.front().strand_anchors[T_IDX].size(); 
    size_t num_strands = read_anchors.size() * 2;
    
    ret.anchored_columns.resize(num_anchors);
    for(size_t ai = 0; ai < num_anchors; ++ai) {
        
        HMMAnchoredColumn& column = ret.anchored_columns[ai];

        for(size_t rai = 0; rai < read_anchors.size(); ++rai) {
            HMMReadAnchorSet& ras = read_anchors[rai];
            assert(ai < ras.strand_anchors[T_IDX].size());
            assert(ai < ras.strand_anchors[C_IDX].size());

            column.anchors.push_back(ras.strand_anchors[T_IDX][ai]);
            column.anchors.push_back(ras.strand_anchors[C_IDX][ai]);
        }
        assert(column.anchors.size() == num_strands);

        // Add sequences except for last anchor
        if(ai != num_anchors - 1) {
            
            // base, these sequences need to overlap by k - 1 bases
            int base_length = stride + k;
            if((int)ai * stride + base_length > fetched_len)
                base_length = fetched_len - ai * stride;

            column.base_sequence = std::string(ref_segment + ai * stride, base_length);
            column.base_contig = contig_name;
            column.base_start_position = start + ai * stride;
            assert(column.base_sequence.back() != '\0');

            // alts
            if(ai < read_substrings.size())
                column.alt_sequences = read_substrings[ai];
        }
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    bam_destroy1(record);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    free(ref_segment);

    return ret;
}

std::vector<AlignedPair> get_aligned_pairs(const bam1_t* record, int read_stride)
{
    std::vector<AlignedPair> out;

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
        } else if(cigar_op == BAM_CDEL || cigar_op == BAM_CREF_SKIP) {
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
                out.push_back({ref_pos, read_pos});
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
