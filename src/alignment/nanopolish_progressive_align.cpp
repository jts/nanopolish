//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_progressive_align -- align squiggle events
// to kmers of a sequence using a banded, 
// progressive algorithm
//
#include "nanopolish_progressive_align.h"
#include "nanopolish_profile_hmm.h"

void progressive_align(SquiggleRead& read,
                       const std::string& sequence)
{
    size_t strand_idx = 0;
    size_t EVENTS_PER_BLOCK = 200;
    size_t BASES_PER_BLOCK = 100;

    // Build a debug event->kmer map
    std::vector<size_t> kmer_for_event(read.events[strand_idx].size());
    for(size_t ki = 0; ki < read.base_to_event_map.size(); ++ki) {
        IndexPair& elem = read.base_to_event_map[ki].indices[0];
        if(elem.start == -1) {
            continue;
        }

        for(size_t j = elem.start; j <= elem.stop; ++j) {
            kmer_for_event[j] = ki;
        }
    }

    // For now we assume the first event in the array matches the first k-mer in the sequence
    size_t curr_k_idx = 0;
    size_t curr_event_idx = 0;

    fprintf(stderr, "aligning events for read %s [shift=%.2lf, scale=%.2lf]\n", read.read_name.substr(0, 6).c_str(),
                                                                                read.pore_model[0].shift,
                                                                                read.pore_model[0].scale);

    /*
    double events_per_base = (double)read.events[strand_idx].size() / sequence.size();
    double events_per_base_per_block_upper_bound = events_per_base * 1.5;
    fprintf(stderr, "parameters -- events_per_base: %.2lf upper bound: %.2lf\n", events_per_base, events_per_base_per_block_upper_bound);
    */

    while(1) {

        size_t end_k_idx = curr_k_idx + BASES_PER_BLOCK;
        size_t end_event_idx = curr_event_idx + EVENTS_PER_BLOCK;
        
        std::string block_seq = sequence.substr(curr_k_idx, end_k_idx - curr_k_idx);
        HMMInputSequence hmm_sequence(block_seq);
        
        HMMInputData input;
        input.read = &read;
        input.event_start_idx = curr_event_idx;
        input.event_stop_idx = end_event_idx;
        
        input.strand = strand_idx;
        input.event_stride = 1;
        input.rc = false;

        std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(hmm_sequence, input, HAF_ALLOW_POST_CLIP);
        
        for(size_t eai = 0; eai < event_alignment.size(); ++eai) {

            HMMAlignmentState& as = event_alignment[eai];
            /*
            if(as.state == 'K') {
                continue;
            }
            */
            size_t k_idx = curr_k_idx + as.kmer_idx;
            fprintf(stderr, "Event %zu aligns to kmer %zu [debug: %zu, state: %c]\n", as.event_idx, k_idx, kmer_for_event[as.event_idx], as.state);
            /*
            EventAlignment ea;
            
            // ref
            ea.ref_name = ref_name;
            ea.ref_position = curr_start_ref + as.kmer_idx;
            ea.ref_kmer = ref_seq.substr(ea.ref_position - ref_offset, k);

            // event
            ea.read_idx = params.read_idx;
            ea.strand_idx = params.strand_idx;
            ea.event_idx = as.event_idx;
            ea.rc = input.rc;

            // hmm
            ea.hmm_state = as.state;

            if(ea.hmm_state != 'B') {
                ea.model_kmer = hmm_sequence.get_kmer(as.kmer_idx, k, input.rc);
            } else {
                ea.model_kmer = std::string(k, 'N');
            }

            // store
            alignment_output.push_back(ea);

            // update
            last_event_output = as.event_idx;
            last_ref_kmer_output = curr_start_ref + as.kmer_idx;

            num_output += 1;
            */
        }

        break;
    }

    exit(0);
}
