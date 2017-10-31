//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
// for R7 data
//
#include <algorithm>
#include "nanopolish_profile_hmm_r7.h"

//#define DEBUG_FILL
//#define PRINT_TRAINING_MESSAGES 1

void profile_hmm_forward_initialize_r7(FloatMatrix& fm)
{
    // initialize forward calculation
    for(uint32_t si = 0; si < fm.n_cols; si++) {
        set(fm, 0, si, -INFINITY);
    }

    for(uint32_t ri = 0; ri < fm.n_rows; ri++) {
        set(fm, ri, PSR7_KMER_SKIP, -INFINITY);
        set(fm, ri, PSR7_EVENT_SPLIT, -INFINITY);
        set(fm, ri, PSR7_MATCH, -INFINITY);
    }
}

// Terminate the forward algorithm by calculating
// the probability of transitioning to the end state
// for all columns and a given row
float profile_hmm_forward_terminate_r7(const FloatMatrix& fm,
                                       const FloatMatrix& tm,
                                       uint32_t row)
{
    assert(false);
    return -INFINITY;
}

float profile_hmm_score_r7(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags)
{
    const uint32_t k = data.pore_model->k;
    uint32_t n_kmers = sequence.length() - k + 1;

    uint32_t n_states = PSR7_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

    uint32_t e_start = data.event_start_idx;
    uint32_t e_end = data.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate a matrix to hold the HMM result
    FloatMatrix fm;
    allocate_matrix(fm, n_rows, n_states);

    profile_hmm_forward_initialize_r7(fm);

    ProfileHMMForwardOutputR7 output(&fm);

    float score = profile_hmm_fill_generic_r7(sequence, data, e_start, flags, output);

    // cleanup
    free_matrix(fm);
    return score;
}

void profile_hmm_viterbi_initialize_r7(FloatMatrix& m)
{
    // Same as forward initialization
    profile_hmm_forward_initialize_r7(m);
}

std::vector<HMMAlignmentState> profile_hmm_align_r7(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags)
{
    std::vector<HMMAlignmentState> alignment;
    const uint32_t k = data.pore_model->k;

    uint32_t n_kmers = sequence.length() - k + 1;
    uint32_t n_states = PSR7_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

    uint32_t e_start = data.event_start_idx;
    uint32_t e_end = data.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;
    assert(n_events >= 2);

    uint32_t n_rows = n_events + 1;
    
    // Allocate matrices to hold the HMM result
    FloatMatrix vm;
    allocate_matrix(vm, n_rows, n_states);
    
    UInt8Matrix bm;
    allocate_matrix(bm, n_rows, n_states);

    ProfileHMMViterbiOutputR7 output(&vm, &bm);

    profile_hmm_viterbi_initialize_r7(vm);
    profile_hmm_fill_generic_r7(sequence, data, e_start, flags, output);

    // Traverse the backtrack matrix to compute the results
    int traversal_stride = data.event_stride;

#if HMM_REVERSE_FIX
    // Hack to support the fixed HMM
    // TODO: clean up
    traversal_stride = 1;
    if(data.event_stride == -1) {
        e_start = data.event_stop_idx;
    }
#endif
    
    // start from the last event matched to the last kmer
    uint32_t row = n_rows - 1;
    uint32_t col = PSR7_NUM_STATES * n_kmers + PSR7_MATCH;

    while(row > 0) {
        
        uint32_t event_idx = e_start + (row - 1) * traversal_stride;
        uint32_t block = col / PSR7_NUM_STATES;
        assert(block > 0);
        assert(get(vm, row, col) != -INFINITY);

        uint32_t kmer_idx = block - 1;
        
        ProfileStateR7 curr_ps = (ProfileStateR7) (col % PSR7_NUM_STATES);

        HMMAlignmentState as;
        as.event_idx = event_idx;
        as.kmer_idx = kmer_idx;
        as.l_posterior = -INFINITY; // not computed
        as.l_fm = get(vm, row, col);
        as.log_transition_probability = -INFINITY; // not computed
        as.state = ps2char(curr_ps);
        alignment.push_back(as);

        // Update the event (row) and k-mer using the current state
        // The next state is encoded in the backtrack matrix for the current cell
        ProfileStateR7 next_ps = (ProfileStateR7)get(bm, row, col);
        
        // If we hit the softclip state we are done aligning
        if(next_ps == PSR7_PRE_SOFT) {
            break;
        }

#if DEBUG_BACKTRACK
        printf("Backtrack [%zu %zu] k: %zu block: %zu curr_ps: %c next_ps: %c\n", row, col, kmer_idx, block, ps2char(curr_ps), ps2char(next_ps));
#endif

        if(curr_ps == PSR7_MATCH) {
            row -= 1;
            kmer_idx -= 1;
        } else if(curr_ps == PSR7_EVENT_SPLIT) {
            row -= 1;
            // kmer stays the same
        } else {
            assert(curr_ps == PSR7_KMER_SKIP);
            // row stays the same
            kmer_idx -= 1;
        }

        col = PSR7_NUM_STATES * (kmer_idx + 1) + next_ps;
    }


#if HMM_REVERSE_FIX
    // change the strand of the kmer indices if we aligned to the reverse strand
    if(data.event_stride == -1) {
        for(size_t ai = 0; ai < alignment.size(); ++ai) {
            size_t k_idx = alignment[ai].kmer_idx;
            alignment[ai].kmer_idx = sequence.length() - k_idx - k;
        }
    } else {
        std::reverse(alignment.begin(), alignment.end());
    }
#else
    //
    std::reverse(alignment.begin(), alignment.end());
#endif

    //
    free_matrix(vm);
    free_matrix(bm);

    return alignment;
}
