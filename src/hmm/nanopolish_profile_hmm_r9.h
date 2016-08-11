//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm_r9 -- Profile Hidden Markov Model
// for R9 data
//
#ifndef NANOPOLISH_PROFILE_HMM_R9_H
#define NANOPOLISH_PROFILE_HMM_R9_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_matrix.h"
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_profile_hmm.h"

//#define HMM_REVERSE_FIX 1
//#define DEBUG_FILL 1

//
// High level algorithms
//

// Calculate the probability of the nanopore events given a sequence
float profile_hmm_score_r9(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);

// Run viterbi to align events to kmers
std::vector<HMMAlignmentState> profile_hmm_align_r9(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);

//
// Forward algorithm
//

// Initialize the forward algorithm
void profile_hmm_forward_initialize_r9(FloatMatrix& fm);

// Terminate the forward algorithm
float profile_hmm_forward_terminate_r9(const FloatMatrix& fm, uint32_t row);

//
// Viterbi
//

// initialize viterbi
void profile_hmm_viterbi_initialize_r9(FloatMatrix& m);

// Convenience enum for keeping track of the states in the profile HMM
enum ProfileStateR9
{
    PSR9_KMER_SKIP = 0,
    PSR9_BAD_EVENT,
    PSR9_MATCH,
    PSR9_NUM_STATES = 3,
    PSR9_PRE_SOFT // intentionally after PS_NUM_STATES
};

enum HMMMovementType
{
    HMT_FROM_SAME_M = 0,
    HMT_FROM_PREV_M,
    HMT_FROM_SAME_B,
    HMT_FROM_PREV_B,
    HMT_FROM_PREV_K,
    HMT_FROM_SOFT,
    HMT_NUM_MOVEMENT_TYPES
};
typedef struct { float x[HMT_NUM_MOVEMENT_TYPES]; } HMMUpdateScores;

// Convert an enumerated state into a symbol
inline char ps2char(ProfileStateR9 ps) { return "KBMNS"[ps]; }

// Pre-computed transitions from the previous block
// into the current block of states. Log-scaled.
struct BlockTransitions
{
    // Transition from m state (match event to k-mer)
    float lp_mm_self;
    float lp_mb;
    float lp_mk;
    float lp_mm_next;

    // Transitions from b state (bad event that should be ignored)
    float lp_bb;
    float lp_bk;
    float lp_bm_next; // movement to next k-mer
    float lp_bm_self; // movement to k-mer that we came from

    // Transitions from k state (no observation from k-mer)
    float lp_kk;
    float lp_km;
};

//
#include "nanopolish_profile_hmm_r9.inl"

#endif
