//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
// for R7 data
//
#ifndef NANOPOLISH_PROFILE_HMM_R7_H
#define NANOPOLISH_PROFILE_HMM_R7_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_matrix.h"
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_profile_hmm.h"

//#define HMM_REVERSE_FIX 1

//
// High level algorithms
//

// Calculate the probability of the nanopore events given a sequence
float profile_hmm_score_r7(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);

// Run viterbi to align events to kmers
std::vector<HMMAlignmentState> profile_hmm_align_r7(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);

//
// Forward algorithm
//

// Initialize the forward algorithm
void profile_hmm_forward_initialize_r7(FloatMatrix& fm);

// Terminate the forward algorithm
float profile_hmm_forward_terminate_r7(const FloatMatrix& fm, uint32_t row);

//
// Viterbi
// 

// initialize viterbi
void profile_hmm_viterbi_initialize_r7(FloatMatrix& m);

// Convenience enum for keeping track of the states in the profile HMM
enum ProfileStateR7
{
    PSR7_KMER_SKIP = 0,
    PSR7_EVENT_SPLIT,
    PSR7_MATCH,
    PSR7_NUM_STATES = 3,
    PSR7_PRE_SOFT // intentionally after PSR7_NUM_STATES
};

// Convert an enumerated state into a symbol
inline char ps2char(ProfileStateR7 ps) { return "KEMNS"[ps]; }

// Pre-computed transitions from the previous block
// into the current block of states. Log-scaled.
struct BlockTransitionsR7
{
    // Transition from m state
    float lp_me;
    float lp_mk;
    float lp_mm;

    // Transitions from e state
    float lp_ee;
    float lp_em;

    // Transitions from k state
    float lp_kk;
    float lp_km;
};

//
#include "nanopolish_profile_hmm_r7.inl"

#endif
