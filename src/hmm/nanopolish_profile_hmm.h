//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
//
#ifndef NANOPOLISH_PROFILE_HMM_H
#define NANOPOLISH_PROFILE_HMM_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_matrix.h"
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"
#include "nanopolish_hmm_input_sequence.h"

//
// High level algorithms
//

// Calculate the probability of the nanopore events given a sequence
float profile_hmm_score(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);
float profile_hmm_score(const HMMInputSequence& sequence, const std::vector<HMMInputData>& data, const uint32_t flags = 0);

// Calculate the probability of the nanopore events given a set of possible sequences (usually methylated alternatives)
float profile_hmm_score_set(const std::vector<HMMInputSequence>& sequence, const HMMInputData& data, const uint32_t flags = 0);

// Run viterbi to align events to kmers
std::vector<HMMAlignmentState> profile_hmm_align(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags = 0);

// Flags to modify the behaviour of the HMM
enum HMMAlignmentFlags
{
    HAF_ALLOW_PRE_CLIP = 1, // allow events to go unmatched before the aligning region
    HAF_ALLOW_POST_CLIP = 2 // allow events to go unmatched after the aligning region
};

#endif
