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

//
// High level algorithms
//

// Calculate the probability of the nanopore events given the consensus sequence
float profile_hmm_score(const std::string& consensus, const HMMInputData& data);

// Run viterbi to align events to kmers
std::vector<AlignmentState> profile_hmm_align(const std::string& sequence, const HMMInputData& data);

//
// Forward algorithm
//

// Initialize the forward algorithm
void profile_hmm_forward_initialize(FloatMatrix& fm);

// Terminate the forward algorithm
float profile_hmm_forward_terminate(const FloatMatrix& fm, uint32_t row);

// Fill in the forward matrix
float profile_hmm_forward_fill(FloatMatrix& fm, // forward matrix
                                const char* sequence,
                                const HMMInputData& data,
                                uint32_t e_start);


//
// Viterbi
// 

// initialize viterbi
void profile_hmm_viterbi_initialize(FloatMatrix& m);

// fill in the score matrix and the backtrack matrix
void profile_hmm_viterbi_fill(FloatMatrix& vm, // viterbi matrix
                              UInt8Matrix& bm, // backtrack matrix
                              const char* sequence,
                              const HMMInputData& data,
                              uint32_t e_start);

//
// Training
//
void profile_hmm_update_training(const std::string& consensus, 
                                 const HMMInputData& data);

#endif
