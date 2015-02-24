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

// Calculate the probability of the nanopore events given the consensus sequence
double profile_hmm_score(const std::string& consensus, const HMMConsReadState& state);

// Initialize the forward algorithm
void profile_hmm_forward_initialize(DoubleMatrix& fm);

// Terminate the forward algorithm
double profile_hmm_forward_terminate(const DoubleMatrix& fm, uint32_t row);

// Fill in the forward matrix
double profile_hmm_forward_fill(DoubleMatrix& fm, // forward matrix
                                const char* sequence,
                                const HMMConsReadState& state,
                                uint32_t e_start);

#endif
