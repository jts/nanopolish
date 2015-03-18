//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_khmm -- Hidden Markov Model with
// a k-mer sequence as the hidden state
//
#ifndef NANOPOLISH_KHMM_H
#define NANOPOLISH_KHMM_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_matrix.h"
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"

// Fill in the transition matrix for the model
void khmm_fill_transitions(DoubleMatrix& matrix, const std::string& consensus, const HMMInputData& data);

// Initialize the forward algorithm
void khmm_forward_initialize(DoubleMatrix& fm);

// Terminate the forward algorithm
double khmm_forward_terminate(const DoubleMatrix& fm, const DoubleMatrix& tm, uint32_t row);

// Fill in the forward matrix
double khmm_forward_fill(DoubleMatrix& fm, // forward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMInputData& data,
                         uint32_t e_start);

// Initialize the backward algorithm
void khmm_backward_initialize(DoubleMatrix& bm, const DoubleMatrix& tm);

// Fill in the backward matrix
void khmm_backward_fill(DoubleMatrix& bm, // backward matrix
                        const DoubleMatrix& tm, //transitions
                        const char* sequence,
                        const HMMInputData& data,
                        uint32_t e_start);

// Calculate the probability of the nanopore events given the consensus sequence
double khmm_score(const std::string& consensus, const HMMInputData& data, AlignmentPolicy policy);

// Perform posterior decoding of the path through the hidden states
std::vector<AlignmentState> khmm_posterior_decode(const std::string& sequence, const HMMInputData& data);

// Update the training data
void khmm_update_training(const std::string& consensus, 
                          const HMMInputData& data);

// Debug the model, printing hopefully useful information
void khmm_debug(const std::string& name,
                uint32_t seq_id,
                uint32_t read_id,
                const std::string& consensus, 
                const HMMInputData& data);


#endif
