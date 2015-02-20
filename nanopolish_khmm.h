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

// TODO: refactor
struct HMMConsReadState;

void fill_khmm_transitions(DoubleMatrix& matrix, const std::string& consensus, const HMMConsReadState& state);

void initialize_forward_khmm(DoubleMatrix& fm);

double forward_khmm_terminate(const DoubleMatrix& fm, const DoubleMatrix& tm, uint32_t row);

double fill_forward_khmm(DoubleMatrix& fm, // forward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMConsReadState& state,
                         uint32_t e_start);

void initialize_backward_khmm(DoubleMatrix& bm, const DoubleMatrix& tm);

void fill_backward_khmm(DoubleMatrix& bm, // backward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMConsReadState& state,
                         uint32_t e_start);

double score_khmm_model(const std::string& consensus, const HMMConsReadState& state, AlignmentPolicy policy);


std::vector<PosteriorState> posterior_decode_khmm(const std::string& sequence, const HMMConsReadState& state);

void update_training_khmm(const std::string& consensus, 
                          const HMMConsReadState& state);

void debug_khmm_model(const std::string& name,
                      uint32_t seq_id,
                      uint32_t read_id,
                      const std::string& consensus, 
                      const HMMConsReadState& state);


#endif
