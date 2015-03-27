//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_khmm_parameters -- parameters for khmm model
//
#ifndef NANOPOLISH_KHMM_PARAMETERS
#define NANOPOLISH_KHMM_PARAMETERS

#include <vector>
#include <stdint.h>
#include "nanopolish_matrix.h"

// 
struct KmerTransitionObservation
{
    double level_1;
    double level_2;
    char state;
};

// This struct holds observations used to learn the parameters
struct TrainingData
{
    uint32_t n_matches;
    uint32_t n_merges;
    uint32_t n_skips;

    std::vector<KmerTransitionObservation> kmer_transitions;
    std::vector<double> emissions_for_matches;
    UInt32Matrix state_transitions;
};

//
struct KHMMParameters
{
    // The probability of staying in the same state in the HMM
    double self_transition;
    
    double trans_m_to_e_not_k;
    double trans_e_to_e;

    // This is a vector that maps from discretized absolute difference
    // between expected signals to a probability that the transition
    // will be observed by the pore. Access to the skip probabilities
    // for a pair of k-mer levels is through get_skip_probability()
    std::vector<double> skip_probabilities;
    double skip_bin_width;

    // Data used to train the model
    TrainingData training_data;

    double fit_quality;
};

// Initialize the parameters to some default values
void khmm_parameters_initialize(KHMMParameters& parameters);

// Clean up parameters
void khmm_parameters_destroy(KHMMParameters& parameters);

// add an observation of a state transition
void add_state_transition(TrainingData& td, char from, char to);

// Train the model
void khmm_parameters_train(KHMMParameters& parameters);

// Get the probability of skipping a kmer observation given the expected signal level
double get_skip_probability(const KHMMParameters& parameters, double k_level1, double k_level2);

#endif
