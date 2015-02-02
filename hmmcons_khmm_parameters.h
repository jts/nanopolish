// TODO: Boilerplate
// This file implements a data structure and functions
// for storing and learning transition probabilities
// for the KHMM model
#ifndef HMMCONS_KHMM_PARAMETERS
#define HMMCONS_KHMM_PARAMETERS

#include <vector>
#include <stdint.h>

// 
struct TransitionObservation
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

    std::vector<TransitionObservation> transitions;
    std::vector<double> emissions_for_matches;
};

//
struct KHMMParameters
{
    // The probability of staying in the same state in the HMM
    double self_transition;

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

// Train the model
void khmm_parameters_train(KHMMParameters& parameters);

// Get the probability of skipping a kmer observation given the expected signal level
double get_skip_probability(const KHMMParameters& parameters, double k_level1, double k_level2);

#endif
