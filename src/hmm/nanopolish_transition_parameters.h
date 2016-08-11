//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_transitions_parameters -- transition
// parameters for profile HMM
//
#ifndef NANOPOLISH_TRANSITION_PARAMETERS
#define NANOPOLISH_TRANSITION_PARAMETERS

#include <vector>
#include <stdint.h>
#include "nanopolish_matrix.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_model_names.h"

// 
struct KmerTransitionObservation
{
    double level_1;
    double level_2;
    char state;
};

// This struct holds observations used to learn the parameters
struct TransitionTrainingData
{
    uint32_t n_matches;
    uint32_t n_merges;
    uint32_t n_skips;

    std::vector<KmerTransitionObservation> kmer_transitions;
    UInt32Matrix state_transitions;
};

//
class TransitionParameters
{
    public:

        //
        // functions
        //
        TransitionParameters();
        ~TransitionParameters();

        void initialize(const ModelMetadata& metadata);

        // update transition parameters from training data
        void train();

        // Get the probability of skipping a kmer observation given the pair of expected levels
        double get_skip_probability(double k_level1, double k_level2) const;

        // add an observation of a state transition to the training data
        void add_transition_observation(char hmm_state_from, char hmm_state_to, bool kmer_move);

        // update the training data using the alignment
        void add_training_from_alignment(const HMMInputSequence& sequence,
                                         const HMMInputData& data,
                                         const std::vector<HMMAlignmentState>& alignment,
                                         size_t ignore_edge_length = 5);

        void print() const;

        //
        // data
        //
        double trans_m_to_e_not_k;
        double trans_e_to_e;

        double trans_start_to_clip;
        double trans_clip_self;
        
        bool is_initialized = false;

        // This is a vector that maps from discretized absolute difference
        // between expected signals to a probability that the transition
        // will be observed by the pore. Access to the skip probabilities
        // for a pair of k-mer levels is through get_skip_probability()
        std::vector<double> skip_probabilities;
        double skip_bin_width;

        // Data used to train the model
        TransitionTrainingData training_data;

    private:

        // Model-specific transition initialization
        void initialize_sqkmap005();
        void initialize_sqkmap006_template();
        void initialize_sqkmap006_complement();
        void initialize_sqkmap007_template();
        void initialize_sqkmap007_complement();

        // Not allowed
        TransitionParameters(const TransitionParameters&) {}

        // Calculate which bin of the skip probability table this level difference falls in
        inline size_t get_skip_bin(double k_level1, double k_level2) const
        {
            assert(!skip_probabilities.empty());

            double d = fabs(k_level1 - k_level2);
            size_t bin = d / skip_bin_width;

            // clamp out-of-range to last value
            bin = bin >= skip_probabilities.size() ? skip_probabilities.size() - 1 : bin;
            return bin;
        }

};


#endif
