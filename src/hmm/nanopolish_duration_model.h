//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_duration_model -- Model the duration
// of bases passing through the pore
//
#ifndef NANOPOLISH_DURATION_MODEL_H
#define NANOPOLISH_DURATION_MODEL_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_common.h"

#define MIN_DURATION 0.00025
#define MAX_INDEX 99

class DurationModel
{
    public:
        DurationModel();

        //
        // Align the events to the reference and accumulate training data for the distribution 
        //
        static void add_training_from_sequence(const std::string& sequence, 
                                               const HMMInputData& data,
                                               const uint32_t alignment_flags);

        //
        // Calculate the probability of the duration observations given the sequence
        //
        static double score_sequence(const std::string& sequence, 
                                     const HMMInputData& data,
                                     const uint32_t alignment_flags);
                                     
        static std::vector<double> score_vector_sequence(const std::string& sequence, 
                                                         const HMMInputData& data,
                                                         const uint32_t alignment_flags);


        //
        // Create a vector containing the total duration of events aligned to each k-mer
        //
        static std::vector<double> generate_aligned_durations(const std::string& sequence,
                                                              const HMMInputData& data,
                                                              const uint32_t alignment_flags);

        //
        // Log of gamma PDF for the sum of n observations
        //
        static double log_gamma_sum(double x, double n);

        //
        // Print the training data to fp
        //
        static void print_training_data(FILE* fp);
        
        //
        // Print the trained model to fp
        //
        static void print_compiled_model(FILE* fp);

    private:
        
        static int get_index(double duration) { 
            int index = duration / MIN_DURATION; 
            return index < MAX_INDEX ? index : MAX_INDEX; 
        }

        // singleton accessor function
        static DurationModel& getInstance()
        {
            static DurationModel instance;
            return instance;
        }

        std::vector<double> training_data;
        std::vector<double> probability_by_index;
};

#endif
