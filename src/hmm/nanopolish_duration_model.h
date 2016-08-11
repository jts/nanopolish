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

struct GammaParameters
{
    double shape;
    double rate;
};

class DurationModel
{
    public:
        DurationModel();

        //
        static std::vector<double> generate_aligned_durations(const std::string& sequence,
                                                              const HMMInputData& data,
                                                              const uint32_t alignment_flags);
        //
        // Log of gamma PDF for the sum of n observations
        //
        static double log_gamma_sum(double x, double n);
        static double log_gamma_sum(double x, const GammaParameters& params, double n);

        //
        // Fit the parameters of the gamma distribution for the input data
        //
        static GammaParameters gamma_fit(const std::vector<double>& input);
        static double gamma_fit_calculate_s(const std::vector<double>& input);

    private:
        
        // singleton accessor function
        static DurationModel& getInstance()
        {
            static DurationModel instance;
            return instance;
        }
};

#endif
