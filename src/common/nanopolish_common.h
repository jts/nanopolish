//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_common -- Data structures and definitions
// shared across files
//
#ifndef NANOPOLISH_COMMON_H
#define NANOPOLISH_COMMON_H

#include <stdint.h>
#include <string>
#include <vector>
#include <math.h>
#include "nanopolish_alphabet.h"
#include "profiler.h"
#include "logsum.h"

#define PACKAGE_NAME "nanopolish"
#define PACKAGE_VERSION "0.3.0"
#define PACKAGE_BUGREPORT "https://github.com/jts/nanopolish/issues"

//
// Enumerated types
//
enum AlignmentPolicy
{
    AP_GLOBAL,
    AP_SEMI_KMER
};

//
// Constants
//
const uint8_t K = 5;

// strands
const uint8_t T_IDX = 0;
const uint8_t C_IDX = 1;
const uint8_t NUM_STRANDS = 2;

//
// Data structures
//

class SquiggleRead;

// This struct is used as input into the HMM
// It tracks where the event stream starts/stops
// for the partial consensus sequence under consideration
struct HMMInputData
{
    SquiggleRead* read;
    uint32_t anchor_index;
    uint32_t event_start_idx;
    uint32_t event_stop_idx;
    uint8_t strand;
    int8_t event_stride;
    uint8_t rc;
};

// A representation of an event->kmer alignment
struct AlignmentState
{
    uint32_t event_idx;
    uint32_t kmer_idx;
    double l_posterior;
    double l_fm;
    double log_transition_probability;
    char state;
};

// The parameters of a gaussian distribution
struct GaussianParameters
{
    GaussianParameters() : mean(0.0f), stdv(1.0f) { log_stdv = log(stdv); }
    GaussianParameters(float m, float s) : mean(m), stdv(s) { log_stdv = log(stdv); }

    float mean;
    float stdv;
    float log_stdv; // == log(stdv), pre-computed for efficiency
};

//
// Functions
//
#define ESL_LOG_SUM 1

// Add the log-scaled values a and b using a transform to avoid precision errors
inline double add_logs(const double a, const double b)
{
#if ESL_LOG_SUM
    return p7_FLogsum(a, b);
#else
    if(a == -INFINITY && b == -INFINITY)
        return -INFINITY;

    if(a > b) {
        double diff = b - a;
        return a + log(1.0 + exp(diff));
    } else {
        double diff = a - b;
        return b + log(1.0 + exp(diff));
    }
#endif
}

// split a string based on a delimiter
std::vector<std::string> split(std::string in, char delimiter);

#endif
