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
#define PACKAGE_VERSION "0.8.4"
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

// strands
const uint8_t T_IDX = 0;
const uint8_t C_IDX = 1;
const uint8_t NUM_STRANDS = 2;

//
// Data structures
//

// Forward declare
class SquiggleRead;
class PoreModel;

// This struct is used as input into the HMM
// It tracks where the event stream starts/stops for the candidate sequence
// that is being evaluated.
struct HMMInputData
{
    SquiggleRead* read;
    const PoreModel* pore_model;
    uint32_t event_start_idx;
    uint32_t event_stop_idx;
    uint8_t strand;
    int8_t event_stride;
    uint8_t rc;
};

// A representation of an event->kmer alignment
struct HMMAlignmentState
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
    GaussianParameters() : mean(0.0f), stdv(1.0f), log_stdv(0.0) { }
    GaussianParameters(float m, float s) : mean(m), stdv(s) { log_stdv = log(stdv); }

    float mean;
    float stdv;
    float log_stdv; // == log(stdv), pre-computed for efficiency
};

//
struct SemVer
{
    int major;
    int minor;
    int patch;
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

// returns true if the provided string ends with the given suffix
bool ends_with(const std::string& str, const std::string& suffix);

// parse a region string (chr:start-end)
void parse_region_string(const std::string& region, std::string& contig, int& start, int& end);

// parse a software version string using the semver convention
SemVer parse_semver_string(const std::string& semver_str);

// from: http://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
size_t nChoosek(size_t n, size_t k);

// print a warning message to stderr a single time
// this is only for debugging, please don't litter the code with them
#define WARN_ONCE(x) static bool _warn_once = true; if(_warn_once) \
                     fprintf(stderr, "WARNING: [%s]\n", (x)); _warn_once = false;

template<class T>
std::string array2str(const T& array)
{
    std::string s(array.data(), array.size());
    size_t null_pos = s.find('\0');
    if(null_pos != std::string::npos) {
        s.erase(null_pos);
    }
    return s;
}

// Print an error message when failing to load a bam index and exit
void bam_index_error_exit(const std::string& bam_filename);

#endif
