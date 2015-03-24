//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_emissions -- Emission distributions and
// related functions for the HMMs
//
#ifndef NANOPOLISH_EMISSIONS_H
#define NANOPOLISH_EMISSIONS_H

#include <math.h>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

//#define MODEL_STDV
//#define DEBUG_HMM_EMISSION 1

// From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
inline double normal_pdf(double x, const GaussianParameters& g)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - g.mean) / g.stdv;
    return inv_sqrt_2pi / g.stdv * exp(-0.5f * a * a);
}

inline double log_normal_pdf(double x, const GaussianParameters& g)
{
    static const double log_inv_sqrt_2pi = log(0.3989422804014327);
    double a = (x - g.mean) / g.stdv;
    return log_inv_sqrt_2pi - g.log_stdv + (-0.5f * a * a);
}

// The probability that a standard normal RV is <= x
// from http://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
inline double log_standard_normal_cdf(double x)
{
    return 0.5 * erfc(-x * M_SQRT1_2);
}

// The probability that a normal RV is <= x
inline double log_normal_cdf(double x, const GaussianParameters& g)
{
    double a = (x - g.mean) / g.stdv;
    return log(0.5 * (1 + erf(a * M_SQRT1_2)));
}

inline double log_probability_match(const SquiggleRead& read,
                                    uint32_t kmer_rank,
                                    uint32_t event_idx, 
                                    uint8_t strand,
                                    double state_scale = 1.0f,
                                    double log_state_scale = 0.0f)
{
    const PoreModel& pm = read.pore_model[strand];

    // event level mean
    double level = read.get_drift_corrected_level(event_idx, strand);

    GaussianParameters model = pm.get_scaled_parameters(kmer_rank);

    // we go to great lengths to avoid calling log() in the inner loop of the HMM
    // for this reason we duplicate data here and require the caller to pass
    // in the scale and log(scale), presumably these are cached
    model.stdv *= state_scale;
    model.log_stdv += log_state_scale;
    double lp = log_normal_pdf(level, model);

#if MODEL_STDV
    // event level stdv
    double stdv = read.events[strand].stdv[event_idx];
    double model_sd_mean = pm.state[kmer_rank].sd_mean * pm.scale_sd;
    double model_sd_stdv = pm.state[kmer_rank].sd_stdv * sqrt(pow(pm.scale_sd, 3.0) / pm.var_sd);
    lp += log_normal_pdf(stdv, model_sd_mean, model_sd_stdv);
#endif

#if DEBUG_HMM_EMISSION
    printf("Event[%d] Kmer: %d -- L:%.1lf m: %.1lf s: %.1lf p: %.3lf p_old: %.3lf\n", event_idx, kmer_rank, level, model.mean, model.stdv, exp(lp), normal_pdf(level, model));
#endif

    return lp;
}

inline double log_probability_event_insert(const SquiggleRead& read,
                                           uint32_t kmer_rank,
                                           uint32_t event_idx, 
                                           uint8_t strand)
{
    static const double scale = 1.75f;
    static const double log_scale = log(scale);

    return log_probability_match(read, kmer_rank, event_idx, strand, scale, log_scale);
}

inline double log_probability_kmer_insert(const SquiggleRead& read,
                                          uint32_t kmer_rank,
                                          uint32_t event_idx, 
                                          uint8_t strand)

{
    return log_probability_match(read, kmer_rank, event_idx, strand);
}

#endif
