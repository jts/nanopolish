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

//#define DEBUG_HMM_EMISSION 1

// From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
static const float inv_sqrt_2pi = 0.3989422804014327;
inline float normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return inv_sqrt_2pi / g.stdv * exp(-0.5f * a * a);
}

inline float normal_pdf(float x, const PoreModelStateParams& s)
{
    float a = (x - s.level_mean) / s.level_stdv;
    return inv_sqrt_2pi / s.level_stdv * exp(-0.5f * a * a);
}

inline float z_score(const SquiggleRead& read,
                     const PoreModel& pore_model,
                     uint32_t kmer_rank,
                     uint32_t event_idx,
                     uint8_t strand)
{
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    return (level - gp.mean) / gp.stdv;
}

static const float log_inv_sqrt_2pi = log(0.3989422804014327);

inline float log_normal_pdf(float x, const PoreModelStateParams& s)
{
    float a = (x - s.level_mean) / s.level_stdv;
    return log_inv_sqrt_2pi - s.level_log_stdv + (-0.5f * a * a);
}

inline float log_normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return log_inv_sqrt_2pi - g.log_stdv + (-0.5f * a * a);
}

inline float log_probability_match_r9(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float lp = log_normal_pdf(level, gp);
    return lp;
}

inline float log_probability_match_r7(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand,
                                      float state_scale = 1.0f,
                                      float log_state_scale = 0.0f)
{
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    gp.stdv *= state_scale;
    gp.log_stdv += log_state_scale;
    float lp = log_normal_pdf(level, gp);
    return lp;
}

inline float log_probability_event_insert_r7(const SquiggleRead& read,
                                             const PoreModel& pore_model,
                                             uint32_t kmer_rank,
                                             uint32_t event_idx,
                                             uint8_t strand)
{
    static const float scale = 1.75f;
    static const float log_scale = log(scale);

    return log_probability_match_r7(read, pore_model, kmer_rank, event_idx, strand, scale, log_scale);
}

inline float log_probability_background(const SquiggleRead&,
                                        uint32_t,
                                        uint8_t)
{
    return -3.0f;
}

#endif
