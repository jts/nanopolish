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
                     uint32_t kmer_rank,
                     uint32_t event_idx,
                     uint8_t strand)
{
    const PoreModel& pm = read.pore_model[strand];
    float level = read.get_drift_corrected_level(event_idx, strand);
    GaussianParameters model = pm.get_scaled_parameters(kmer_rank);
    return (level - model.mean) / model.stdv;
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

inline float log_invgauss_pdf(float x, float log_x, const PoreModelStateParams& s)
{
    static const float log_2pi = log(2 * M_PI);
    float a = (x - s.sd_mean) / s.sd_mean;
    return (s.sd_log_lambda - log_2pi - 3 * log_x - s.sd_lambda * a * a / x) / 2;
}

inline bool& model_stdv()
{
    static bool _model_stdv = false;
    return _model_stdv;
}

inline float log_probability_match(const SquiggleRead& read,
                                   uint32_t kmer_rank,
                                   uint32_t event_idx,
                                   uint8_t strand,
                                   float state_scale = 1.0f,
                                   float log_state_scale = 0.0f)
{
    const PoreModel& pm = read.pore_model[strand];

    // event level mean
    float level = read.get_drift_corrected_level(event_idx, strand);

    PoreModelStateParams state = pm.get_scaled_state(kmer_rank);

    // we go to great lengths to avoid calling log() in the inner loop of the HMM
    // for this reason we duplicate data here and require the caller to pass
    // in the scale and log(scale), presumably these are cached
    state.level_stdv *= state_scale;
    state.level_log_stdv += log_state_scale;
    float lp = log_normal_pdf(level, state);

    if(model_stdv())
    {
        float stdv = read.get_stdv(event_idx, strand);
        float log_stdv = read.get_log_stdv(event_idx, strand);
        float lp_stdv = log_invgauss_pdf(stdv, log_stdv, state);
        lp += lp_stdv;
    }

#if DEBUG_HMM_EMISSION
    printf("Event[%d] Kmer: %d -- L:%.1lf m: %.1lf s: %.1lf p: %.3lf p_old: %.3lf\n", event_idx, kmer_rank, level, state.level_mean, state.level_stdv, exp(lp), normal_pdf(level, state));
#endif

    return lp;
}

inline float log_probability_event_insert(const SquiggleRead& read,
                                          uint32_t kmer_rank,
                                          uint32_t event_idx,
                                          uint8_t strand)
{
    static const float scale = 1.0f;
    static const float log_scale = log(scale);

    return log_probability_match(read, kmer_rank, event_idx, strand, scale, log_scale);
}

inline float log_probability_background(const SquiggleRead&,
                                        uint32_t,
                                        uint8_t)
{
    return -3.0f;
}


inline float log_probability_kmer_insert(const SquiggleRead& read,
                                         uint32_t kmer_rank,
                                         uint32_t event_idx,
                                         uint8_t strand)

{
    return log_probability_match(read, kmer_rank, event_idx, strand);
}

#endif
