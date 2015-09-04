//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#ifndef NANOPOLISH_POREMODEL_H
#define NANOPOLISH_POREMODEL_H

#include <assert.h>
#include "nanopolish_common.h"

#define PORE_MODEL_STATES 1024

//
struct PoreModelStateParams
{
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;

    double level_log_stdv;
    double sd_lambda;
    double sd_log_lambda;
};

//
class PoreModel
{
    public:

        //
        PoreModel() : is_scaled(false) {}

        inline GaussianParameters get_scaled_parameters(const uint32_t kmer_rank) const
        {
            assert(is_scaled);
            return scaled_params[kmer_rank];
        }

        inline PoreModelStateParams get_parameters(const uint32_t kmer_rank) const
        {
            return state[kmer_rank];
        }
        
        // Pre-compute the GaussianParameters to avoid
        // taking numerous logs in the emission calculations
        void bake_gaussian_parameters();
    
    friend SquiggleRead;

    private:

        double scale;
        double shift;
        double drift;
        double var;
        double scale_sd;
        double var_sd;
        bool is_scaled;

        //
        PoreModelStateParams state[PORE_MODEL_STATES];
        GaussianParameters scaled_params[PORE_MODEL_STATES];
};

#endif
