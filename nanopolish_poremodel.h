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

// The pore model is defined by global scale/shift parameters
// and a mean/stddev per k-mer. These are parameterize the
// Gaussian PDF.
struct PoreModelStateParams
{
    double level_mean;
    double level_stdv;
    
    double sd_mean;
    double sd_stdv;
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
            GaussianParameters ret = { scaled_state[kmer_rank].level_mean,
                                       scaled_state[kmer_rank].level_stdv };
            return ret;
        }

        // Pre-compute the scaled and shifted pore model states
        // to avoid doing it in the loop of the HMM
        void transform();
    
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
        PoreModelStateParams scaled_state[PORE_MODEL_STATES];
};

#endif
