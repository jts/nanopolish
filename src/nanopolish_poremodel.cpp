//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#include "nanopolish_poremodel.h"

void PoreModel::bake_gaussian_parameters()
{
    for(int i = 0; i < PORE_MODEL_STATES; ++i) {

        // as per ONT documents
        scaled_state[i].level_mean = state[i].level_mean * scale + shift;
        scaled_state[i].level_stdv = state[i].level_stdv * var;

        scaled_state[i].sd_mean = state[i].sd_mean * scale_sd;
        scaled_state[i].sd_lambda = state[i].sd_lambda * var_sd;
        scaled_state[i].sd_stdv = sqrt(pow(scaled_state[i].sd_mean, 3) / scaled_state[i].sd_lambda);

        // for efficiency
        scaled_state[i].level_log_stdv = log(scaled_state[i].level_stdv);
        scaled_state[i].sd_log_lambda = log(scaled_state[i].sd_lambda);

        // for compatibility
        scaled_params[i].mean = scaled_state[i].level_mean;
        scaled_params[i].stdv = scaled_state[i].level_stdv;
        scaled_params[i].log_stdv = scaled_state[i].level_log_stdv;
    }
    is_scaled = true;
}
