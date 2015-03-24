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

        // these functions are provided by ONT
        scaled_params[i].mean = state[i].level_mean * scale + shift;
        scaled_params[i].stdv = state[i].level_stdv * var;
        scaled_params[i].log_stdv = log(scaled_params[i].stdv); // pre-computed for efficiency

        // These are not used, for now
        //scaled_state[i].sd_mean = state[i].sd_mean * scale_sd;
        //scaled_state[i].sd_stdv = state[i].sd_stdv * sqrt(pow(scale_sd, 3.0) / var_sd);
    }
    is_scaled = true;
}
