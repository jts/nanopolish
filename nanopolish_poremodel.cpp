//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#include "nanopolish_poremodel.h"

void PoreModel::transform()
{
    for(int i = 0; i < PORE_MODEL_STATES; ++i) {
        // these are provided by ONT
        scaled_state[i].level_mean = state[i].level_mean * scale + shift;
        scaled_state[i].level_stdv = state[i].level_stdv * var;
        scaled_state[i].sd_mean = state[i].sd_mean * scale_sd;
        scaled_state[i].sd_stdv = state[i].sd_stdv * sqrt(pow(scale_sd, 3.0) / var_sd);
    }
    is_scaled = true;
}
