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

//
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
struct PoreModel
{
    double scale;
    double shift;
    double drift;
    double var;
    double scale_sd;
    double var_sd;
    PoreModelStateParams state[1024];
};

#endif
