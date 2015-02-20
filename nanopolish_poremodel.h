// TODO: Boilerplate

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
    PoreModelStateParams state[1024];
};

#endif
