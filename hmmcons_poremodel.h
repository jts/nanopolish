// TODO: Boilerplate

#ifndef HMMCONS_POREMODEL_H
#define HMMCONS_POREMODEL_H

//
// The pore model is defined by global scale/shift parameters
// and a mean/stddev per k-mer. These are parameterize the
// Gaussian PDF.
struct CPoreModelStateParams
{
    double mean;
    double sd;
};

//
struct CPoreModel
{
    double scale;
    double shift;
    CPoreModelStateParams state[1024];
};

#endif
