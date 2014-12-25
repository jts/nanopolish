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

    /* These are provided by the hdf5 file but currently unused
    double drift;
    double scale_sd;
    double var;
    double var_2d;
    */
    CPoreModelStateParams state[1024];
};


#endif
