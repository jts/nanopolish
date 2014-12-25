#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include "hmmcons_poremodel.h"

struct CSquiggleRead
{
    // unique identifier of the read
    uint32_t read_id;

    // one model for each strand
    CPoreModel pore_model[2];
};

// These structures are initialized by python
// and passed into the compiled code to prepare
// the data structures
struct CPoreModelParameters
{
    // Number of k-mers that are modelled
    int n_states;
    double* mean;
    double* sd;
    double scale;
    double shift;
};

struct CSquiggleReadParameters
{
    CPoreModelParameters pore_model[2];
};

// A global vector used to store data we've received from the python code
struct HmmConsData
{
    std::vector<CSquiggleRead> reads;
};
HmmConsData g_data;

/*
extern "C"
void print_reads()
{
    for(size_t i = 0; i < g_data.reads.size(); ++i)
        printf("read[%zu]: %s\n", i, g_data.reads[i].c_str());
}
*/

extern "C"
void add_read(CSquiggleReadParameters params)
{
    printf("Received %zu values from python\n", params.pore_model[0].n_states);
    for(size_t i = 0; i < params.pore_model[0].n_states; ++i)
        printf("state[%zu] model: %lf %lf %lf %lf\n", i, params.pore_model[0].mean[i], 
                                                         params.pore_model[0].sd[i],
                                                         params.pore_model[1].mean[i],
                                                         params.pore_model[1].sd[i]);
    printf("Scale: %lf %lf Shift: %lf %lf\n", params.pore_model[0].scale,
                                              params.pore_model[1].scale,
                                              params.pore_model[0].shift,
                                              params.pore_model[1].shift);
}

int main(int argc, char** argv)
{

}
