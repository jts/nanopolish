//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_methyltrain -- train a methylation model
//
#ifndef NANOPOLISH_METHYLTRAIN_H
#define NANOPOLISH_METHYLTRAIN_H

//#include <string>
//#include <map>
//#include "nanopolish_poremodel.h"
#include <vector>
#include "nanopolish_eventalign.h"
#include "nanopolish_squiggle_read.h"

// recalculate shift, scale, drift, scale_sd from an alignment and the read
// returns true if the recalibration was performed
// in either case, sets residual to the L1 norm of the residual
bool recalibrate_model(SquiggleRead &sr,
                       const PoreModel& pore_model,
                       const int strand_idx,
                       const std::vector<EventAlignment> &alignment_output,
                       bool scale_var=true,
                       bool scale_drift=true);

int methyltrain_main(int argc, char** argv);

#endif
