//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_progressive_align -- align squiggle events
// to kmers of a sequence using a banded,
// progressive algorithm
//
#ifndef NANOPOLISH_PROGRESSIVE_ALIGN_H
#define NANOPOLISH_PROGRESSIVE_ALIGN_H

#include "nanopolish_squiggle_read.h"
#include "scrappie_structures.h"

void estimate_scalings_using_mom(const std::string& sequence,
                                 const PoreModel& pore_model,
                                 const event_table& et,
                                 double& out_shift,
                                 double& out_scale);

void progressive_align(SquiggleRead& read,
                       const std::string& sequence);

#endif
