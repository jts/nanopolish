//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_raw_loader - utilities and helpers for loading
// data directly from raw nanopore files without events
//
#ifndef NANOPOLISH_RAW_LOADER_H
#define NANOPOLISH_RAW_LOADER_H

#include "nanopolish_squiggle_read.h"
#include "scrappie_structures.h"

void estimate_scalings_using_mom(const std::string& sequence,
                                 const PoreModel& pore_model,
                                 const event_table& et,
                                 double& out_shift,
                                 double& out_scale);

void banded_simple_event_align(SquiggleRead& read,
                               const std::string& sequence);

#endif
