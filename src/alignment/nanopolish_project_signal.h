//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_project_signal -- project signal data onto a reference genome
//
#ifndef NANOPOLISH_PROJECT_SIGNAL_H
#define NANOPOLISH_PROJECT_SIGNAL_H

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "nanopolish_alphabet.h"
#include "nanopolish_common.h"

// Entry point from nanopolish.cpp
int project_signal_main(int argc, char** argv);

#endif
