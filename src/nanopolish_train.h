//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_train -- train k-mer models
//
#ifndef NANOPOLISH_TRAIN_H
#define NANOPOLISH_TRAIN_H

#include <vector>
#include "nanopolish_eventalign.h"
#include "nanopolish_squiggle_read.h"

int train_main(int argc, char** argv);

#endif
