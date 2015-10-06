//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_methyltrain -- train a methylation model
//
#ifndef NANOPOLISH_METHYLTRAIN_H
#define NANOPOLISH_METHYLTRAIN_H

#include <string>
#include <map>
#include "nanopolish_poremodel.h"

//
// Typedefs
//
typedef std::map<std::string, PoreModel> ModelMap;

// read models from a file
ModelMap read_models_fofn(const std::string& fofn_name);

int methyltrain_main(int argc, char** argv);

#endif
