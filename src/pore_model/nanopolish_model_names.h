//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_model_names -- Get metadata for ONT model names
//
#ifndef NANOPOLISH_MODEL_NAMES_H
#define NANOPOLISH_MODEL_NAMES_H

#include <stdint.h>
#include <string>
#include <vector>
#include <math.h>

enum KitVersion
{
    KV_SQK005 = 0,
    KV_SQK006,
    KV_R9_250BPS,
    KV_R9_4_450BPS,
    KV_R9_4_70BPS,
    NUM_KITS
};

// The parameters of a gaussian distribution
struct ModelMetadata
{
    uint8_t strand_idx;
    uint8_t model_idx; // template = 0, pop1 = 1, pop2 = 2
    KitVersion kit;

    // functions
    std::string get_short_name() const;
    std::string get_kit_name() const;
    std::string get_strand_model_name() const;

    bool is_r9() const { return kit >= KV_R9_250BPS; }
    bool is_r7() const { return kit >= KV_SQK005 && kit <= KV_SQK006; }
};

//
// Functions
//
ModelMetadata get_model_metadata_from_name(const std::string& name);


#endif
