//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_model_names -- Get metadata for ONT model names
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include "nanopolish_model_names.h"
#include "nanopolish_common.h"

static std::string strand_by_idx[] = { "template", "complement.pop1", "complement.pop2" };
static std::string short_strand_by_idx[] = { "t", "c.p1", "c.p2" };
static std::string kit_name_by_idx[] = { "SQK005", "SQK006", "r9_250bps", "r9.4_450bps", "r9.4_70bps" };

static std::map< std::string, ModelMetadata > known_models = {

    // SQK005 models
    { "r7.3_template_median68pA.model", {T_IDX, 0, KV_SQK005 } },
    { "r7.3_complement_median68pA_pop1.model", {C_IDX, 1, KV_SQK005 } },
    { "r7.3_complement_median68pA_pop2.model", {C_IDX, 2, KV_SQK005 } },
    
    // SQK006 models
    { "r7.3_e6_70bps_6mer_template_median68pA.model", {T_IDX, 0, KV_SQK006 } },
    { "r7.3_e6_70bps_6mer_complement_median68pA_pop1.model", {C_IDX, 1, KV_SQK006 } },
    { "r7.3_e6_70bps_6mer_complement_median68pA_pop2.model", {C_IDX, 2, KV_SQK006 } },
};

std::string ModelMetadata::get_short_name() const
{
    return get_kit_name() + "." + short_strand_by_idx[this->model_idx];
}

std::string ModelMetadata::get_kit_name() const
{
    assert(this->kit < NUM_KITS);
    return kit_name_by_idx[this->kit];
}

std::string ModelMetadata::get_strand_model_name() const
{
    assert(this->model_idx < 3);
    return strand_by_idx[this->model_idx];
}

ModelMetadata get_model_metadata_from_name(const std::string& name)
{
    auto iter = known_models.find(name);
    if (iter != known_models.end()) {
        return iter->second;
    } else {
        std::cerr << "Error: unknown model: " << name << "\n";
        exit(EXIT_FAILURE);
    }
}
