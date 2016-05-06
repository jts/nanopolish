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

static std::map< std::string, ModelMetadata > known_models = {

    // SQK005 models
    { "r7.3_template_median68pA.model", {T_IDX, 0, KV_SQK005, "t.005"} },
    { "r7.3_complement_median68pA_pop1.model", {C_IDX, 1, KV_SQK005, "c.p1.005"} },
    { "r7.3_complement_median68pA_pop2.model", {C_IDX, 2, KV_SQK005, "c.p2.005"} },
    
    // SQK006 models
    { "r7.3_e6_70bps_6mer_template_median68pA.model", {T_IDX, 0, KV_SQK006, "t.006"} },
    { "r7.3_e6_70bps_6mer_complement_median68pA_pop1.model", {C_IDX, 1, KV_SQK006, "c.p1.006"} },
    { "r7.3_e6_70bps_6mer_complement_median68pA_pop2.model", {C_IDX, 2, KV_SQK006, "c.p2.006"} },

    // SQK007 models
    { "r9.template.model", {T_IDX, 0, KV_SQK007, "t.007"} },
    { "r9.template.5mer.base.model", {T_IDX, 0, KV_SQK007, "t.007"} }

};

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
