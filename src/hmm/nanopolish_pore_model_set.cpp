//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_pore_model_set -- A class that maintains
// a collection of pore models that SquiggleReads
// can load during initialization.
//
#include "nanopolish_pore_model_set.h"

PoreModelSet::~PoreModelSet()
{
}

PoreModel PoreModelSet::get_model_by_name(const std::string& name)
{
    PoreModelSet& model_set = getInstance();

    // look up the model in the cache
    auto const& iter = model_set.m_pore_model_cache.find(name);
    if(iter == model_set.m_pore_model_cache.end()) {
        // load the model from disk
        // this will intentially exit with an error if the model cannot be found
        std::cerr << "Loading model from disk: " << name << "\n";
        PoreModel model(name);
        model_set.m_pore_model_cache.insert(std::make_pair(name, model));
        return model;
    } else {
        return iter->second;
    }
}
