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

//
PoreModelSet::~PoreModelSet()
{
}

//
void PoreModelSet::initialize(const std::string& fofn_filename)
{
    // grab singleton instance
    PoreModelSet& model_set = getInstance();
    
    // open the fofn file reader
    std::ifstream fofn_reader(fofn_filename);
    if(!fofn_reader.is_open()) {
        fprintf(stderr, "Error: could not read %s\n", fofn_filename.c_str());
        exit(EXIT_FAILURE);
    }

    std::string model_filename;
    while(getline(fofn_reader, model_filename)) {

        // read the model
        PoreModel p(model_filename);
        assert(!p.name.empty());
        assert(!p.type.empty());

        model_set.model_type_sets[p.type][p.metadata.get_short_name()] = p;

        fprintf(stderr, "registering model %s-%s\n", p.metadata.get_short_name().c_str(), p.type.c_str());
    }
}

//
bool PoreModelSet::has_model(const std::string& type, const std::string& short_name)
{
    PoreModelSet& model_set = getInstance();
    
    auto iter_type = model_set.model_type_sets.find(type);
    if(iter_type == model_set.model_type_sets.end()) {
        return false;
    }
    
    auto iter_short_name = iter_type->second.find(short_name);
    return iter_short_name != iter_type->second.end();
}

//
const PoreModel& PoreModelSet::get_model(const std::string& type, const std::string& short_name)
{
    PoreModelSet& model_set = getInstance();
    
    auto iter_type = model_set.model_type_sets.find(type);
    if(iter_type == model_set.model_type_sets.end()) {
        fprintf(stderr, "Error: cannot find model type %s\n", type.c_str());
        exit(EXIT_FAILURE);
    }
    
    auto iter_short_name = iter_type->second.find(short_name);
    if(iter_short_name == iter_type->second.end()) {
        fprintf(stderr, "Error: cannot find model %s for type %s\n", short_name.c_str(), type.c_str());
        exit(EXIT_FAILURE);
    }
    
    return iter_short_name->second;
}

//
const PoreModelMap& PoreModelSet::get_models(const std::string& type)
{
    PoreModelSet& model_set = getInstance();
    
    auto iter_type = model_set.model_type_sets.find(type);
    if(iter_type == model_set.model_type_sets.end()) {
        fprintf(stderr, "Error: cannot find model type %s\n", type.c_str());
        exit(EXIT_FAILURE);
    }
    
    return iter_type->second;
}

void PoreModelSet::insert_model(const std::string& type, const PoreModel& model)
{
    #pragma omp critical
    {
        PoreModelSet& model_set = getInstance();
        std::string key = model.metadata.get_short_name();
        model_set.model_type_sets[type][key] = model;
    }
}
