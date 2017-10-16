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
#include "nanopolish_builtin_models.h"

//
PoreModelSet::PoreModelSet()
{
    // Copy the built-in models into the map
    for(auto p : builtin_models) {
        register_model(p);
    }
}

//
PoreModelSet::~PoreModelSet()
{

}

//
std::vector<std::string> PoreModelSet::initialize(const std::string& fofn_filename)
{
    // grab singleton instance
    PoreModelSet& model_set = getInstance();

    std::vector<std::string> out;
    
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

        // add the model and push the new reference to it to the output
        const PoreModel& imported = model_set.register_model(p);
        out.push_back(model_set.get_model_key(imported));
    }
    return out;
}

PoreModel& PoreModelSet::register_model(const PoreModel& p)
{
    // Check that this model doesn't exist already
    std::string key = get_model_key(p);
    auto iter = model_map.find(key);
    if(iter != model_map.end()) {
        fprintf(stderr, "Warning: overwriting model %s\n", key.c_str());
    }
    //fprintf(stderr, "registered model with key %s\n", key.c_str());

    #pragma omp critical
    model_map[key] = p;

    return model_map[key];
}

void PoreModelSet::add_model(const PoreModel& p)
{
    PoreModelSet& model_set = getInstance();
    model_set.register_model(p);
}

bool PoreModelSet::has_model(const PoreModel& p)
{
    PoreModelSet& model_set = getInstance();
    std::string model_key = model_set.get_model_key(p);
    auto iter = model_set.model_map.find(model_key);
    return iter != model_set.model_map.end();
}

//
bool PoreModelSet::has_model(const std::string& kit_name,
                             const std::string& alphabet,
                             const std::string& strand,
                             size_t k)
{
    PoreModelSet& model_set = getInstance();
    std::string model_key = model_set.get_model_key(kit_name, alphabet, strand, k);
    auto iter = model_set.model_map.find(model_key);
    return iter != model_set.model_map.end();
}

//
const PoreModel& PoreModelSet::get_model(const std::string& kit_name,
                                         const std::string& alphabet,
                                         const std::string& strand,
                                         size_t k)
{
    PoreModelSet& model_set = getInstance();
    return get_model_by_key(model_set.get_model_key(kit_name, alphabet, strand, k));
}

const PoreModel& PoreModelSet::get_model_by_key(const std::string& key)
{
    PoreModelSet& model_set = getInstance();
    auto iter = model_set.model_map.find(key);
    if(iter == model_set.model_map.end()) {
        fprintf(stderr, "Error: cannot find model with key %s\n", key.c_str());
        exit(EXIT_FAILURE);
    }
    return iter->second;
}

//
PoreModelMap PoreModelSet::copy_strand_models(const std::string& kit_name,
                                              const std::string& alphabet,
                                              size_t k)
{
    PoreModelMap out;
    PoreModelSet& model_set = getInstance();
    for(const auto& kv : model_set.model_map) {
        const PoreModel& model = kv.second;
        if(model.metadata.get_kit_name() == kit_name &&
           model.pmalphabet->get_name() == alphabet &&
           model.k == k) {
            out.insert(kv);
        }
    }
    return out;
}

std::string PoreModelSet::get_model_key(const PoreModel& model)
{
    return PoreModelSet::get_model_key(model.metadata.get_kit_name(),
                                       model.pmalphabet->get_name(),
                                       model.metadata.get_strand_model_name(),
                                       model.k);
}

//
std::string PoreModelSet::get_model_key(const std::string& kit_name,
                                        const std::string& alphabet,
                                        const std::string& strand,
                                        size_t k)
{
    std::stringstream ss;
    ss << kit_name << "." << alphabet << "." << k << "mer." << strand;
    return ss.str();
}
