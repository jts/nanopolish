//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_pore_model_set -- A class that maintains
// a collection of pore models that SquiggleReads
// can load during initialization.
//
#ifndef NANOPOLISH_PORE_MODEL_SET_H
#define NANOPOLISH_PORE_MODEL_SET_H

#include <map>
#include "nanopolish_poremodel.h"

#define DEFAULT_MODEL_TYPE "ONT"

class PoreModelSet
{
    public:

        //
        // initialize the model set from a .fofn file
        // returns pointers to the imported models
        //
        static std::vector<const PoreModel*> initialize(const std::string& fofn_filename);

        //
        // check if a model with the same key already exists
        //
        static bool has_model(const PoreModel& model);

        //
        // check if a model with the specification exists
        //
        static bool has_model(const std::string& kit_name, 
                              const std::string& alphabet,
                              const std::string& strand,
                              size_t k);
        
        //
        // Build a descriptive key from this pore model
        //
        static std::string get_model_key(const PoreModel& model);

        //
        // get a model from the set, returning a pointer
        // this pointer will always be valid throughput the lifetime of the program
        // but it may be NULL if the requested model does not exist.
        //
        static const PoreModel* get_model(const std::string& kit_name, 
                                          const std::string& alphabet,
                                          const std::string& strand,
                                          size_t k);

        static const PoreModel* get_model_by_key(const std::string& key);

        //
        // get all the models for the combination of parameters
        //
        static std::map<std::string, PoreModel> copy_strand_models(const std::string& kit_name,
                                                                   const std::string& alphabet,
                                                                   size_t k);

        //
        // Add a copy of this model to the collection. Returns a pointer
        // owned by the PoreModelSet that the caller can use to refer to the new model
        //
        static const PoreModel* add_model(const PoreModel& p);

        // destructor
        ~PoreModelSet();

    private:

        // singleton accessor function
        static PoreModelSet& getInstance()
        {
            static PoreModelSet instance;
            return instance;
        }

        // Internal function for adding this model into the collection
        // Returns a reference to the model in the map
        PoreModel* register_model(const PoreModel& p);

        // Build a unique identify string
        static std::string get_model_key(const std::string& kit_name,
                                         const std::string& alphabet,
                                         const std::string& strand,
                                         size_t k);

        PoreModelSet();

        // do not allow copies of this classs
        PoreModelSet(PoreModelSet const&) = delete;
        void operator=(PoreModelSet const&) = delete;

        // map from a string key to a pointer to the actual model
        std::map<std::string, PoreModel*> model_map;
};

#endif
