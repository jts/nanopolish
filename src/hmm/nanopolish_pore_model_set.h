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

typedef std::map<std::string, PoreModel> PoreModelMap;

class PoreModelSet
{
    public:

        //
        // initialize the model set from a .fofn file
        // returns the keys of the imported models
        //
        static std::vector<std::string> initialize(const std::string& fofn_filename);

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
        // get a model from the set
        //
        static const PoreModel& get_model(const std::string& kit_name, 
                                          const std::string& alphabet,
                                          const std::string& strand,
                                          size_t k);

        static const PoreModel& get_model_by_key(const std::string& key);

        //
        // get all the models for the combination of parameters
        //
        static PoreModelMap copy_strand_models(const std::string& kit_name,
                                               const std::string& alphabet,
                                               size_t k);

        //
        // 
        //
        static void add_model(const PoreModel& p);

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
        PoreModel& register_model(const PoreModel& p);

        // Build a unique identify string
        static std::string get_model_key(const std::string& kit_name,
                                         const std::string& alphabet,
                                         const std::string& strand,
                                         size_t k);

        PoreModelSet();

        // do not allow copies of this classs
        PoreModelSet(PoreModelSet const&) = delete;
        void operator=(PoreModelSet const&) = delete;

        // map from a string representing a pore model to the actual model
        std::map<std::string, PoreModel> model_map;
};

#endif
