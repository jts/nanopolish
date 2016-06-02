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

#define DEFAULT_MODEL_TYPE "reftrained"

typedef std::map<std::string, PoreModel> PoreModelMap;

class PoreModelSet
{
    public:

        //
        // initialize the model set from a .fofn file
        //
        static void initialize(const std::string& fofn_filename);

        //
        // get a model from the set using its type and short name
        //
        static PoreModel get_model(const std::string& type, const std::string& short_name);

        //
        // get all the models for this type
        //
        static const PoreModelMap& get_models(const std::string& type);

        //
        // insert the new model into the specified type
        //
        static void insert_model(const std::string& type, const PoreModel& model);

        // destructor
        ~PoreModelSet();

    private:

        // singleton accessor function
        static PoreModelSet& getInstance()
        {
            static PoreModelSet instance;
            return instance;
        }

        // do not allow copies of this classs
        PoreModelSet(PoreModelSet const&) = delete;
        void operator=(PoreModelSet const&) = delete;
        PoreModelSet() {}; // public constructor not allowed

        // this is a map from a pore model type (like "base" or "derived"
        // to a map of models indexed by their short name
        // for example m_model_type_sets["base"]["t.007"]
        std::map<std::string, PoreModelMap> model_type_sets;
};

#endif
