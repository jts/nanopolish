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

class PoreModelSet
{
    public:

        //
        // get a model from the set using its name
        // usage: PoreModelSet::get_model_by_name(name);
        //
        static PoreModel get_model_by_name(const std::string& name);

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

        std::map<std::string, PoreModel> m_pore_model_cache;
};

#endif
