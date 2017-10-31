//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#ifndef NANOPOLISH_POREMODEL_H
#define NANOPOLISH_POREMODEL_H

#include <assert.h>
#include "nanopolish_common.h"
#include <inttypes.h>
#include <string>
#include <map>
#include "nanopolish_model_names.h"
#include <fast5.hpp>

//
struct PoreModelStateParams
{
    //
    // Data
    //

    // primary data, loaded from a file, fast5 or built-in model
    double level_mean;
    double level_stdv;
    double sd_mean;
    double sd_stdv;

    // data calculated from the above, for convenience
    double level_log_stdv;
    double sd_lambda;
    double sd_log_lambda;

    //
    // Constructors, initializors, etc
    //
    PoreModelStateParams() {};

    PoreModelStateParams(double lm, double ls, double sm, double ss)
    {
        level_mean = lm;
        level_stdv = ls;
        sd_mean = sm;
        sd_stdv = ss;
        update_sd_lambda();
        update_logs();
    }

    PoreModelStateParams& operator =(const fast5::Basecall_Model_State& e)
    {
        level_mean = e.level_mean;
        level_stdv = e.level_stdv;
        sd_mean = e.sd_mean;
        sd_stdv = e.sd_stdv;
        update_logs();
        update_sd_lambda();
        return *this;
    }

    // update sd_lambda based on sd_mean & sd_stdv
    void update_sd_lambda()
    {
        sd_lambda = pow(sd_mean, 3.0) / pow(sd_stdv, 2.0);
    }
    // update sd_stdv based on sd_mean & sd_lambda
    void update_sd_stdv()
    {
        sd_stdv = pow(pow(sd_mean, 3.0) / sd_lambda, .5);
    }
    void update_logs()
    {
        level_log_stdv = log(level_stdv);
        sd_log_lambda = log(sd_lambda);
    }
};

//
class PoreModel
{
    public:
        PoreModel(uint32_t _k=5) : k(_k), pmalphabet(&gDNAAlphabet) {}

        // These constructors and the output routine take an alphabet 
        // so that kmers are inserted/written in order
        // nicer might be to store the states as a map from kmer -> state
        PoreModel(const std::string filename, const Alphabet *alphabet=NULL);
        PoreModel(fast5::File *f_p,
                  const size_t strand,
                  const std::string& bc_gr = std::string(),
                  const Alphabet *alphabet=NULL);

        void write(const std::string filename, const std::string modelname="") const;

        inline PoreModelStateParams get_parameters(const uint32_t kmer_rank) const
        {
            return states[kmer_rank];
        }
        
        inline size_t get_num_states() const { return states.size(); }

        // update states with those given, or from another model
        void update_states( const PoreModel &other );
        void update_states( const std::vector<PoreModelStateParams> &otherstates );

        // Set the metadata from a kit/strand string
        void set_metadata(const std::string& kit, const std::string& strand);

        //
        // Data
        //

        // model metadata
        std::string model_filename;
        std::string name;
        std::string type;
        ModelMetadata metadata;
        uint32_t k;
        const Alphabet *pmalphabet; 

        // model parameters, one per k-mer
        std::vector<PoreModelStateParams> states;
};

#endif
