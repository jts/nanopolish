#ifndef __TRAINING_CORE_HPP
#define __TRAINING_CORE_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "nanopolish_eventalign.h"
#include "nanopolish_squiggle_read.h"

// The state training data comes in two different
// sizes Full and Minimal. The model training functions
// only actually need the Minimal data but for exploration
// the Full data is useful so left as an option.

struct MinimalStateTrainingData
{
    //
    // Functions
    //
    MinimalStateTrainingData() = default;
    MinimalStateTrainingData(const SquiggleRead& sr,
                             const EventAlignment& ea,
                             uint32_t,
                             const std::string&,
                             const std::string&)
    {
        // scale the observation to the expected pore model
        this->level_mean = sr.get_fully_scaled_level(ea.event_idx, ea.strand_idx);
        this->log_level_mean = std::log(this->level_mean);
        this->level_stdv = sr.get_scaled_stdv(ea.event_idx, ea.strand_idx);
        this->log_level_stdv = std::log(this->level_stdv);
        this->read_var = sr.pore_model[ea.strand_idx].var;
        this->log_read_var = std::log(this->read_var);
        this->read_scale_sd = sr.pore_model[ea.strand_idx].scale_sd;
        this->log_read_scale_sd = std::log(this->read_scale_sd);
        this->read_var_sd = sr.pore_model[ea.strand_idx].var_sd;
        this->log_read_var_sd = std::log(this->read_var_sd);
    }

    static void write_header(std::ostream& os)
    {
        write_header_nonl(os);
        os << std::endl;
    }
    static void write_header_nonl(std::ostream& os)
    {
        os << "model\tmodel_kmer\tlevel_mean\tlevel_stdv\tread_var\tread_scale_sd\tread_var_sd";
    }

    void write_tsv(std::ostream& os, const std::string& model_name, const std::string& kmer) const
    {
        write_tsv_nonl(os, model_name, kmer);
        os << std::endl;
    }
    void write_tsv_nonl(std::ostream& os, const std::string& model_name, const std::string& kmer) const
    {
        os << model_name << '\t'
           << kmer << '\t'
           << std::fixed << std::setprecision(2) << level_mean << '\t'
           << level_stdv << '\t'
           << read_var << '\t'
           << read_scale_sd << '\t'
           << read_var_sd;
    }

    //
    // Data
    //
    float level_mean;
    float log_level_mean;
    float level_stdv;
    float log_level_stdv;
    float read_var;
    float log_read_var;
    float read_scale_sd;
    float log_read_scale_sd;
    float read_var_sd;
    float log_read_var_sd;
}; // struct MinimalStateTrainingData

struct FullStateTrainingData
    : public MinimalStateTrainingData
{
    //
    // Functions
    //
    FullStateTrainingData() = default;
    FullStateTrainingData(const SquiggleRead& sr,
                          const EventAlignment& ea,
                          uint32_t rank,
                          const std::string& prev_kmer,
                          const std::string& next_kmer)
        : MinimalStateTrainingData(sr, ea, rank, prev_kmer, next_kmer)
    {
        this->duration = sr.events[ea.strand_idx][ea.event_idx].duration;
        this->ref_position = ea.ref_position;
        this->ref_strand = ea.rc;
        GaussianParameters model = sr.pore_model[ea.strand_idx].get_scaled_parameters(rank);
        this->z = (sr.get_drift_corrected_level(ea.event_idx, ea.strand_idx) -  model.mean ) / model.stdv;
        this->prev_kmer = prev_kmer;
        this->next_kmer = next_kmer;
    }

    static void write_header(std::ostream& os)
    {
        write_header_nonl(os);
        os << std::endl;
    }
    static void write_header_nonl(std::ostream& os)
    {
        MinimalStateTrainingData::write_header_nonl(os);
        os << "\tduration\tref_pos\tref_strand\tz\tprev_kmer\tnext_kmer";
    }

    void write_tsv(std::ostream& os, const std::string& model_name, const std::string& kmer) const
    {
        write_tsv_nonl(os, model_name, kmer);
        os << std::endl;
    }
    void write_tsv_nonl(std::ostream& os, const std::string& model_name, const std::string& kmer) const
    {
        MinimalStateTrainingData::write_tsv_nonl(os, model_name, kmer);
        os << '\t' << duration << '\t'
           << ref_position << '\t'
           << ref_strand << '\t'
           << z << '\t'
           << prev_kmer << '\t'
           << next_kmer;
    }

    //
    // Data
    //
    float duration;
    int ref_position;
    int ref_strand;
    float z;
    std::string prev_kmer;
    std::string next_kmer;
}; // struct FullStateTrainingData

typedef MinimalStateTrainingData StateTrainingData;
//typedef FullStateTrainingData StateTrainingData;

struct ParamMixture
{
    std::vector< float > log_weights;
    std::vector< PoreModelStateParams > params;
}; // struct ParamMixture

// training functions
ParamMixture train_gaussian_mixture   (const std::vector< StateTrainingData >& data, const ParamMixture& input_mixture);
ParamMixture train_invgaussian_mixture(const std::vector< StateTrainingData >& data, const ParamMixture& input_mixture);

#endif
