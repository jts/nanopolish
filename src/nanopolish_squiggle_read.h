//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#ifndef NANOPOLISH_SQUIGGLE_READ_H
#define NANOPOLISH_SQUIGGLE_READ_H

#include "nanopolish_common.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_read_db.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_fast5_io.h"
#include "nanopolish_fast5_loader.h"
#include <string>

enum PoreType
{
    PORETYPE_UNKNOWN = 0,
    PORETYPE_R7,
    PORETYPE_R9,
};

// the type of the read
// do not change as template must always be 0 and complement 1
enum SquiggleReadType
{
    SRT_TEMPLATE = 0,
    SRT_COMPLEMENT,
    SRT_2D
};

// the type of nucleotides that passed through the pore
enum SquiggleReadNucleotideType
{
    SRNT_DNA,
    SRNT_RNA
};

// Flags to control the behaviour of the read
enum SquiggleReadFlags
{
    SRF_NO_MODEL = 1, // do not load a model
    SRF_LOAD_RAW_SAMPLES = 2
};

// The raw event data for a read
struct SquiggleEvent
{
    float mean;        // current level mean in picoamps
    float stdv;        // current level stdv
    double start_time; // start time of the event in seconds
    float duration;    // duration of the event in seconds
    float log_stdv;    // precompute for efficiency
};

// Scaling parameters to account for per-read variations from the model
struct SquiggleScalings
{
    SquiggleScalings() : scale(1.0), shift(0.0), drift(0.0), var(1.0), scale_sd(1.0), var_sd(1.0) {}

    // Set scaling parameteters. Theses function should be used rather than
    // setting values directly so the cached values are updated.
    void set4(double _shift,
              double _scale,
              double _drift,
              double _var);

    void set6(double _shift,
              double _scale,
              double _drift,
              double _var,
              double _scale_sd,
              double _var_sd);

    // direct parameters that must be set
    double scale;
    double shift;
    double drift;
    double var;
    double scale_sd;
    double var_sd;

    // derived parameters that are cached for efficiency
    double log_var;
    double scaled_var;
    double log_scaled_var;
};

struct IndexPair
{
    IndexPair() : start(-1), stop(-1) {}
    int32_t start;
    int32_t stop; // inclusive
};

// This struct maps from base-space k-mers to the range of template/complement
// events used to call it.
struct EventRangeForBase
{
    IndexPair indices[2]; // one per strand
};

//
class SquiggleRead
{
    public:

        SquiggleRead() {} // legacy TODO remove
        SquiggleRead(const std::string& name, const ReadDB& read_db, const uint32_t flags = 0);
        SquiggleRead(const ReadDB& read_db, const Fast5Data& data, const uint32_t flags = 0);
        SquiggleRead(const std::string& sequence, const Fast5Data& data, const uint32_t flags = 0);
        ~SquiggleRead();

        //
        // I/O
        //
        void init(const std::string& read_sequence, const Fast5Data& data, const uint32_t flags);

        //
        // Access to data
        //

        // Return the duration of the specified event for one strand
        inline float get_duration(uint32_t event_idx, uint32_t strand) const
        {
            assert(event_idx < events[strand].size());
            return events[strand][event_idx].duration;
        }

        // Return the current stdv for the given event
        inline float get_stdv(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].stdv;
        }

        // Return log of the current stdv for the given event
        inline float get_log_stdv(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].log_stdv;
        }

        // Return the observed current level corrected for drift
        inline float get_drift_scaled_level(uint32_t event_idx, uint32_t strand) const
        {
            float level = get_unscaled_level(event_idx, strand);
            float time = get_time(event_idx, strand);
            return level - time * this->scalings[strand].drift;
        }

        // Return the observed current level after correcting for drift, shift and scale
        inline float get_fully_scaled_level(uint32_t event_idx, uint32_t strand) const
        {
            return (get_drift_scaled_level(event_idx, strand) - scalings[strand].shift) / scalings[strand].scale;
        }

        // Return the observed current level stdv, after correcting for scale
        inline float get_scaled_stdv(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].stdv / scalings[strand].scale_sd;
        }

        inline float get_time(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].start_time - events[strand][0].start_time;
        }

        // Return the observed current level after correcting for drift
        inline float get_unscaled_level(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].mean;
        }

        // Return k-mer sized used by the pore model for this read strand
        size_t get_model_k(uint32_t strand) const
        {
            assert(this->base_model[strand] != NULL);
            return this->base_model[strand]->k;
        }

        // Return name of the pore model kit for this read strand
        std::string get_model_kit_name(uint32_t strand) const
        {
            assert(this->base_model[strand] != NULL);
            return this->base_model[strand]->metadata.get_kit_name();
        }

        // Return name of the strand model for this read strand
        std::string get_model_strand_name(uint32_t strand) const
        {
            assert(this->base_model[strand] != NULL);
            return this->base_model[strand]->metadata.get_strand_model_name();
        }

        // Get the cached pointer to the nucleotide model for this read
        const PoreModel* get_base_model(uint32_t strand) const
        {
            assert(this->base_model[strand] != NULL);
            return this->base_model[strand];
        }

        // Get the pore model that should be used for this read, for a given alphabet
        const PoreModel* get_model(uint32_t strand, const std::string& alphabet) const
        {
            return PoreModelSet::get_model(this->get_model_kit_name(strand),
                                           alphabet,
                                           this->get_model_strand_name(strand),
                                           this->get_model_k(strand));
        }

        // Get the parameters to the gaussian PDF scaled to this read
        inline GaussianParameters get_scaled_gaussian_from_pore_model_state(const PoreModel& pore_model, size_t strand_idx, size_t rank) const
        {
            const SquiggleScalings& scalings = this->scalings[strand_idx];
            const PoreModelStateParams& params = pore_model.states[rank];
            GaussianParameters gp;
            gp.mean = scalings.scale * params.level_mean + scalings.shift;
            gp.stdv = params.level_stdv * scalings.var;
            gp.log_stdv = params.level_log_stdv + scalings.log_var;
            return gp;
        }

        // Calculate the index of this k-mer on the other strand
        inline int32_t flip_k_strand(int32_t k_idx, uint32_t k) const
        {
            assert(!read_sequence.empty());
            return read_sequence.size() - k_idx - k;
        }

        // get the index of the event that is nearest to the given kmer
        int get_closest_event_to(int k_idx, uint32_t strand) const;

        // returns true if this read has events for this strand
        bool has_events_for_strand(size_t strand_idx) const { return !this->events[strand_idx].empty(); }

        // Create an eventalignment between the events of this read and its 1D basecalled sequence
        std::vector<EventAlignment> get_eventalignment_for_1d_basecalls(const std::string& read_sequence_1d,
                                                                        const std::string& alphabet_name,
                                                                        const std::vector<EventRangeForBase>& base_to_event_map_1d,
                                                                        const size_t k,
                                                                        const size_t strand_idx,
                                                                        const int label_shift) const;

        // Sample-level access
        size_t get_sample_index_at_time(size_t sample_time) const;
        std::vector<float> get_scaled_samples_for_event(size_t strand_idx, size_t event_idx) const;
        std::pair<size_t, size_t> get_event_sample_idx(size_t strand_idx, size_t event_idx) const;

        // print the scaling parameters for this strand
        void print_scaling_parameters(FILE* fp, size_t strand_idx) const
        {
            fprintf(fp, "shift: %.2lf scale: %.2lf drift: %.2lf var: %.2lf\n", this->scalings[strand_idx].shift,
                                                                               this->scalings[strand_idx].scale,
                                                                               this->scalings[strand_idx].drift,
                                                                               this->scalings[strand_idx].var);
        }

        //
        // Data
        //

        // unique identifier of the read
        std::string read_name;
        SquiggleReadType read_type;
        SquiggleReadNucleotideType nucleotide_type;
        PoreType pore_type;
        std::string fast5_path;
        uint32_t read_id;
        std::string read_sequence;

        // one event sequence for each strand
        std::vector<SquiggleEvent> events[2];

        // scaling parameters for each strand
        SquiggleScalings scalings[2];

        // pointers to the base model (DNA or RNA) for this read
        const PoreModel* base_model[2];

        // optional fields holding the raw data
        // this is not split into strands so there is only one vector, unlike events
        std::vector<float> samples;
        double sample_rate;
        int64_t sample_start_time;
        int channel_id;

        // summary stats
        double events_per_base[2];

        //
        std::vector<EventRangeForBase> base_to_event_map;

        // one set of parameters per strand
        TransitionParameters parameters[2];

    private:

        // disallowed
        SquiggleRead(const SquiggleRead&) {}

        // Load all read data from raw samples
        void load_from_raw(const Fast5Data& fast5_data, const uint32_t flags);

        // helper for get_closest_event_to
        int get_next_event(int start, int stop, int stride, uint32_t strand) const;
};

#endif
