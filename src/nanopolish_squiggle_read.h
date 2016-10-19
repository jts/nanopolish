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
#include <string>

enum PoreType
{
    PT_UNKNOWN = 0,
    PT_R7,
    PT_R9,
};

// the type of the read
// do not change as template must always be 0 and complement 1
enum SquiggleReadType
{
    SRT_TEMPLATE = 0,
    SRT_COMPLEMENT,
    SRT_2D
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
    float mean;       // current level mean in picoamps
    float stdv;       // current level stdv
    double start_time; // start time of the event in seconds
    float duration;     // duration of the event in seconds
    float log_stdv;   // precompute for efficiency
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

        SquiggleRead() : drift_correction_performed(false) {} // legacy TODO remove
        SquiggleRead(const std::string& name, const std::string& path, const uint32_t flags = 0);
        ~SquiggleRead();

        //
        // I/O
        // 

        //
        // Access to data
        //

        // Return the duration of the specified event for one strand
        inline float get_duration(uint32_t event_idx, uint32_t strand) const
        {
            assert(event_idx < events[strand].size());
            return events[strand][event_idx].duration;
        }

        // Return the observed current level after correcting for drift
        inline float get_drift_corrected_level(uint32_t event_idx, uint32_t strand) const
        {
            assert(drift_correction_performed);
            return events[strand][event_idx].mean;
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

        // Return the observed current level after correcting for drift, shift and scale
        inline float get_fully_scaled_level(uint32_t event_idx, uint32_t strand) const
        {
            assert(drift_correction_performed);
            float level = get_drift_corrected_level(event_idx, strand);
            return (level - pore_model[strand].shift) / pore_model[strand].scale;
        }

        // Return the observed current level stdv, after correcting for scale
        inline float get_scaled_stdv(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].stdv / pore_model[strand].scale_sd;
        }

        inline float get_time(uint32_t event_idx, uint32_t strand) const
        {
            return events[strand][event_idx].start_time - events[strand][0].start_time;
        }

        // Return the observed current level after correcting for drift
        inline float get_uncorrected_level(uint32_t event_idx, uint32_t strand) const
        {
            if (!drift_correction_performed)
                return events[strand][event_idx].mean;
            else {
                double time = get_time(event_idx, strand);
                return events[strand][event_idx].mean + (time * pore_model[strand].drift);
            }
        }
        
        // Calculate the index of this k-mer on the other strand
        inline int32_t flip_k_strand(int32_t k_idx) const
        {
            assert(!read_sequence.empty());
            return read_sequence.size() - k_idx - pore_model[T_IDX].k;
        }

        // Transform each event by correcting for current drift
        void transform();

        // get the index of the event that is nearest to the given kmer 
        int get_closest_event_to(int k_idx, uint32_t strand) const;

        // replace the pore models with the models specified in the map or by a string
        void replace_models(const std::string& model_type);
        void replace_model(size_t strand_idx, const std::string& model_type);
        void replace_model(size_t strand_idx, const PoreModel& model);

        // returns true if this read has events for this strand
        bool has_events_for_strand(size_t strand_idx) { return !this->events[strand_idx].empty(); }

        // Create an eventalignment between the events of this read and its 1D basecalled sequence
        std::vector<EventAlignment> get_eventalignment_for_1d_basecalls(const std::string& read_sequence_1d,
                                                                        const std::vector<EventRangeForBase>& base_to_event_map_1d, 
                                                                        const size_t k, 
                                                                        const size_t strand_idx) const;

        // Sample-level access
        size_t get_sample_index_at_time(size_t sample_time) const;
        std::vector<float> get_scaled_samples_for_event(size_t strand_idx, size_t event_idx) const;

        // print the scaling parameters for this strand
        void print_scaling_parameters(FILE* fp, size_t strand_idx) const
        {
            fprintf(fp, "shift: %.2lf scale: %.2lf drift: %.2lf var: %.2lf\n", this->pore_model[strand_idx].shift,
                                                                               this->pore_model[strand_idx].scale,
                                                                               this->pore_model[strand_idx].drift,
                                                                               this->pore_model[strand_idx].var);
        }

        //
        // Data
        //

        // unique identifier of the read
        std::string read_name;
        SquiggleReadType read_type;
        PoreType pore_type;
        std::string fast5_path;
        uint32_t read_id;
        std::string read_sequence;
        bool drift_correction_performed;

        // one model for each strand
        PoreModel pore_model[2];

        // one event sequence for each strand
        std::vector<SquiggleEvent> events[2];
        
        // optional fields holding the raw data
        // this is not split into strands so there is only one vector, unlike events
        std::vector<float> samples;
        double sample_rate;
        int64_t sample_start_time;

        //
        std::vector<EventRangeForBase> base_to_event_map;

        // one set of parameters per strand
        TransitionParameters parameters[2];

    private:
        // private data
        fast5::File* f_p;
        std::string basecall_group;

        // warning triggers
        static bool& show_warning_read_name_not_extract() { static bool val = true; return val; }

        SquiggleRead(const SquiggleRead&) {}

        // Load all the read data from a fast5 file
        void load_from_fast5(const uint32_t flags);

        // Version-specific intialization functions
        void _load_R7(uint32_t si);
        void _load_R9(uint32_t si,
                      const std::string& read_sequence_1d,
                      const std::vector<EventRangeForBase>& event_map_1d,
                      const std::vector<double>& p_model_states,
                      const uint32_t flags);

        // make a map from a base of the 1D read sequence to the range of events supporting that base
        std::vector<EventRangeForBase> build_event_map_1d(const std::string& read_sequence_1d,
                                                          uint32_t strand, 
                                                          std::vector<fast5::Event_Entry>& f5_events);

        // as above but for the 2D sequence. this fills in both the template and complete event indices
        void build_event_map_2d_r7();
        void build_event_map_2d_r9();

        // helper for get_closest_event_to
        int get_next_event(int start, int stop, int stride, uint32_t strand) const;
        // detect pore_type
        void detect_pore_type();
        // detect basecall_group and read_type
        void detect_basecall_group();
        // check basecall_group and read_type
        bool check_basecall_group() const;
};

#endif
