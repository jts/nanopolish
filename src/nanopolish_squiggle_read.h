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
#include "nanopolish_khmm_parameters.h"
#include <string>

// The raw event data for a read
struct SquiggleEvent
{
    float mean;       // current level mean in picoamps
    float stdv;       // current level stdv
    float start_time; // start time of the event in seconds
    float duration;     // duration of the event in seconds
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
        SquiggleRead(const std::string& name, const std::string& path);
        ~SquiggleRead();

        //
        // I/O
        // 

        // Load all the read data from a fast5 file
        void load_from_fast5(const std::string& fast5_path);

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
        
        // Calculate the index of this k-mer on the other strand
        inline int32_t flip_k_strand(int32_t k_idx) const
        {
            assert(!read_sequence.empty());
            return read_sequence.size() - k_idx - pore_model[T_IDX].k;
        }

        // Transform each event by correcting for current drift
        void transform();

        // get the index of the event tht is nearest to the given kmer 
        int get_closest_event_to(int k_idx, uint32_t strand) const;

        //
        // Data
        //

        // unique identifier of the read
        std::string read_name;
        uint32_t read_id;
        std::string read_sequence;
        bool drift_correction_performed;

        // one model for each strand
        PoreModel pore_model[2];

        // one event sequence for each strand
        std::vector<SquiggleEvent> events[2];

        //
        std::vector<EventRangeForBase> base_to_event_map;

        // one set of parameters per strand
        KHMMParameters parameters[2];

    private:
        
        SquiggleRead(const SquiggleRead& other) {}

        // helper for get_closest_event_to
        int get_next_event(int start, int stop, int stride, uint32_t strand) const;
};

#endif
