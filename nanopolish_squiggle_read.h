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

#include "nanopolish_poremodel.h"
#include "nanopolish_khmm_parameters.h"
#include <string>

// The raw event data for a read
struct EventSequence
{
    uint32_t n_events;
    const double* level;
    const double* stdv;
    const double* start;
    const double* duration;
};

//
class SquiggleRead
{
    public:

        SquiggleRead() {} // legacy TODO remove
        SquiggleRead(const std::string& name, const std::string& path);

        //
        // I/O
        // 

        void load_from_fast5(const std::string& fast5_path);

        //
        // Access to data
        //

        // Return the duration of the specified event for one strand
        inline double get_duration(uint32_t event_idx, uint32_t strand) const 
        {
            assert(event_idx < events[strand].n_events);
            return events[strand].duration[event_idx];
        }

        // Return the observed current level after correcting for drift
        inline double get_drift_corrected_level(uint32_t event_idx, uint32_t strand) const
        {
            double level = events[strand].level[event_idx];
            // correct level by drift
            double read_start = events[strand].start[0];
            double time = events[strand].start[event_idx] - read_start;
            return level - (time * pore_model[strand].drift);
        }

        //
        // Data
        //

        // unique identifier of the read
        std::string read_name;
        uint32_t read_id;
        std::string twod_sequence;

        // one model for each strand
        PoreModel pore_model[2];

        // one event sequence for each strand
        EventSequence events[2];

        // one set of parameters per strand
        KHMMParameters parameters[2];
};

#endif
