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

        // unique identifier of the read
        uint32_t read_id;

        // one model for each strand
        PoreModel pore_model[2];

        // one event sequence for each strand
        EventSequence events[2];

        // one set of parameters per strand
        KHMMParameters parameters[2];
};

#endif
