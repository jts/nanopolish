//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_interface -- data structures for the 
// C++/python interface
//
#ifndef NANOPOLISH_INTERFACE_H
#define NANOPOLISH_INTERFACE_H

// These structures are initialized by python
// and passed into the compiled code to prepare
// the data structures. These definitions can't
// be modified without changing the corresponding
// python code.

// Pore Model
struct CPoreModelInterface
{
    // Number of k-mers that are modelled
    int n_states;
    
    //
    double* level_mean;
    double* level_stdv;
    double* sd_mean;
    double* sd_stdv;
   
    // transformation parameters
    double scale;
    double shift;
    double drift;
    double var;
};

// Events
struct CEventSequenceInterface
{
    int n_events;
    double* level;
    double* stdv;
    double* time;
};

// SquiggleRead
struct CSquiggleReadInterface
{
    CPoreModelInterface pore_model[2];
    CEventSequenceInterface events[2];
};

// ReadState
struct CReadStateInterface
{
    int read_idx;
    int event_start_idx;
    int event_stop_idx;
    int strand;
    int stride;
    int rc;
};

// ReadAnchor
struct CReadAnchorInterface
{
    int event_idx;
    int rc;
};

#endif

