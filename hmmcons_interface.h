// TODO: Boilerplate

#ifndef HMMCONS_INTERFACE_H
#define HMMCONS_INTERFACE_H

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
    double* mean;
    double* sd;
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

#endif

