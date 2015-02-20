//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_common -- Data structures and definitions
// shared across files
//
#include "nanopolish_common.h"

//
double get_duration(const SquiggleRead& read, uint32_t event_idx, uint32_t strand)
{
    double e_start = read.events[strand].time[event_idx];
    double e_end = read.events[strand].time[event_idx + 1];
    return e_end - e_start;
}

//
double get_drift_corrected_level(const SquiggleRead& read, uint32_t event_idx, uint32_t strand)
{
    double level = read.events[strand].level[event_idx];
    // correct level by drift
    double start = read.events[strand].time[0];
    double time = read.events[strand].time[event_idx] - start;
    return level - (time * read.pore_model[strand].drift);
}

// Increment the input string to be the next DNA sequence in lexicographic order
void lexicographic_next(std::string& str)
{
    int carry = 1;
    int i = str.size() - 1;
    do {
        uint32_t r = base_rank[str[i]] + carry;
        str[i] = "ACGT"[r % 4];
        carry = r / 4;
        i -= 1;
    } while(carry > 0 && i >= 0);
}
