//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_raw_loader - utilities and helpers for loading
// data directly from raw nanopore files without events
//
#ifndef NANOPOLISH_RAW_LOADER_H
#define NANOPOLISH_RAW_LOADER_H

#include "nanopolish_squiggle_read.h"
#include "nanopolish_anchor.h"
#include "scrappie_structures.h"
#include "nanopolish_alignment_db.h"

SquiggleScalings estimate_scalings_using_mom(const std::string& sequence,
                                             const PoreModel& pore_model,
                                             const event_table& et);

struct AdaBandedParameters
{
    double min_average_log_emission = -5.0;
    int max_gap_threshold = 50;
    int max_stay_threshold = 10000;
    int bandwidth = 100;
    float p_skip = 1e-10;
    float p_trim = 1e-2;
    int verbose = 0;
};

// Align events to k-mers of a sequence using Suzuki's adaptive banded algorithm
// see: https://www.biorxiv.org/content/early/2017/09/07/130633
std::vector<AlignedPair> adaptive_banded_simple_event_align(SquiggleRead& read,
                                                            const PoreModel& pore_model,
                                                            const std::string& sequence,
                                                            const AdaBandedParameters parameters = AdaBandedParameters());

std::vector<AlignedPair> adaptive_banded_generic_simple_event_align(SquiggleRead& read,
                                                                    const PoreModel& pore_model,
                                                                    const std::string& sequence,
                                                                    const AdaBandedParameters parameters = AdaBandedParameters());


std::vector<AlignedPair> guide_banded_generic_simple_event_align(SquiggleRead& read,
                                                                 const PoreModel& pore_model,
                                                                 const Haplotype& haplotype,
                                                                 const EventAlignmentRecord& event_align_record,
                                                                 const AdaBandedParameters parameters = AdaBandedParameters());

void guide_banded_generic_simple_posterior(SquiggleRead& read,
                                           const PoreModel& pore_model,
                                           const Haplotype& haplotype,
                                           const EventAlignmentRecord& event_align_record,
                                           const AdaBandedParameters parameters);


// Simple banded alignmend algorithm
// Deprecated, use the above
std::vector<AlignedPair> banded_simple_event_align(SquiggleRead& read,
                                                   const PoreModel& pore_model,
                                                   const std::string& sequence);


#include "nanopolish_raw_loader.inl"

#endif
