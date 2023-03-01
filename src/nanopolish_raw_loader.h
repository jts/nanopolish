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

#include <limits>
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
    float min_posterior = 1e-3;
    int verbose = 0;
};

// Align events of a read to k-mers of some sequence (either a basecalled read or the reference)
// There are two variations on this algorithm:
//    -Suzuki's adaptive banded algorithm (https://www.biorxiv.org/content/early/2017/09/07/130633)
//     changes the position of the band along the anti-diagonal based on the scores in the corners of the band
//    -the "guide banded" algorithm uses an existing alignment to determine the placement of the band.
//
std::vector<AlignedPair> adaptive_banded_simple_event_align(SquiggleRead& read,
                                                                    const PoreModel& pore_model,
                                                                    const std::string& sequence,
                                                                    const AdaBandedParameters parameters = AdaBandedParameters());

//
std::vector<AlignedPair> guide_banded_simple_event_align(SquiggleRead& read,
                                                                 const PoreModel& pore_model,
                                                                 const Haplotype& haplotype,
                                                                 const EventAlignmentRecord& event_align_record,
                                                                 const AdaBandedParameters parameters = AdaBandedParameters());

//
// Calculate the posterior probability of an event being generate from some k-mer state

//
struct EventKmerPosterior
{
    int event_idx;
    int kmer_idx;
    float log_posterior;
};

//
std::vector<EventKmerPosterior> guide_banded_simple_posterior(SquiggleRead& read,
                                                                      const PoreModel& pore_model,
                                                                      const Haplotype& haplotype,
                                                                      const EventAlignmentRecord& event_align_record,
                                                                      const AdaBandedParameters parameters);

#include "nanopolish_raw_loader.inl"

#endif
