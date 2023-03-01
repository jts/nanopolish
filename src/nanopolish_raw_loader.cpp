//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_raw_loader - utilities and helpers for loading
// data directly from raw nanopore files without events
//
#include "nanopolish_raw_loader.h"
#include "nanopolish_profile_hmm.h"

//#define DEBUG_BANDED 1
//#define DEBUG_ADAPTIVE 1
//#define DEBUG_PRINT_STATS 1

//
SquiggleScalings estimate_scalings_using_mom(const std::string& sequence,
                                             const PoreModel& pore_model,
                                             const event_table& et)
{
    SquiggleScalings out;
    size_t k = pore_model.k;
    size_t n_kmers = sequence.size() - k + 1;
    const Alphabet* alphabet = pore_model.pmalphabet;

    // Calculate summary statistics over the events and
    // the model implied by the read
    double event_level_sum = 0.0f;
    for(size_t i = 0; i < et.n; ++i) {
        event_level_sum += et.event[i].mean;
    }

    double kmer_level_sum = 0.0f;
    double kmer_level_sq_sum = 0.0f;
    for(size_t i = 0; i < n_kmers; ++i) {
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
        double l = pore_model.get_parameters(kmer_rank).level_mean;
        kmer_level_sum += l;
        kmer_level_sq_sum += pow(l, 2.0f);
    }

    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    double event_level_sq_sum = 0.0f;
    for(size_t i = 0; i < et.n; ++i) {
        event_level_sq_sum += pow(et.event[i].mean - shift, 2.0);
    }

    double scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    out.set4(shift, scale, 0.0, 1.0);

#if DEBUG_PRINT_STATS
    fprintf(stderr, "event mean: %.2lf kmer mean: %.2lf shift: %.2lf\n", event_level_sum / et.n, kmer_level_sum / n_kmers, out.shift);
    fprintf(stderr, "event sq-mean: %.2lf kmer sq-mean: %.2lf scale: %.2lf\n", event_level_sq_sum / et.n, kmer_level_sq_sum / n_kmers, out.scale);
    fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
#endif
    return out;
}

//
bool qc_simple_event_alignment(SquiggleRead& read,
                               const PoreModel& pore_model,
                               const std::string& sequence,
                               const AdaBandedParameters parameters,
                               const std::vector<AlignedPair>& alignment, const char* short_name)
{
    // Calculate QC metrics
    size_t strand_idx = 0;
    size_t k = pore_model.k;
    const Alphabet* alphabet = pore_model.pmalphabet;
    size_t n_aligned_events = 0;
    double sum_emission = 0.0f;

    for(size_t i = 0; i < alignment.size(); ++i) {
        size_t event_idx = alignment[i].read_pos;
        size_t kmer_idx = alignment[i].ref_pos;
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(kmer_idx, k).c_str(), k);
        sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
        n_aligned_events += 1;
    }

    double avg_log_emission = sum_emission / n_aligned_events;
    bool spanned = !alignment.empty() && alignment.front().ref_pos == 0 && alignment.back().ref_pos == sequence.length() - k;
    bool passed = avg_log_emission > parameters.min_average_log_emission && spanned;
    int first_aligned_kmer = !alignment.empty() ? alignment.front().read_pos : -1;
    double events_per_kmer = 0.0f;
    if(parameters.verbose) {
        fprintf(stderr, "%s\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\n", short_name,
            read.read_name.substr(0, 6).c_str(), passed ? "OK" : "FAILED", events_per_kmer, sequence.size(), avg_log_emission, first_aligned_kmer);
    }
    return passed;
}

//
std::vector<AlignedPair> adaptive_banded_simple_event_align(SquiggleRead& read, const PoreModel& pore_model, const std::string& sequence, const AdaBandedParameters parameters)
{
    size_t strand_idx = 0;

    AdaptiveBandedViterbi abv;
    abv.initialize(read, sequence, pore_model.k, strand_idx, parameters);

    generic_banded_simple_hmm(read, pore_model, sequence, parameters, abv);
    std::vector<AlignedPair> alignment = adaptive_banded_backtrack(abv);

    // qc
    bool qc_pass = qc_simple_event_alignment(read, pore_model, sequence, parameters, alignment, abv.get_short_name());

    if(!qc_pass) {
        alignment.clear();
    }

    return alignment;
}

//
std::vector<AlignedPair> guide_banded_simple_event_align(SquiggleRead& read,
                                                                 const PoreModel& pore_model,
                                                                 const Haplotype& haplotype,
                                                                 const EventAlignmentRecord& event_align_record,
                                                                 const AdaBandedParameters parameters)
{
    size_t strand_idx = 0;
    std::vector<AlignedPair> alignment;

    EventBandedViterbi ebv;
    ebv.initialize(read, haplotype, event_align_record, pore_model.k, strand_idx, parameters);
    if(!ebv.are_bands_continuous()) {
        return alignment;
    }

    // fill in DP matrix
    generic_banded_simple_hmm(read, pore_model, haplotype.get_sequence(), parameters, ebv);
    alignment = adaptive_banded_backtrack(ebv);

    // qc
    bool qc_pass = qc_simple_event_alignment(read, pore_model, haplotype.get_sequence(), parameters, alignment, ebv.get_short_name());

    if(!qc_pass) {
        alignment.clear();
    }

    return alignment;
}

//
std::vector<EventKmerPosterior> guide_banded_simple_posterior(SquiggleRead& read,
                                                                      const PoreModel& pore_model,
                                                                      const Haplotype& haplotype,
                                                                      const EventAlignmentRecord& event_align_record,
                                                                      const AdaBandedParameters parameters)
{
    size_t strand_idx = 0;
    size_t k = pore_model.k;
    std::vector<EventKmerPosterior> assignment;
    const Alphabet* alphabet = pore_model.pmalphabet;

    // Forward
    EventBandedForward ebf;
    ebf.initialize(read, haplotype, event_align_record, pore_model.k, strand_idx, parameters);
    if(!ebf.are_bands_continuous()) {
        return assignment;
    }

    generic_banded_simple_hmm(read, pore_model, haplotype.get_sequence(), parameters, ebf);

    // Backward
    EventBandedForward ebb;
    ebb.initialize(read, haplotype, event_align_record, pore_model.k, strand_idx, parameters);
    assert(ebb.are_bands_continuous()); // if forward is OK backwards must be too
    generic_banded_simple_hmm_backwards(read, pore_model, haplotype.get_sequence(), parameters, ebb);

    float f = ebf.get_by_event_kmer(ebf.get_num_events(), ebf.get_num_kmers());
    float b = ebb.get_by_event_kmer(-1, -1);

    if(parameters.verbose) {
        fprintf(stderr, "%s\t%s\t%s\t%.2f\t%zu\t%.2f\t%.8f\t%.8f\n", "ebf",
            read.read_name.substr(0, 6).c_str(), "OK", (float)ebf.get_num_events() / ebf.get_num_kmers(), haplotype.get_sequence().length(), f / ebf.get_num_events(), f, b);
    }


    // Pass 1: calculate normalization term by event
    std::vector<float> normalization(ebf.get_num_events(), -INFINITY);
    for(size_t band_idx = 1; band_idx < ebf.get_num_bands() - 1; ++band_idx) {
        // trim start
        int trim_event_idx = band_idx - 1;
        float lp = ebf.get_by_event_kmer(trim_event_idx, -1) + ebb.get_by_event_kmer(trim_event_idx, -1) - f;
        normalization[trim_event_idx] = logsumexpf(normalization[trim_event_idx], lp);
#ifdef DEBUG_POSTERIOR_NORMALIZATION
        fprintf(stderr, "[normalization] %d %d %.8f %.8f\n", trim_event_idx, -1, lp, normalization[trim_event_idx]);
#endif
        // normal event, kmer pairs
        int min_offset, max_offset;
        ebf.get_offset_range_for_band(band_idx, min_offset, max_offset);
        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = ebf.get_event_at_band_offset(band_idx, offset);
            int kmer_idx = ebf.get_kmer_at_band_offset(band_idx, offset);
            float fke = ebf.get_by_event_kmer(event_idx, kmer_idx);
            float bke = ebb.get_by_event_kmer(event_idx, kmer_idx);
            lp = fke + bke - f;
            normalization[event_idx] = logsumexpf(normalization[event_idx], lp);
#ifdef DEBUG_POSTERIOR_NORMALIZATION
            fprintf(stderr, "[normalization] %d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", event_idx, kmer_idx, fke, bke, fke + bke, f, lp, expf(lp), normalization[event_idx]);
#endif
        }

        // trim end
        int trim_end_kmer_idx = ebf.get_num_kmers();
        lp = ebf.get_by_event_kmer(trim_event_idx, trim_end_kmer_idx) + ebb.get_by_event_kmer(trim_event_idx, trim_end_kmer_idx) - f;
        normalization[trim_event_idx] = logsumexpf(normalization[trim_event_idx], lp);
#ifdef DEBUG_POSTERIOR_NORMALIZATION
        fprintf(stderr, "[normalization] %d %d %.8f %.8f\n", trim_event_idx, ebf.get_num_kmers(), lp, normalization[trim_event_idx]);
#endif
    }

    // Pass 2: calculate posteriors
    float log_min_posterior = log(parameters.min_posterior);
    for(size_t band_idx = 1; band_idx < ebf.get_num_bands() - 1; ++band_idx) {
        int min_offset, max_offset;
        ebf.get_offset_range_for_band(band_idx, min_offset, max_offset);
        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = ebf.get_event_at_band_offset(band_idx, offset);
            int kmer_idx = ebf.get_kmer_at_band_offset(band_idx, offset);

            float fke = ebf.get_by_event_kmer(event_idx, kmer_idx);
            float bke = ebb.get_by_event_kmer(event_idx, kmer_idx);
            float n = normalization[event_idx];
            float log_w = fke + bke - f - n;
            EventKmerPosterior ekp = { event_idx, kmer_idx, log_w };
            if(log_w > log_min_posterior) {
                assignment.push_back(ekp);
                /*
                std::string kmer = haplotype.get_sequence().substr(ekp.kmer_idx, k);
                size_t kmer_rank = alphabet->kmer_rank(kmer.c_str(), k);
                fprintf(stderr, "[posterior] %d %d %s %.8f %.2f %.2f\n",
                    ekp.event_idx, ekp.kmer_idx, kmer.c_str(), exp(log_w), read.get_fully_scaled_level(event_idx, strand_idx), pore_model.states[kmer_rank].level_mean);
                //fprintf(stderr, "[posterior] \t(%.2lf %.2lf %.8lf %.8lf %.8lf)\n", fke, bke, fke + bke, n, fke + bke - n);
                */
            }
        }
    }

    return assignment;
}
