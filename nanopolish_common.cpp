//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_common -- Data structures and definitions
// shared across files
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

// Print the alignment between the read-strand and a sequence
void print_alignment(const std::string& name,
                     uint32_t seq_id,
                     uint32_t read_id,
                     const std::string& consensus, 
                     const HMMConsReadState& state,
                     const std::vector<AlignmentState>& alignment)
{
    size_t n_matches = 0;
    size_t n_merges = 0;
    size_t n_skips = 0;
    size_t n_mergeskips = 0;
    
    char prev_s = '\0';
    for(size_t pi = 0; pi < alignment.size(); ++pi) {

        uint32_t ei = alignment[pi].event_idx;
        uint32_t ki = alignment[pi].kmer_idx;
        char s = alignment[pi].state;
    
        double level = get_drift_corrected_level(*state.read, ei, state.strand);
        double sd = state.read->events[state.strand].stdv[ei];
        double duration = get_duration(*state.read, ei, state.strand);
        uint32_t rank = get_rank(state, consensus.c_str(), ki);
        
        const PoreModel& pm = state.read->pore_model[state.strand];
        double model_m = (pm.state[rank].level_mean + pm.shift) * pm.scale;
        double model_s = pm.state[rank].level_stdv * pm.scale;
        double norm_level = (level - model_m) / model_s;
        
        double model_sd_mean = pm.state[rank].sd_mean;
        double model_sd_stdv = pm.state[rank].sd_stdv;

        n_matches += (s == 'M');
        n_merges += (s == 'E');
        n_skips += (s == 'K');
        n_mergeskips += (s == 'K' && prev_s == 'E');

        double lp_diff = 0.0f;
        if(pi > 0) {
            lp_diff = alignment[pi].l_fm - alignment[pi - 1].l_fm;
        } else {
            lp_diff = alignment[pi].l_fm;
        }
        std::string kmer = consensus.substr(ki, K);
 
        printf("DEBUG\t%s\t%d\t%d\t%c\t", name.c_str(), read_id, state.rc, state.strand ? 't' : 'c');
        printf("%c\t%d\t%d\t", s, ei, ki);
        printf("%s\t%.3lf\t", kmer.c_str(), duration);
        printf("%.1lf\t%.1lf\t%.1lf\t", level, model_m, norm_level);
        printf("\t%.1lf\t%.1lf\t%.1lf\t", sd, model_sd_mean, (sd - model_sd_mean) / model_sd_stdv);
        printf("%.2lf\t%.2lf\t%.2lf\n", exp(alignment[pi].l_posterior), alignment[pi].l_fm, lp_diff);
        prev_s = s;
    }

    // Summarize alignment
    double time_start = state.read->events[state.strand].time[state.event_start_idx];
    double time_end = state.read->events[state.strand].time[state.event_stop_idx];
    double total_duration = fabs(time_start - time_end);
    double num_events = abs(state.event_start_idx - state.event_stop_idx) + 1;
    double final_lp = alignment[alignment.size() - 1].l_fm;
    double mean_lp = final_lp / num_events;

    // Print summary header on first entry
    static int once = 1;
    if(once) {
        printf("SUMMARY\tseq_name\tseq_id\tread_id\tis_rc\tstrand\t");
        printf("lp\tmean_lp\tnum_events\t");
        printf("n_matches\tn_merges\tn_skips\tn_mergeskips\ttotal_duration\n");
        once = 0;
    }

    printf("SUMMARY\t%s\t%d\t%d\t%d\t%c\t", name.c_str(), seq_id, read_id, state.rc, state.strand ? 't' : 'c');
    printf("%.2lf\t%.2lf\t%.0lf\t", final_lp, mean_lp, num_events);
    printf("%zu\t%zu\t%zu\t%zu\t%.2lf\n", n_matches, n_merges, n_skips, n_mergeskips, total_duration);
}

