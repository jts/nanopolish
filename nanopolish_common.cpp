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
#include "nanopolish_squiggle_read.h"

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
                     const HMMInputData& data,
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
    
        double level = data.read->get_drift_corrected_level(ei, data.strand);
        double sd = data.read->events[data.strand][ei].stdv;
        double duration = data.read->get_duration(ei, data.strand);
        uint32_t rank = get_rank(data, consensus.c_str(), ki);
        
        const PoreModel& pm = data.read->pore_model[data.strand];
        GaussianParameters model = pm.get_scaled_parameters(rank);

        double norm_level = (level - model.mean) / model.stdv;
        
        double model_sd_mean = 0.0f;
        double model_sd_stdv = 0.0f;

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
 
        printf("DEBUG\t%s\t%d\t%d\t%c\t", name.c_str(), read_id, data.rc, "tc"[data.strand]);
        printf("%c\t%d\t%d\t", s, ei, ki);
        printf("%s\t%.3lf\t", kmer.c_str(), duration);
        printf("%.1lf\t%.1lf\t%.1lf\t", level, model.mean, norm_level);
        printf("\t%.1lf\t%.1lf\t%.1lf\t", sd, model_sd_mean, (sd - model_sd_mean) / model_sd_stdv);
        printf("%.2lf\t%.2lf\t%.2lf\n", exp(alignment[pi].l_posterior), alignment[pi].l_fm, lp_diff);
        prev_s = s;
    }

    // Summarize alignment
    double time_start = data.read->events[data.strand][data.event_start_idx].start_time;
    double time_end = data.read->events[data.strand][data.event_stop_idx].start_time;
    double total_duration = fabs(time_start - time_end);
    double num_events = abs(data.event_start_idx - data.event_stop_idx) + 1;
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

    printf("SUMMARY\t%s\t%d\t%d\t%d\t%c\t", name.c_str(), seq_id, read_id, data.rc, data.strand ? 't' : 'c');
    printf("%.2lf\t%.2lf\t%.0lf\t", final_lp, mean_lp, num_events);
    printf("%zu\t%zu\t%zu\t%zu\t%.2lf\n", n_matches, n_merges, n_skips, n_mergeskips, total_duration);
}

