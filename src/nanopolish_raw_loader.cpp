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

//
void estimate_scalings_using_mom(const std::string& sequence,
                                 const PoreModel& pore_model,
                                 const event_table& et,
                                 double& out_shift,
                                 double& out_scale)
{
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
    out_shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    double event_level_sq_sum = 0.0f;
    for(size_t i = 0; i < et.n; ++i) {
        event_level_sq_sum += pow(et.event[i].mean - out_shift, 2.0);
    }

    out_scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

#if DEBUG_PRINT_STATS
    fprintf(stderr, "event mean: %.2lf kmer mean: %.2lf shift: %.2lf\n", event_level_sum / et.n, kmer_level_sum / n_kmers, out_shift);
    fprintf(stderr, "event sq-mean: %.2lf kmer sq-mean: %.2lf scale: %.2lf\n", event_level_sq_sum / et.n, kmer_level_sq_sum / n_kmers, out_scale);
    fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
#endif
}

std::vector<AlignedPair> banded_simple_event_align(SquiggleRead& read, const std::string& sequence)
{
    size_t strand_idx = 0;
    size_t k = read.pore_model[strand_idx].k;
    const Alphabet* alphabet = read.pore_model[strand_idx].pmalphabet;

#if DEBUG_PRINT_STATS
    // Build a debug event->kmer map
    std::vector<size_t> kmer_for_event(read.events[strand_idx].size());
    for(size_t ki = 0; ki < read.base_to_event_map.size(); ++ki) {
        IndexPair& elem = read.base_to_event_map[ki].indices[0];
        if(elem.start == -1) {
            continue;
        }

        for(int j = elem.start; j <= elem.stop; ++j) {
            if(j >= 0 && j < kmer_for_event.size()) {
                kmer_for_event[j] = ki;
            }
        }
    }
#endif

    const uint8_t FROM_D = 0;
    const uint8_t FROM_U = 1;
    const uint8_t FROM_L = 2;
    
    // qc
    double min_average_log_emission = -5.0;

    // banding
    int bandwidth = 500;
    int half_band = bandwidth / 2;

    // transitions
    double lp_skip = log(0.001);
    double lp_stay = log(0.5);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.1);

    size_t n_events = read.events[strand_idx].size();
    size_t n_kmers = sequence.size() - k + 1;

    // Calculate the minimum event index that is within the band for each read kmer
    // We determine this using the expected number of events observed per kmer
    double events_per_kmer = (double)n_events / n_kmers;
    std::vector<int> min_event_idx_by_kmer(n_kmers);
    for(size_t ki = 0; ki < n_kmers; ++ki) {
        int expected_event_idx = (double)(ki * events_per_kmer);
        min_event_idx_by_kmer[ki] = std::max(expected_event_idx - half_band, 0);
    }
    // Initialize DP matrices
    DoubleMatrix viterbi_matrix;
    UInt8Matrix backtrack_matrix;

    size_t n_rows = bandwidth;
    size_t n_cols = n_kmers + 1;
    allocate_matrix(viterbi_matrix, n_rows, n_cols);
    allocate_matrix(backtrack_matrix, n_rows, n_cols);

    for(size_t i = 0; i < bandwidth; ++i) {
        set(viterbi_matrix, i, 0, i * lp_trim);
        set(backtrack_matrix, i, 0, 0);
    }

    // Fill in the matrix
    for(int col = 1; col < n_cols; ++col) {
        int kmer_idx = col - 1;
        int min_event_idx = min_event_idx_by_kmer[kmer_idx];
        int min_event_idx_prev_col = kmer_idx > 0 ? min_event_idx_by_kmer[kmer_idx - 1] : 0;
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(kmer_idx, k).c_str(), k);

        for(int row = 0; row < n_rows; ++row) {
            
            int event_idx = min_event_idx + row;
            if(event_idx >= n_events) {
                set(viterbi_matrix, row, col, -INFINITY);
                continue;
            }

            // dp update
            // here we are calculating whether the event for each neighboring cell is within the band
            // and calculating its position within the column
            int row_up = event_idx - min_event_idx - 1;
            int row_diag = event_idx - min_event_idx_prev_col - 1;
            int row_left = event_idx - min_event_idx_prev_col;

            double up = row_up >= 0 && row_up < n_rows ?        get(viterbi_matrix, row_up, col) : -INFINITY;
            double diag = row_diag >= 0 && row_diag < n_rows ?  get(viterbi_matrix, row_diag, col - 1) : -INFINITY;
            double left = row_left >= 0 && row_left < n_rows ?  get(viterbi_matrix, row_left, col - 1) : -INFINITY;
            
            float lp_emission = log_probability_match_r9(read, kmer_rank, event_idx, strand_idx);
     
            double score_d = diag + lp_step + lp_emission;
            double score_u = up + lp_stay + lp_emission;
            double score_l = left + lp_skip;

            double max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;
     
            set(viterbi_matrix, row, col, max_score);
            set(backtrack_matrix, row, col, from);
        }
    }

    // Backtrack

    // Initialize by finding best alignment between an event and the last kmer
    int curr_k_idx = n_kmers - 1;
    int curr_event_idx = 0;
    double max_score = -INFINITY;
    for(size_t row = 0; row < n_rows; ++row) {
        size_t col = curr_k_idx + 1;
        int ei = row + min_event_idx_by_kmer[curr_k_idx];
        double s = get(viterbi_matrix, row, col) + (n_events - ei - 1) * lp_trim;
        if(s > max_score && ei < n_events) {
            max_score = s;
            curr_event_idx = ei;
        }
    }

    // debug stats
    double sum_emission = 0;
    double n_aligned_events = 0;
    std::vector<AlignedPair> out;

#if DEBUG_PRINT_STATS
    int sum_distance_from_debug = 0;
    int max_distance_from_debug = 0;
    int num_exact_matches = 0;
    int max_distance_from_expected = 0;
    int sum_distance_from_expected = 0;
    std::vector<int> debug_event_for_kmer(n_kmers, -1);
#endif

    while(curr_k_idx >= 0) {
        // emit alignment
        out.push_back({curr_k_idx, curr_event_idx});
        
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(curr_k_idx, k).c_str(), k);
        sum_emission += log_probability_match_r9(read, kmer_rank, curr_event_idx, strand_idx);
        n_aligned_events += 1;

#if DEBUG_PRINT_STATS
        // update debug stats
        debug_event_for_kmer[curr_k_idx] = curr_event_idx;
        
        int kd = abs(curr_k_idx - kmer_for_event[curr_event_idx]);
        sum_distance_from_debug += kd;
        max_distance_from_debug = std::max(kd, max_distance_from_debug);
        num_exact_matches += curr_k_idx == kmer_for_event[curr_event_idx];

        int expected_event = (double)(curr_k_idx * events_per_kmer);
        int ed = curr_event_idx - expected_event;
        max_distance_from_expected = std::max(abs(ed), max_distance_from_expected);
        sum_distance_from_expected += ed;
#endif

        // update indices using backtrack pointers
        int row = curr_event_idx - min_event_idx_by_kmer[curr_k_idx];
        int col = curr_k_idx + 1;

        uint8_t from = get(backtrack_matrix, row, col);
        if(from == FROM_D) {
            curr_k_idx -= 1;
            curr_event_idx -= 1;
        } else if(from == FROM_U) {
            curr_event_idx -= 1;
        } else {
            curr_k_idx -= 1;
        }   
    }
    std::reverse(out.begin(), out.end());

#if DEBUG_PRINT_MATRIX
    // Columm-wise debugging
    for(size_t ci = 1; ci < n_cols; ++ci) {
        size_t curr_k_idx = ci - 1;
        size_t expected_event = (double)(curr_k_idx * events_per_kmer);

        // Find the maximum value in this column
        double max_row = -INFINITY;
        size_t max_event_idx = 0;

        std::stringstream debug_out;
        for(size_t ri = 0; ri < n_rows; ++ri) {
            size_t curr_event_idx = ri + min_event_idx_by_kmer[curr_k_idx];
            double s = get(viterbi_matrix, ri, ci);
            if(s > max_row) {
                max_row = s;
                max_event_idx = curr_event_idx;
            }
            debug_out << s << ",";
        }

        std::string debug_str = debug_out.str();
        fprintf(stderr, "DEBUG_MATRIX\t%s\n", debug_str.substr(0, debug_str.size() - 1).c_str());

        size_t aligned_event = debug_event_for_kmer[curr_k_idx];
        fprintf(stderr, "DEBUG BAND: k: %zu ee: %zu me: %zu ae: %zu d: %d\n", 
            curr_k_idx, expected_event, max_event_idx, aligned_event, (int)max_event_idx - (int)aligned_event); 
    }
#endif

    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    bool spanned = out.front().ref_pos == 0 &&  out.back().ref_pos == n_kmers - 1;
    
    if(avg_log_emission < min_average_log_emission || !spanned) {
        out.clear();
    }
    
#if DEBUG_PRINT_STATS
    fprintf(stderr, "events per base: %.2lf\n", events_per_kmer);
    fprintf(stderr, "truth stats -- avg: %.2lf max: %d exact: %d\n", (double)sum_distance_from_debug / n_kmers, max_distance_from_debug, num_exact_matches);
    fprintf(stderr, "event stats -- avg: %.2lf max: %d\n", (double)sum_distance_from_expected / n_events, max_distance_from_expected);
    fprintf(stderr, "emission stats -- avg: %.2lf\n", sum_emission / n_aligned_events);
#endif
    free_matrix(viterbi_matrix);
    free_matrix(backtrack_matrix);
    return out;
}
