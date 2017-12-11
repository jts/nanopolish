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

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

std::vector<AlignedPair> adaptive_banded_simple_event_align(SquiggleRead& read, const PoreModel& pore_model, const std::string& sequence)
{
    size_t strand_idx = 0;
    size_t k = pore_model.k;
    const Alphabet* alphabet = pore_model.pmalphabet;
    size_t n_events = read.events[strand_idx].size();
    size_t n_kmers = sequence.size() - k + 1;

    // backtrack markers
    const uint8_t FROM_D = 0;
    const uint8_t FROM_U = 1;
    const uint8_t FROM_L = 2;
 
    // qc
    double min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int bandwidth = 100;
    int half_bandwidth = bandwidth / 2;
 
    // transition penalties
    double events_per_kmer = (double)n_events / n_kmers;
    double p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);
 
    // dp matrix
    size_t n_rows = n_events + 1;
    size_t n_cols = n_kmers + 1;
    size_t n_bands = n_rows + n_cols;
 
    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
    std::vector<size_t> kmer_ranks(n_kmers);
    for(size_t i = 0; i < n_kmers; ++i) {
        kmer_ranks[i] = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
    }

    typedef std::vector<float> bandscore;
    typedef std::vector<uint8_t> bandtrace;

    std::vector<bandscore> bands(n_bands);
    std::vector<bandtrace> trace(n_bands);

    for(size_t i = 0; i < n_bands; ++i) {
        bands[i].resize(bandwidth, -INFINITY);
        trace[i].resize(bandwidth, 0);
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two bands have their coordinates initialized, the rest are computed adaptively
    struct EventKmerPair
    {
        int event_idx;
        int kmer_idx;
    };

    std::vector<EventKmerPair> band_lower_left(n_bands);
 
    // initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    // band 0: score zero in the central cell
    int start_cell_offset = band_kmer_to_offset(0, -1);
    assert(is_offset_valid(start_cell_offset));
    assert(band_event_to_offset(0, -1) == start_cell_offset);
    bands[0][start_cell_offset] = 0.0f;
    
    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0);
    assert(kmer_at_offset(1, first_trim_offset) == -1);
    assert(is_offset_valid(first_trim_offset));
    bands[1][first_trim_offset] = lp_trim;
    trace[1][first_trim_offset] = FROM_U;

    int fills = 0;
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, bands[1][first_trim_offset]);
#endif

    // fill in remaining bands
    for(int band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = bands[band_idx - 1][0];
        float ur = bands[band_idx - 1][bandwidth - 1];
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;
        
        bool right = false;
        if(ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if(right) {
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
        }

/*
        float max_score = -INFINITY;
        int tmp_max_offset = 0;
        for(int tmp = 0; tmp < bandwidth; ++tmp) {
            float s = bands[band_idx - 1][tmp];
            if(s > max_score) {
                max_score = s;
                tmp_max_offset = tmp;
            }
        }
        fprintf(stderr, "bi: %d ll: %.2f up: %.2f [%d %d] [%d %d] max: %.2f [%d %d] move: %s\n", 
            band_idx, bands[band_idx - 1][0], bands[band_idx - 1][bandwidth - 1], 
            band_lower_left[band_idx - 1].event_idx, band_lower_left[band_idx - 1].kmer_idx,
            event_at_offset(band_idx - 1, bandwidth - 1), kmer_at_offset(band_idx - 1, bandwidth - 1),
            max_score, event_at_offset(band_idx - 1, tmp_max_offset), kmer_at_offset(band_idx - 1, tmp_max_offset),
            (right ? "RIGHT" : "DOWN"));
*/

        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1);
        if(is_offset_valid(trim_offset)) {
            int event_idx = event_at_offset(band_idx, trim_offset);
            if(event_idx >= 0 && event_idx < n_events) {
                bands[band_idx][trim_offset] = lp_trim * (event_idx + 1);
                trace[band_idx][trim_offset] = FROM_U;
            } else {
                bands[band_idx][trim_offset] = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int event_max_offset = band_event_to_offset(band_idx, -1);

        int min_offset = std::max(kmer_min_offset, event_min_offset);
        min_offset = std::max(min_offset, 0);

        int max_offset = std::min(kmer_max_offset, event_max_offset);
        max_offset = std::min(max_offset, bandwidth);

        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset);
            int kmer_idx = kmer_at_offset(band_idx, offset);

            size_t kmer_rank = kmer_ranks[kmer_idx];
 
            int offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1); 
            int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);
 
#ifdef DEBUG_ADAPTIVE
            // verify loop conditions
            assert(kmer_idx >= 0 && kmer_idx < n_kmers);
            assert(event_idx >= 0 && event_idx < n_events);
            assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
            assert(offset_up - offset_left == 1);
            assert(offset >= 0 && offset < bandwidth);
#endif

            float up   = is_offset_valid(offset_up)   ? bands[band_idx - 1][offset_up]   : -INFINITY;
            float left = is_offset_valid(offset_left) ? bands[band_idx - 1][offset_left] : -INFINITY;
            float diag = is_offset_valid(offset_diag) ? bands[band_idx - 2][offset_diag] : -INFINITY;
 
            float lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
            float score_d = diag + lp_step + lp_emission;
            float score_u = up + lp_stay + lp_emission;
            float score_l = left + lp_skip;

            float max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;
    
#ifdef DEBUG_ADAPTIVE
            fprintf(stderr, "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left);
            fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
            fprintf(stderr, "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", band_idx, offset, event_idx, kmer_idx, max_score, from, lp_emission);
#endif
            bands[band_idx][offset] = max_score;
            trace[band_idx][offset] = from;
            fills += 1;
        }
    }
    
    /*
    // Debug, print some of the score matrix
    for(int col = 0; col <= 10; ++col) {
        for(int row = 0; row < 100; ++row) {
            int kmer_idx = col - 1;
            int event_idx = row - 1;
            int band_idx = event_kmer_to_band(event_idx, kmer_idx);
            int offset = band_kmer_to_offset(band_idx, kmer_idx);
            assert(offset == band_event_to_offset(band_idx, event_idx));
            assert(event_idx == event_at_offset(band_idx, offset));
            fprintf(stdout, "ei: %d ki: %d bi: %d o: %d s: %.2f\n", event_idx, kmer_idx, band_idx, offset, bands[band_idx][offset]);
        }
    }
    */

    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;
    std::vector<AlignedPair> out;
    
    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = n_kmers -1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for(int event_idx = 0; event_idx < n_events; ++event_idx) {
        int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
        assert(band_idx < bands.size());
        int offset = band_event_to_offset(band_idx, event_idx);
        if(is_offset_valid(offset)) {
            float s = bands[band_idx][offset] + (n_events - event_idx) * lp_trim;
            if(s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }
    
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score);
#endif

    int curr_gap = 0;
    int max_gap = 0;
    while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        
        // emit alignment
        out.push_back({curr_kmer_idx, curr_event_idx});
#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx);
#endif       
        // qc stats
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(curr_kmer_idx, k).c_str(), k);
        sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, curr_event_idx, strand_idx);
        n_aligned_events += 1;

        int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        uint8_t from = trace[band_idx][offset];
        if(from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if(from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = std::max(curr_gap, max_gap);
        }   
    }
    std::reverse(out.begin(), out.end());
    
    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    
    bool failed = false;
    if(avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        failed = true;
        out.clear();
    }

    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    return out;
}

std::vector<AlignedPair> banded_simple_event_align(SquiggleRead& read, const PoreModel& pore_model, const std::string& sequence)
{
    size_t strand_idx = 0;
    size_t k = pore_model.k;
    const Alphabet* alphabet = pore_model.pmalphabet;

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
    int bandwidth = 1000;
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

    for(size_t i = 0; i < n_cols; ++i) {
        set(viterbi_matrix, 0, i, -INFINITY);
        set(backtrack_matrix, 0, i, 0);
    }

    for(size_t i = 0; i < bandwidth; ++i) {
        set(viterbi_matrix, i, 0, i * lp_trim);
        set(backtrack_matrix, i, 0, 0);
    }

    // Fill in the matrix
    int fills = 0;
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
            
            float lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
     
            double score_d = diag + lp_step + lp_emission;
            double score_u = up + lp_stay + lp_emission;
            double score_l = left + (kmer_idx > 0 ? lp_skip : lp_step + lp_emission);

            double max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;
     
            //fprintf(stderr, "[orgfill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
#ifdef DEBUG_BANDED
            fprintf(stderr, "[orgfill] e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", event_idx, kmer_idx, max_score, from, lp_emission);
#endif
            set(viterbi_matrix, row, col, max_score);
            set(backtrack_matrix, row, col, from);
            fills += 1;
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
#ifdef DEBUG_BANDED
    fprintf(stderr, "[orgback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_k_idx, max_score);
#endif

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
#ifdef DEBUG_BANDED
        fprintf(stderr, "[orgback] ei: %d ki: %d\n", curr_event_idx, curr_k_idx);
#endif
        
        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(curr_k_idx, k).c_str(), k);
        sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, curr_event_idx, strand_idx);
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
    bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    
    bool failed = false;
    if(avg_log_emission < min_average_log_emission || !spanned) {
        failed = true;
        out.clear();
    }
    //fprintf(stderr, "org\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, fills);
    
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
