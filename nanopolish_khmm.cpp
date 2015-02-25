//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_khmm -- Hidden Markov Model with
// a k-mer sequence as the hidden state
//
#include <algorithm>
#include "nanopolish_khmm.h"

//
// Constants 
//
const static uint32_t KHMM_MAX_JUMP = 5;
const static uint32_t KHMM_MAX_MERGE = 10;

void khmm_fill_transitions(DoubleMatrix& matrix, const std::string& consensus, const HMMConsReadState& state)
{
    PROFILE_FUNC("fill_khmm_transitions")

    const PoreModel& pm = state.read->pore_model[state.strand];
    const KHMMParameters& parameters = state.read->parameters[state.strand];

    uint32_t n_kmers = consensus.size() - K + 1;
    uint32_t n_states = n_kmers + 2;
    uint32_t terminal_state = n_states - 1;

    assert(matrix.n_rows == n_states && matrix.n_cols == n_states);

    // Initialize the transition matrix to -INFINITY for all states
    for(size_t si = 0; si < n_states; ++si)
        for(size_t sj = 0; sj < n_states; ++sj)
            set(matrix, si, sj, -INFINITY);

    // Start state transitions -- only allowed to go to k_0
    set(matrix, 0, 1, 0.0f);
    
    // TODO: calculate in log space
    for(size_t si = 1; si < n_states - 1; si++) {
        size_t ki = si - 1;
        double sum = 0.0f;

        uint32_t last_valid_state = si + KHMM_MAX_JUMP;
        if(last_valid_state >= terminal_state)
            last_valid_state = terminal_state - 1;

        for(size_t sj = si; sj <= last_valid_state; ++sj) {
            
            size_t kj = sj - 1;
                        
            // transition probability
            double p_i_j = 0.0f;

            if(ki == kj) {
                p_i_j = parameters.self_transition;
            } else {
        
                uint32_t rank_i = get_rank(state, consensus.c_str(), ki);
                uint32_t rank_j = get_rank(state, consensus.c_str(), kj);

                double level_i = (pm.state[rank_i].level_mean + pm.shift) * pm.scale;
                double level_j = (pm.state[rank_j].level_mean + pm.shift) * pm.scale;
                
                double p_skip = get_skip_probability(parameters, level_i, level_j);
                p_i_j = (1 - sum) * (1 - p_skip);
                assert(p_i_j >= 0.0f && p_i_j <= 1.0f);

#ifdef DEBUG_TRANSITION
                printf("\t\t%zu -> %zu %.2lf %.2lf p_skip: %.4lf p: %.2lf\n", ki, kj, level_i, level_j, p_skip, p_i_j);
#endif
            }

            sum += p_i_j;
            set(matrix, si, sj, log(p_i_j));
        }
    }

    // Transition to end state -- only the last k-mer can go to the end state
    // TODO: allow the last k-mer to be skipped??
    set(matrix, n_states - 2, n_states - 1, 0.0f);
}

void khmm_forward_initialize(DoubleMatrix& fm)
{
    // initialize forward calculation
    for(uint32_t si = 0; si < fm.n_cols; si++)
        set(fm, 0, si, -INFINITY);
    for(uint32_t ri = 0; ri < fm.n_rows; ri++)
        set(fm, ri, 0, -INFINITY);

    set(fm, 0, 0, 0.0f); // probability 1 in the start state for the null row
}

// Terminate the forward algorithm by calculating
// the probability of transitioning to the end state
// for all columns and a given row
double khmm_forward_terminate(const DoubleMatrix& fm,
                              const DoubleMatrix& tm,
                              uint32_t row)
{
    double sum = -INFINITY;
    uint32_t tcol = fm.n_cols - 1;
    for(uint32_t sk = 0; sk < fm.n_cols - 1; sk++) {
        // transition probability from state k to state l
        double t_kl = get(tm, sk, tcol);
        double fm_k = get(fm, row, sk);
        sum = add_logs(sum, t_kl + fm_k);
    }
    return sum;
}

double khmm_forward_fill(DoubleMatrix& fm, // forward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMConsReadState& state,
                         uint32_t e_start)
{
    PROFILE_FUNC("fill_forward_khmm")

    // Fill in matrix
    for(uint32_t row = 1; row < fm.n_rows; row++) {
        for(uint32_t sl = 1; sl < fm.n_cols - 1; sl++) {

            // cell indices
            //uint32_t c = cell(matrix, row, col);
            //uint32_t event_i = e_start + (row - 1) * state.stride;
            //uint32_t kmer_idx = k_start + col - 1;

            // Sum over states for previous row
            double sum = -INFINITY;

            // Only look back as far as the first state that can jump here
            uint32_t first_possible_state = 0;
            if(sl >= KHMM_MAX_JUMP)
                first_possible_state = sl - KHMM_MAX_JUMP;

            for(uint32_t sk = first_possible_state; sk <= sl; sk++) {

                // transition probability from state k to state l
                double t_kl = get(tm, sk, sl);
                double fm_k = get(fm, row - 1, sk);
                sum = add_logs(sum, t_kl + fm_k);
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf f: %.2lf s: %.2lf\n", row, sl, sk, t_kl, fm_k, sum);
#endif
            }

            // Emission probability for event i in state sl
            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t kmer_idx = sl - 1;
            uint32_t rank = get_rank(state, sequence, kmer_idx);
            double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);
            
            set(fm, row, sl, lp_e + sum);

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) ei: %zu ki: %zu\n", row, sl, event_idx, kmer_idx);
            printf("(%d %d) sum: %.2lf lp_e: %.2lf fm: %.2lf\n", row, sl, sum, lp_e, get(fm, row, sl));
#endif
        }
    }

    // terminate by summing the last row and transitioning to end state
    double sum = -INFINITY;
    uint32_t tcol = fm.n_cols - 1;
    uint32_t lrow = fm.n_rows - 1;
    for(uint32_t sk = 0; sk < fm.n_cols - 1; sk++) {

        // transition probability from state k to state l
        double t_kl = get(tm, sk, tcol);
        double fm_k = get(fm, lrow, sk);
        sum = add_logs(sum, t_kl + fm_k);
    }
    return sum;
}

void khmm_backward_initialize(DoubleMatrix& bm, const DoubleMatrix& tm)
{
    // initialize forward calculation
    uint32_t tcol = tm.n_cols - 1;
    uint32_t row = bm.n_rows - 1;

    for(uint32_t si = 0; si < bm.n_cols; si++)
        set(bm, row, si, get(tm, si, tcol));
}

void khmm_backward_fill(DoubleMatrix& bm, // backward matrix
                        const DoubleMatrix& tm, //transitions
                        const char* sequence,
                        const HMMConsReadState& state,
                        uint32_t e_start)
{
    // Fill in matrix
    for(uint32_t row = bm.n_rows - 2; row > 0; row--) {
        for(uint32_t sk = 1; sk < bm.n_cols - 1; sk++) {

            // Sum over states for next row
            double sum = -INFINITY;
            for(uint32_t sl = 1; sl < bm.n_cols - 1; sl++) {

                // transition probability from state k to state l
                double t_kl = get(tm, sk, sl);
                double bm_l = get(bm, row + 1, sl);

                // Emit E_(i+1) in state sl
                uint32_t event_idx = e_start + row * state.stride; // for i + 1
                uint32_t kmer_idx = sl - 1;
                uint32_t rank = get_rank(state, sequence, kmer_idx);
                double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);

                sum = add_logs(sum, lp_e + t_kl + bm_l);
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf b: %.2lf e: %.2lf s: %.2lf\n", row, sk, sl, t_kl, bm_l, lp_e, sum);
#endif
            }
            
            set(bm, row, sk, sum);

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) bm: %.2lf\n", row, sk, get(bm, row, sk));
#endif
        }
    }
}

double khmm_score(const std::string& consensus, const HMMConsReadState& state, AlignmentPolicy policy)
{
    uint32_t n_kmers = consensus.size() - K + 1;
    uint32_t n_states = n_kmers + 2; // one start and one end state

    DoubleMatrix tm;
    allocate_matrix(tm, n_states, n_states);

    khmm_fill_transitions(tm, consensus, state);
    
    uint32_t e_start = state.event_start_idx;
    uint32_t e_end = state.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate a matrix to hold the HMM result
    DoubleMatrix fm;
    allocate_matrix(fm, n_rows, n_states);

    khmm_forward_initialize(fm);
    khmm_forward_fill(fm, tm, consensus.c_str(), state, e_start);

    double score = 0.0f;
    if(policy == AP_GLOBAL) {
        // score by the bottom-right cell
        uint32_t last_row = fm.n_rows - 1;
        score = khmm_forward_terminate(fm, tm, last_row);
    } else if(policy == AP_SEMI_KMER) {

        // score by the best cell in the last column
        double best_score = -INFINITY;
        uint32_t best_row = 0;
        for(size_t row = 1; row < fm.n_rows - 1; ++row) {
            double s = khmm_forward_terminate(fm, tm, row);
            if(s > best_score) {
                best_score = s;
                best_row = row;
            }
        }

        score = best_score;
    } else {
        assert(false);
    }

    free_matrix(tm);
    free_matrix(fm);
    return score;
}


std::vector<AlignmentState> khmm_posterior_decode(const std::string& sequence, const HMMConsReadState& state)
{
    uint32_t n_kmers = sequence.size() - K + 1;
    uint32_t n_states = n_kmers + 2; // one start and one end state

    DoubleMatrix tm;
    allocate_matrix(tm, n_states, n_states);

    khmm_fill_transitions(tm, sequence, state);
    
    uint32_t e_start = state.event_start_idx;
    uint32_t e_end = state.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate and compute forward matrix
    DoubleMatrix fm;
    allocate_matrix(fm, n_rows, n_states);

    khmm_forward_initialize(fm);
    double lf = khmm_forward_fill(fm, tm, sequence.c_str(), state, e_start);

    // Allocate and compute backward matrix
    DoubleMatrix bm;
    allocate_matrix(bm, n_rows, n_states);

    khmm_backward_initialize(bm, tm);
    khmm_backward_fill(bm, tm, sequence.c_str(), state, e_start);

    // posterior decode
    std::vector<AlignmentState> output;
    
    uint32_t row = fm.n_rows - 1;
    uint32_t col = fm.n_cols - 1;

    while(row > 0) {

        // Calculate posterior probability that e_i is matched to k_j
        double max_posterior = -INFINITY;
        uint32_t max_s = 0;

        // Only check states that are possible to transition given the previous match col
        uint32_t first_possible_col = 1;
        if(col >= KHMM_MAX_JUMP)
            first_possible_col = col - KHMM_MAX_JUMP;
        
        for(uint32_t si = first_possible_col; si <= col; ++si) {
            double lp = get(fm, row, si) + get(bm, row, si) - lf;
            if(lp > max_posterior) {
                max_posterior = lp;
                max_s = si;
            }
        }
    
        uint32_t event_idx = e_start + (row - 1) * state.stride;
        uint32_t kmer_idx = max_s - 1;
 
        double lpfm = get(fm, row, max_s);

        AlignmentState ps = { event_idx, kmer_idx, max_posterior, lpfm, 0.0f, 'N' };
        output.push_back(ps);

        //
        row -= 1;
        col = max_s;
    }

    std::reverse(output.begin(), output.end());

    // First state is always a match
    output[0].state = 'M';
    uint32_t prev_ei = output[0].event_idx;
    uint32_t prev_ki = output[0].kmer_idx;

    // store the transition probability to this state
    // the + 1 is to convert a k-mer index to a column
    output[0].log_transition_probability = get(tm, 0, output[0].kmer_idx + 1);

    for(uint32_t pi = 1; pi < output.size(); ++pi) {
        uint32_t ei = output[pi].event_idx;
        uint32_t ki = output[pi].kmer_idx;

        output[pi].log_transition_probability = get(tm, prev_ki + 1, ki + 1);
        assert(abs(ei - prev_ei) == 1);

        if(ki == prev_ki) {
            output[pi].state = 'E';
        } else if(ki - prev_ki == 1) {
            output[pi].state = 'M';
        } else {
            assert(ki - prev_ki > 1);
            output[pi].state = 'K';
        }

        prev_ei = ei;
        prev_ki = ki;
    }

    free_matrix(tm);
    free_matrix(fm);
    free_matrix(bm);
    return output;
}

void khmm_update_training(const std::string& consensus, 
                          const HMMConsReadState& state)
{
    std::vector<AlignmentState> pstates = khmm_posterior_decode(consensus, state);

    const PoreModel& pm = state.read->pore_model[state.strand];
    TrainingData& training_data = state.read->parameters[state.strand].training_data;
    size_t n_kmers = consensus.size() - K + 1;
    uint32_t strand_idx = get_strand_idx(state);

    for(size_t pi = 0; pi < pstates.size(); ++pi) {

        uint32_t ei = pstates[pi].event_idx;
        uint32_t ki = pstates[pi].kmer_idx;
        char s = pstates[pi].state;
    
        // Record transition observations
        // We do not record observations for merge states as there was no kmer transitions
        // We also do not record observations for the beginning of the matches as the
        // alignment may be poor due to edge effects
        if(pi > 5 && pi < pstates.size() - 5) {
 
            // transition           
            if(s != 'E') {
                uint32_t transition_kmer_from = pstates[pi - 1].kmer_idx;
                uint32_t transition_kmer_to = pstates[pi].kmer_idx;

                // Specially handle skips
                // We only want to record the first k-mer skipped if multiple were skipped
                if(s == 'K') {
                    transition_kmer_from = pstates[pi - 1].kmer_idx;
                    transition_kmer_to = transition_kmer_from + 1;
                }
                
                assert(transition_kmer_from < n_kmers && transition_kmer_to < n_kmers);

                uint32_t rank1 = get_rank(state, consensus.c_str(), transition_kmer_from);
                uint32_t rank2 = get_rank(state, consensus.c_str(), transition_kmer_to);
            
                double ke1 = (pm.state[rank1].level_mean + pm.shift) * pm.scale;
                double ke2 = (pm.state[rank2].level_mean + pm.shift) * pm.scale;

#ifdef PRINT_TRAINING_MESSAGES
                printf("TRAIN_SKIP\t%d\t%.3lf\t%.3lf\t%c\n", strand_idx, ke1, ke2, s);
#endif
                KmerTransitionObservation to = { ke1, ke2, s };
                training_data.kmer_transitions.push_back(to);
            }

            // emission
            double level = get_drift_corrected_level(*state.read, ei, state.strand);
            double sd = state.read->events[state.strand].stdv[ei];
            double start_time = state.read->events[state.strand].time[ei];
            double end_time = state.read->events[state.strand].time[ei + 1];
            if(ki >= n_kmers)
                printf("%zu %d %d %zu %.2lf %c\n", pi, ei, ki, n_kmers, pstates[pi].l_fm, s);
            
            assert(ki < n_kmers);
            uint32_t rank = get_rank(state, consensus.c_str(), ki);
        
            double model_m = (pm.state[rank].level_mean + pm.shift) * pm.scale;
            double model_s = pm.state[rank].level_stdv * pm.scale;
            double norm_level = (level - model_m) / model_s;

            if(s == 'M')
                training_data.emissions_for_matches.push_back(norm_level);

#ifdef PRINT_TRAINING_MESSAGES
            printf("TRAIN_EMISSION\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%c\n", strand_idx, ei, level, sd, model_m, model_s, norm_level, end_time - start_time, s);
#endif
        }

        // summary
        training_data.n_matches += (s == 'M');
        training_data.n_merges += (s == 'E');
        training_data.n_skips += (s == 'K');
    }
}
