//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
//
#include <algorithm>
#include "nanopolish_profile_hmm.h"

//#define DEBUG_FILL

enum ProfileState
{
    PS_KMER_SKIP = 0,
    PS_EVENT_SPLIT,
    PS_MATCH,
    PS_NUM_STATES = 3
};

void profile_hmm_forward_initialize(DoubleMatrix& fm)
{
    // initialize forward calculation
    for(uint32_t si = 0; si < fm.n_cols; si++) {
        set(fm, 0, si, -INFINITY);
    }

    for(uint32_t ri = 0; ri < fm.n_rows; ri++) {
        set(fm, ri, PS_KMER_SKIP, -INFINITY);
        set(fm, ri, PS_EVENT_SPLIT, -INFINITY);
        set(fm, ri, PS_MATCH, -INFINITY);
    }

    set(fm, 0, PS_KMER_SKIP, -INFINITY);
    set(fm, 0, PS_EVENT_SPLIT, -INFINITY);
    set(fm, 0, PS_MATCH, 0.0f);
}

// Terminate the forward algorithm by calculating
// the probability of transitioning to the end state
// for all columns and a given row
double profile_hmm_forward_terminate(const DoubleMatrix& fm,
                                     const DoubleMatrix& tm,
                                     uint32_t row)
{
    assert(false);
    return -INFINITY;

    /*
    double sum = -INFINITY;
    uint32_t tcol = fm.n_cols - 1;
    for(uint32_t sk = 0; sk < fm.n_cols - 1; sk++) {
        // transition probability from state k to state l
        double t_kl = get(tm, sk, tcol);
        double fm_k = get(fm, row, sk);
        sum = add_logs(sum, t_kl + fm_k);
    }
    return sum;
    */
}

inline double calculate_skip_probability(const char* sequence, 
                                         const HMMConsReadState& state, 
                                         uint32_t ki, 
                                         uint32_t kj)
{
    const PoreModel& pm = state.read->pore_model[state.strand];
    const KHMMParameters& parameters = state.read->parameters[state.strand];

    uint32_t rank_i = get_rank(state, sequence, ki);
    uint32_t rank_j = get_rank(state, sequence, kj);

    double level_i = (pm.state[rank_i].level_mean + pm.shift) * pm.scale;
    double level_j = (pm.state[rank_j].level_mean + pm.shift) * pm.scale;
    
    return get_skip_probability(parameters, level_i, level_j);
}

// Pre-computed transitions from the previous block
// into the current block of states. Log-scaled.
struct BlockTransitions
{
    // Transition from m state
    double lp_me;
    double lp_mk;
    double lp_mm;

    // Transitions from e state
    double lp_ee;
    double lp_em;
    
    // Transitions from e state
    double lp_kk;
    double lp_km;
};

double profile_hmm_forward_fill(DoubleMatrix& fm, // forward matrix
                                const char* sequence,
                                const HMMConsReadState& state,
                                uint32_t e_start)
{
    PROFILE_FUNC("profile_hmm_fill_forward")

    const KHMMParameters& parameters = state.read->parameters[state.strand];

    // Calculate number of blocks
    // A block of the HMM is a set of PS_KMER_SKIP, PS_EVENT_SPLIT, PS_MATCH
    // events for one kmer
    uint32_t num_blocks = fm.n_cols / PS_NUM_STATES;
    
    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks

    std::vector<BlockTransitions> transitions(num_kmers);
    
    for(uint32_t ki = 0; ki < num_kmers; ++ki) {

        // probability of skipping k_i from k_(i - 1)
        double p_skip = ki > 0 ? calculate_skip_probability(sequence, state, ki - 1, ki) : 0.0f;

        // transitions from match state in previous block
        double p_me = parameters.self_transition;
        double p_mk = (1.0f - p_me) * p_skip;
        double p_mm = 1.0f - p_me - p_mk;

        // transitions from event split state in previous block
        double p_ee = parameters.self_transition;
        double p_em = 1.0f - p_ee;
        // p_ie not allowed

        // transitions from kmer skip state in previous block
        double p_kk = p_skip;
        double p_km = 1 - p_skip;
        // p_ei not allowed    

        // log-transform and store
        BlockTransitions& bt = transitions[ki];

        bt.lp_me = log(p_me);
        bt.lp_mk = log(p_mk);
        bt.lp_mm = log(p_mm);

        bt.lp_ee = log(p_ee);
        bt.lp_em = log(p_em);
        
        bt.lp_kk = log(p_kk);
        bt.lp_km = log(p_km);
    }

    // Fill in matrix
    for(uint32_t row = 1; row < fm.n_rows; row++) {

        // Skip the first block which is the start state, it was initialized above
        // Similarily skip the last block, which is calculated in the terminate() function
        for(uint32_t block = 1; block < num_blocks - 1; block++) {

            // retrieve transitions
            uint32_t kmer_idx = block - 1;
            BlockTransitions& bt = transitions[kmer_idx];

            uint32_t prev_block = block - 1;
            uint32_t prev_block_offset = PS_NUM_STATES * prev_block;
            uint32_t curr_block_offset = PS_NUM_STATES * block;
            
            // state PS_MATCH
            double m1 = bt.lp_mm + get(fm, row - 1, prev_block_offset + PS_MATCH);
            double m2 = bt.lp_em + get(fm, row - 1, prev_block_offset + PS_EVENT_SPLIT);
            double m3 = bt.lp_km + get(fm, row - 1, prev_block_offset + PS_KMER_SKIP);
            double sum_m = add_logs(add_logs(m1, m2), m3);

#ifdef DEBUG_FILL    
            printf("Row %zu block %zu\n", row, block);
            printf("\tTransitions: p_mx [%.3lf %.3lf %.3lf]\n", p_mm, p_me, p_mk);
            printf("\t             p_ex [%.3lf %.3lf %.3lf]\n", p_em, p_ee, 0.0f);
            printf("\t             p_lx [%.3lf %.3lf %.3lf]\n", p_km, 0.0, p_kk);

            printf("\tPS_MATCH -- Transitions: [%.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf] sum: %.2lf\n", 
                    p_mm, p_em, p_km, 
                    get(fm, row - 1, prev_block_offset + PS_MATCH),
                    get(fm, row - 1, prev_block_offset + PS_EVENT_SPLIT),
                    get(fm, row - 1, prev_block_offset + PS_KMER_SKIP),
                    sum_m);
#endif


            // state PS_EVENT_SPLIT
            double e1 = bt.lp_me + get(fm, row - 1, curr_block_offset + PS_MATCH);
            double e2 = bt.lp_ee + get(fm, row - 1, curr_block_offset + PS_EVENT_SPLIT);
            
            double sum_e = add_logs(e1, e2);
#ifdef DEBUG_FILL    
            printf("\tPS_EVENT_SPLIT -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    p_me, p_ee,
                    get(fm, row - 1, curr_block_offset + PS_MATCH),
                    get(fm, row - 1, curr_block_offset + PS_EVENT_SPLIT),
                    sum_e);
#endif


            // state PS_KMER_SKIP
            double k1 = bt.lp_mk + get(fm, row, prev_block_offset + PS_MATCH);
            double k2 = bt.lp_kk + get(fm, row, prev_block_offset + PS_KMER_SKIP);
            double sum_k = add_logs(k1, k2);

#ifdef DEBUG_FILL    
            printf("\tPS_KMER_SKIP -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    p_mk, p_kk,
                    get(fm, row, prev_block_offset + PS_MATCH),
                    get(fm, row, prev_block_offset + PS_KMER_SKIP),
                    sum_k);
#endif

            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t rank = get_rank(state, sequence, kmer_idx);
            double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);

#ifdef DEBUG_FILL    
            printf("\tEMISSION: %.2lf\n", lp_e);
#endif
            set(fm, row, curr_block_offset + PS_MATCH, sum_m + lp_e);
            set(fm, row, curr_block_offset + PS_EVENT_SPLIT, sum_e + lp_e);
            set(fm, row, curr_block_offset + PS_KMER_SKIP, sum_k);
        }
    }

    uint32_t last_aligned_block = num_blocks - 2;
    return get(fm, fm.n_rows - 2, PS_NUM_STATES * last_aligned_block + PS_MATCH);
}

double profile_hmm_score(const std::string& sequence, const HMMConsReadState& state)
{
    uint32_t n_kmers = sequence.size() - K + 1;

    uint32_t n_states = PS_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

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

    profile_hmm_forward_initialize(fm);
    double score = profile_hmm_forward_fill(fm, sequence.c_str(), state, e_start);

    // cleanup
    free_matrix(fm);
    return score;
}

