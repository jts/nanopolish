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
#define PRINT_TRAINING_MESSAGES 1

enum ProfileState
{
    PS_KMER_SKIP = 0,
    PS_EVENT_SPLIT,
    PS_MATCH,
    PS_NUM_STATES = 3
};

char ps2char(ProfileState ps) { return "KEMN"[ps]; }

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
                                         const HMMInputData& data,
                                         uint32_t ki,
                                         uint32_t kj)
{
    const PoreModel& pm = data.read->pore_model[data.strand];
    const KHMMParameters& parameters = data.read->parameters[data.strand];

    uint32_t rank_i = get_rank(data, sequence, ki);
    uint32_t rank_j = get_rank(data, sequence, kj);

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

std::vector<BlockTransitions> calculate_transitions(uint32_t num_kmers, const char* sequence, const HMMInputData& data)
{
    const KHMMParameters& parameters = data.read->parameters[data.strand];

    std::vector<BlockTransitions> transitions(num_kmers);
    
    for(uint32_t ki = 0; ki < num_kmers; ++ki) {

        // probability of skipping k_i from k_(i - 1)
        double p_skip = ki > 0 ? calculate_skip_probability(sequence, data, ki - 1, ki) : 0.0f;

        // transitions from match state in previous block
        double p_mk = p_skip;
        double p_me = (1 - p_skip) * parameters.trans_m_to_e_not_k;
        double p_mm = 1.0f - p_me - p_mk;

        // transitions from event split state in previous block
        double p_ee = parameters.trans_e_to_e;
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

    return transitions;
}

double profile_hmm_forward_fill(DoubleMatrix& fm, // forward matrix
                                const char* sequence,
                                const HMMInputData& data,
                                uint32_t e_start)
{
    PROFILE_FUNC("profile_hmm_fill_forward")

    const KHMMParameters& parameters = data.read->parameters[data.strand];

    // Calculate number of blocks
    // A block of the HMM is a set of PS_KMER_SKIP, PS_EVENT_SPLIT, PS_MATCH
    // events for one kmer
    uint32_t num_blocks = fm.n_cols / PS_NUM_STATES;
    
    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks

    std::vector<BlockTransitions> transitions = calculate_transitions(num_kmers, sequence, data);
    
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
            printf("\tTransitions: p_mx [%.3lf %.3lf %.3lf]\n", bt.lp_mm, bt.lp_me, bt.lp_mk);
            printf("\t             p_ex [%.3lf %.3lf %.3lf]\n", bt.lp_em, bt.lp_ee, 0.0f);
            printf("\t             p_lx [%.3lf %.3lf %.3lf]\n", bt.lp_km, 0.0, bt.lp_kk);

            printf("\tPS_MATCH -- Transitions: [%.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mm, bt.lp_em, bt.lp_km, 
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
                    bt.lp_me, bt.lp_ee,
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
                    bt.lp_mk, bt.lp_kk,
                    get(fm, row, prev_block_offset + PS_MATCH),
                    get(fm, row, prev_block_offset + PS_KMER_SKIP),
                    sum_k);
#endif

            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * data.stride;
            uint32_t rank = get_rank(data, sequence, kmer_idx);
            double lp_emission_m = log_probability_match(*data.read, rank, event_idx, data.strand);
            double lp_emission_e = log_probability_event_insert(*data.read, rank, event_idx, data.strand);

#ifdef DEBUG_FILL    
            printf("\tEMISSION: %.2lf\n", lp_emission_m);
#endif
            set(fm, row, curr_block_offset + PS_MATCH, sum_m + lp_emission_m);
            set(fm, row, curr_block_offset + PS_EVENT_SPLIT, sum_e + lp_emission_e);
            set(fm, row, curr_block_offset + PS_KMER_SKIP, sum_k);
        }
    }
    uint32_t last_event_row = fm.n_rows - 1;
    uint32_t last_aligned_block = num_blocks - 2;
    uint32_t match_state_last_block = PS_NUM_STATES * last_aligned_block + PS_MATCH;
    return get(fm, last_event_row, match_state_last_block);
}

double profile_hmm_score(const std::string& sequence, const HMMInputData& data)
{
    uint32_t n_kmers = sequence.size() - K + 1;

    uint32_t n_states = PS_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

    uint32_t e_start = data.event_start_idx;
    uint32_t e_end = data.event_stop_idx;
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
    double score = profile_hmm_forward_fill(fm, sequence.c_str(), data, e_start);

    // cleanup
    free_matrix(fm);
    return score;
}

void profile_hmm_viterbi_initialize(DoubleMatrix& m)
{
    // Same as forward initialization
    profile_hmm_forward_initialize(m);
}

std::vector<AlignmentState> profile_hmm_align(const std::string& sequence, const HMMInputData& data)
{
    std::vector<AlignmentState> alignment;

    uint32_t n_kmers = sequence.size() - K + 1;
    uint32_t n_states = PS_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

    uint32_t e_start = data.event_start_idx;
    uint32_t e_end = data.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;
    assert(n_events >= 2);

    uint32_t n_rows = n_events + 1;
    
    // Allocate matrices to hold the HMM result
    DoubleMatrix vm;
    allocate_matrix(vm, n_rows, n_states);
    
    UInt8Matrix bm;
    allocate_matrix(bm, n_rows, n_states);

    profile_hmm_viterbi_initialize(vm);
    profile_hmm_viterbi_fill(vm, bm, sequence.c_str(), data, e_start);

    // Traverse the backtrack matrix to compute the results
    
    // start from the last event matched to the last kmer
    uint32_t row = n_rows - 1;
    uint32_t col = PS_NUM_STATES * n_kmers + PS_MATCH;

    while(row > 0) {
        
        uint32_t event_idx = e_start + (row - 1) * data.stride;
        uint32_t block = col / PS_NUM_STATES;
        assert(block > 0);
        assert(get(vm, row, col) != -INFINITY);

        uint32_t kmer_idx = block - 1;
        
        ProfileState curr_ps = (ProfileState) (col % PS_NUM_STATES);

        AlignmentState as;
        as.event_idx = event_idx;
        as.kmer_idx = kmer_idx;
        as.l_posterior = -INFINITY; // not computed
        as.l_fm = get(vm, row, col);
        as.log_transition_probability = -INFINITY; // not computed
        as.state = ps2char(curr_ps);
        alignment.push_back(as);

        // Update the event (row) and k-mer using the current state
        // The next state is encoded in the backtrack matrix for the current cell
        ProfileState next_ps = (ProfileState)get(bm, row, col);

#if DEBUG_BACKTRACK
        printf("Backtrack [%zu %zu] k: %zu block: %zu curr_ps: %c next_ps: %c\n", row, col, kmer_idx, block, ps2char(curr_ps), ps2char(next_ps));
#endif

        if(curr_ps == PS_MATCH) {
            row -= 1;
            kmer_idx -= 1;
        } else if(curr_ps == PS_EVENT_SPLIT) {
            row -= 1;
            // kmer stays the same
        } else {
            assert(curr_ps == PS_KMER_SKIP);
            // row stays the same
            kmer_idx -= 1;
        }

        col = PS_NUM_STATES * (kmer_idx + 1) + next_ps;
    }

    //
    std::reverse(alignment.begin(), alignment.end());

    //
    free_matrix(vm);
    free_matrix(bm);

    return alignment;
}

void profile_hmm_viterbi_fill(DoubleMatrix& vm, // viterbi matrix
                              UInt8Matrix& bm, // backtrack matrix
                              const char* sequence,
                              const HMMInputData& data,
                              uint32_t e_start)
{
    PROFILE_FUNC("profile_hmm_viterbi_forward")

    const KHMMParameters& parameters = data.read->parameters[data.strand];

    // Calculate number of blocks
    // A block of the HMM is a set of PS_KMER_SKIP, PS_EVENT_SPLIT, PS_MATCH
    // events for one kmer
    uint32_t num_blocks = vm.n_cols / PS_NUM_STATES;
    
    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks

    std::vector<BlockTransitions> transitions = calculate_transitions(num_kmers, sequence, data);
    
    // Fill in matrix
    for(uint32_t row = 1; row < vm.n_rows; row++) {

        // Skip the first block which is the start state, it was initialized above
        for(uint32_t block = 1; block < num_blocks - 1; block++) {

            // retrieve transitions
            uint32_t kmer_idx = block - 1;
            BlockTransitions& bt = transitions[kmer_idx];

            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * data.stride;
            uint32_t rank = get_rank(data, sequence, kmer_idx);
            
            double lp_emission_m = log_probability_match(*data.read, rank, event_idx, data.strand);
            double lp_emission_e = log_probability_event_insert(*data.read, rank, event_idx, data.strand);
            
            uint32_t prev_block = block - 1;
            uint32_t prev_block_offset = PS_NUM_STATES * prev_block;
            uint32_t curr_block_offset = PS_NUM_STATES * block;
            
            // state PS_MATCH
            double m1 = bt.lp_mm + get(vm, row - 1, prev_block_offset + PS_MATCH);
            double m2 = bt.lp_em + get(vm, row - 1, prev_block_offset + PS_EVENT_SPLIT);
            double m3 = bt.lp_km + get(vm, row - 1, prev_block_offset + PS_KMER_SKIP);
            double max_m = std::max(std::max(m1, m2), m3);
            set(vm, row, curr_block_offset + PS_MATCH, max_m + lp_emission_m);
            
            ProfileState m_from = PS_NUM_STATES;
            if(max_m == m1)
                m_from = PS_MATCH;
            else if(max_m == m2)
                m_from = PS_EVENT_SPLIT;
            else if(max_m == m3)
                m_from = PS_KMER_SKIP;
            assert(m_from != PS_NUM_STATES);

            set(bm, row, curr_block_offset + PS_MATCH, m_from);
            assert(get(bm, row, curr_block_offset + PS_MATCH) == m_from);

            // state PS_EVENT_SPLIT
            double e1 = bt.lp_me + get(vm, row - 1, curr_block_offset + PS_MATCH);
            double e2 = bt.lp_ee + get(vm, row - 1, curr_block_offset + PS_EVENT_SPLIT);
            double max_e = std::max(e1, e2);
            ProfileState e_from = max_e == e1 ? PS_MATCH : PS_EVENT_SPLIT;

            set(vm, row, curr_block_offset + PS_EVENT_SPLIT, max_e + lp_emission_e);
            set(bm, row, curr_block_offset + PS_EVENT_SPLIT, e_from);

            // state PS_KMER_SKIP
            double k1 = bt.lp_mk + get(vm, row, prev_block_offset + PS_MATCH);
            double k2 = bt.lp_kk + get(vm, row, prev_block_offset + PS_KMER_SKIP);
            double max_k = std::max(k1, k2);
            ProfileState k_from = max_k == k1 ? PS_MATCH : PS_KMER_SKIP;
            set(vm, row, curr_block_offset + PS_KMER_SKIP, max_k);
            set(bm, row, curr_block_offset + PS_KMER_SKIP, k_from);

#if DEBUG_VITERBI
            printf("\tM: [%zu %zu] block: %zu [%.2lf %.2lf %.2lf] from: %c v: %.2lf\n", row, curr_block_offset + PS_MATCH, block, m1, m2, m3, ps2char(m_from), max_m + lp_e);
            printf("\tE: [%zu %zu] block: %zu [%.2lf %.2lf] from: %c v: %.2lf\n", row, curr_block_offset + PS_EVENT_SPLIT, block, e1, e2, ps2char(e_from), max_e + lp_e);
            printf("\tK: [%zu %zu] block: %zu [%.2lf %.2lf] from: %c v: %.2lf\n", row, curr_block_offset + PS_KMER_SKIP, block, k1, k2, ps2char(k_from), max_k + lp_e);
#endif
        }
    }
    
    /*
    uint32_t last_event_row = vm.n_rows - 1;
    uint32_t last_aligned_block = num_blocks - 2;
    uint32_t match_state_last_block = PS_NUM_STATES * last_aligned_block + PS_MATCH;

    // force match -> terminal transition
    set(bm, last_event_row, match_state_last_block, PS_MATCH);
    */
}

void profile_hmm_update_training(const std::string& consensus, 
                                 const HMMInputData& data)
{
    std::vector<AlignmentState> alignment = profile_hmm_align(consensus, data);

    const PoreModel& pm = data.read->pore_model[data.strand];
    TrainingData& training_data = data.read->parameters[data.strand].training_data;

    size_t n_kmers = consensus.size() - K + 1;
    uint32_t strand_idx = get_strand_idx(data);
    char prev_s = 'M';

    for(size_t pi = 0; pi < alignment.size(); ++pi) {

        uint32_t ei = alignment[pi].event_idx;
        uint32_t ki = alignment[pi].kmer_idx;
        char s = alignment[pi].state;
    
        // Record transition observations
        // We do not record observations for merge states as there was no kmer transitions
        // We also do not record observations for the beginning of the matches as the
        // alignment may be poor due to edge effects
        if(pi > 5 && pi < alignment.size() - 5) {
 
            // skip transition training data
            // we do not process the E state here as no k-mer move was made
            if(s != 'E') {
                uint32_t transition_kmer_from = alignment[pi - 1].kmer_idx;
                uint32_t transition_kmer_to = alignment[pi].kmer_idx;

                // Specially handle skips
                // We only want to record the first k-mer skipped if multiple were skipped
                if(s == 'K') {
                    transition_kmer_from = alignment[pi - 1].kmer_idx;
                    transition_kmer_to = transition_kmer_from + 1;
                }
                
                assert(transition_kmer_from < n_kmers && transition_kmer_to < n_kmers);

                uint32_t rank1 = get_rank(data, consensus.c_str(), transition_kmer_from);
                uint32_t rank2 = get_rank(data, consensus.c_str(), transition_kmer_to);
            
                double ke1 = (pm.state[rank1].level_mean + pm.shift) * pm.scale;
                double ke2 = (pm.state[rank2].level_mean + pm.shift) * pm.scale;

#ifdef PRINT_TRAINING_MESSAGES
                printf("TRAIN_SKIP\t%d\t%.3lf\t%.3lf\t%c\t%c\n", strand_idx, ke1, ke2, s, prev_s);
#endif
                KmerTransitionObservation to = { ke1, ke2, s };
                training_data.kmer_transitions.push_back(to);
            }

            // State-to-state transition
            add_state_transition(training_data, prev_s, s);

            // emission
            double level = get_drift_corrected_level(*data.read, ei, data.strand);
            double sd = data.read->events[data.strand].stdv[ei];
            double duration = get_duration(*data.read, ei, data.strand);
            if(ki >= n_kmers)
                printf("%zu %d %d %zu %.2lf %c\n", pi, ei, ki, n_kmers, alignment[pi].l_fm, s);
            
            assert(ki < n_kmers);
            uint32_t rank = get_rank(data, consensus.c_str(), ki);
        
            double model_m = (pm.state[rank].level_mean + pm.shift) * pm.scale;
            double model_s = pm.state[rank].level_stdv * pm.scale;
            double norm_level = (level - model_m) / model_s;

            if(s == 'M')
                training_data.emissions_for_matches.push_back(norm_level);
            prev_s = s;
#ifdef PRINT_TRAINING_MESSAGES
            printf("TRAIN_EMISSION\t%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%c\n", strand_idx, ei, level, sd, model_m, model_s, norm_level, duration, s);
#endif
        }

        // summary
        training_data.n_matches += (s == 'M');
        training_data.n_merges += (s == 'E');
        training_data.n_skips += (s == 'K');
    }
}

