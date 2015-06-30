//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
//
inline float calculate_skip_probability(const char* sequence,
                                        const HMMInputData& data,
                                        uint32_t ki,
                                        uint32_t kj)
{
    const PoreModel& pm = data.read->pore_model[data.strand];
    const KHMMParameters& parameters = data.read->parameters[data.strand];

    uint32_t rank_i = get_rank(data, sequence, ki);
    uint32_t rank_j = get_rank(data, sequence, kj);

    GaussianParameters level_i = pm.get_scaled_parameters(rank_i);
    GaussianParameters level_j = pm.get_scaled_parameters(rank_j);

    return get_skip_probability(parameters, level_i.mean, level_j.mean);
}

inline std::vector<BlockTransitions> calculate_transitions(uint32_t num_kmers, const char* sequence, const HMMInputData& data)
{
    const KHMMParameters& parameters = data.read->parameters[data.strand];

    std::vector<BlockTransitions> transitions(num_kmers);
    
    for(uint32_t ki = 0; ki < num_kmers; ++ki) {

        // probability of skipping k_i from k_(i - 1)
        float p_skip = ki > 0 ? calculate_skip_probability(sequence, data, ki - 1, ki) : 0.0f;

        // transitions from match state in previous block
        float p_mk = p_skip;
        float p_me = (1 - p_skip) * parameters.trans_m_to_e_not_k;
        float p_mm = 1.0f - p_me - p_mk;

        // transitions from event split state in previous block
        float p_ee = parameters.trans_e_to_e;
        float p_em = 1.0f - p_ee;
        // p_ie not allowed

        // transitions from kmer skip state in previous block
        float p_kk = p_skip;
        float p_km = 1 - p_skip;
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

// Output writer for the Forward Algorithm
class ProfileHMMForwardOutput
{
    public:
        ProfileHMMForwardOutput(FloatMatrix* p) : p_fm(p) {}
        
        //
        inline void update_m(uint32_t row, uint32_t col, float m, float e, float k, float lp_emission)
        {
            float sum = add_logs(m, add_logs(e, k)) + lp_emission;
            set(*p_fm, row, col, sum);
        }

        //
        inline void update_e(uint32_t row, uint32_t col, float m, float e, float lp_emission)
        {
            float sum = add_logs(m, e) + lp_emission;
            set(*p_fm, row, col, sum);
        }
        
        //
        inline void update_k(uint32_t row, uint32_t col, float m, float k, float lp_emission)
        {
            float sum = add_logs(m, k) + lp_emission;
            set(*p_fm, row, col, sum);
        }

        // get the log probability stored at a particular row/column
        inline float get(uint32_t row, uint32_t col) const
        {
            return ::get(*p_fm, row, col);
        }

        inline size_t get_num_columns() const
        {
            return p_fm->n_cols;
        }

        inline size_t get_num_rows() const
        {
            return p_fm->n_rows;
        }
    
    private:
        ProfileHMMForwardOutput(); // not allowed
        FloatMatrix* p_fm;
};

// Output writer for the Viterbi Algorithm
class ProfileHMMViterbiOutput
{
    public:
        ProfileHMMViterbiOutput(FloatMatrix* pf, UInt8Matrix* pb) : p_fm(pf), p_bm(pb) {}
        
        inline void update_m(uint32_t row, uint32_t col, float m, float e, float k, float lp_emission)
        {
            // probability update
            float max = std::max(m, std::max(e, k));
            set(*p_fm, row, col, max + lp_emission);

            // backtrack update
            uint8_t from;
            if(max == m)
                from = PS_MATCH;
            else if(max == e)
                from = PS_EVENT_SPLIT;
            else if(max == k)
                from = PS_KMER_SKIP;
            set(*p_bm, row, col, from);
        }

        inline void update_e(uint32_t row, uint32_t col, float m, float e, float lp_emission)
        {
            // probability update
            float max = std::max(m, e);
            set(*p_fm, row, col, max + lp_emission);
            
            // backtrack update
            ProfileState from = max == m ? PS_MATCH : PS_EVENT_SPLIT;
            set(*p_bm, row, col, from);
        }

        inline void update_k(uint32_t row, uint32_t col, float m, float k, float lp_emission)
        {
            // probability update
            float max = std::max(m, k);
            set(*p_fm, row, col, max + lp_emission);
            
            // backtrack update
            ProfileState from = max == m ? PS_MATCH : PS_KMER_SKIP;
            set(*p_bm, row, col, from);
        }

        // get the log probability stored at a particular row/column
        inline float get(uint32_t row, uint32_t col) const
        {
            return ::get(*p_fm, row, col);
        }

        inline size_t get_num_columns() const
        {
            return p_fm->n_cols;
        }

        inline size_t get_num_rows() const
        {
            return p_fm->n_rows;
        }
    
    private:
        ProfileHMMViterbiOutput(); // not allowed

        FloatMatrix* p_fm;
        UInt8Matrix* p_bm;
};

// This function fills in a matrix with the result of running the HMM.
// The templated ProfileHMMOutput class allows one to run either Viterbi
// or the Forward algorithm.
template<class ProfileHMMOutput>
inline float profile_hmm_fill_generic(const char* sequence,
                                      const HMMInputData& data,
                                      const uint32_t e_start,
                                      ProfileHMMOutput& output)
{
    PROFILE_FUNC("profile_hmm_fill_generic")

    const KHMMParameters& parameters = data.read->parameters[data.strand];

    // Calculate number of blocks
    // A block of the HMM is a set of PS_KMER_SKIP, PS_EVENT_SPLIT, PS_MATCH
    // events for one kmer
    uint32_t num_blocks = output.get_num_columns() / PS_NUM_STATES;

    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks

    std::vector<BlockTransitions> transitions = calculate_transitions(num_kmers, sequence, data);
 
    // Precompute kmer ranks
    std::vector<uint32_t> kmer_ranks(num_kmers);
    for(size_t ki = 0; ki < num_kmers; ++ki)
        kmer_ranks[ki] = get_rank(data, sequence, ki);

    // Fill in matrix
    for(uint32_t row = 1; row < output.get_num_rows(); row++) {

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
            float m1 = bt.lp_mm + output.get(row - 1, prev_block_offset + PS_MATCH);
            float m2 = bt.lp_em + output.get(row - 1, prev_block_offset + PS_EVENT_SPLIT);
            float m3 = bt.lp_km + output.get(row - 1, prev_block_offset + PS_KMER_SKIP);

            // state PS_EVENT_SPLIT
            float e1 = bt.lp_me + output.get(row - 1, curr_block_offset + PS_MATCH);
            float e2 = bt.lp_ee + output.get(row - 1, curr_block_offset + PS_EVENT_SPLIT);

            // state PS_KMER_SKIP
            float k1 = bt.lp_mk + output.get(row, prev_block_offset + PS_MATCH);
            float k2 = bt.lp_kk + output.get(row, prev_block_offset + PS_KMER_SKIP);

            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * data.event_stride;
            uint32_t rank = kmer_ranks[kmer_idx];
            float lp_emission_m = log_probability_match(*data.read, rank, event_idx, data.strand);
            float lp_emission_e = log_probability_event_insert(*data.read, rank, event_idx, data.strand);

            // These functions either sum over the previous three states (forward algorithm)
            // or take the maximum and set the backtracking matrix (viterbi).
            output.update_m(row, curr_block_offset + PS_MATCH, m1, m2, m3, lp_emission_m);
            output.update_e(row, curr_block_offset + PS_EVENT_SPLIT, e1, e2, lp_emission_e);
            output.update_k(row, curr_block_offset + PS_KMER_SKIP, k1, k2, 0.0f); // no emission

#ifdef DEBUG_FILL    
            printf("Row %u block %u\n", row, block);
            printf("\tTransitions: p_mx [%.3lf %.3lf %.3lf]\n", bt.lp_mm, bt.lp_me, bt.lp_mk);
            printf("\t             p_ex [%.3lf %.3lf %.3lf]\n", bt.lp_em, bt.lp_ee, 0.0f);
            printf("\t             p_lx [%.3lf %.3lf %.3lf]\n", bt.lp_km, 0.0, bt.lp_kk);

            printf("\tPS_MATCH -- Transitions: [%.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mm, bt.lp_em, bt.lp_km, 
                    output.get(row - 1, prev_block_offset + PS_MATCH),
                    output.get(row - 1, prev_block_offset + PS_EVENT_SPLIT),
                    output.get(row - 1, prev_block_offset + PS_KMER_SKIP),
                    0.0f);
            printf("\tPS_EVENT_SPLIT -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_me, bt.lp_ee,
                    output.get(row - 1, curr_block_offset + PS_MATCH),
                    output.get(row - 1, curr_block_offset + PS_EVENT_SPLIT),
                    0.0f);

            printf("\tPS_KMER_SKIP -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mk, bt.lp_kk,
                    output.get(row, prev_block_offset + PS_MATCH),
                    output.get(row, prev_block_offset + PS_KMER_SKIP),
                    0.0f);

            printf("\tEMISSION: %.2lf %.2lf\n", lp_emission_m, lp_emission_e);
#endif
        }
    }
    
    uint32_t last_event_row = output.get_num_rows() - 1;
    uint32_t last_aligned_block = num_blocks - 2;
    uint32_t match_state_last_block = PS_NUM_STATES * last_aligned_block + PS_MATCH;
    return output.get(last_event_row, match_state_last_block);
}
