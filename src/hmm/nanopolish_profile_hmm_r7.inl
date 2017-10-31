//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
// for R7 data
//
inline float calculate_skip_probability_r7(const HMMInputSequence& sequence,
                                           const HMMInputData& data,
                                           uint32_t ki,
                                           uint32_t kj)
{
    const PoreModel& pm = *data.pore_model;
    const TransitionParameters& parameters = data.read->parameters[data.strand];

    uint32_t rank_i = sequence.get_kmer_rank(ki, pm.k, data.rc);
    uint32_t rank_j = sequence.get_kmer_rank(kj, pm.k, data.rc);

    GaussianParameters level_i = data.read->get_scaled_gaussian_from_pore_model_state(*data.pore_model, data.strand, rank_i);
    GaussianParameters level_j = data.read->get_scaled_gaussian_from_pore_model_state(*data.pore_model, data.strand, rank_j);

    return parameters.get_skip_probability(level_i.mean, level_j.mean);
}

inline std::vector<BlockTransitionsR7> calculate_transitions_r7(uint32_t num_kmers, const HMMInputSequence& sequence, const HMMInputData& data)
{
    const TransitionParameters& parameters = data.read->parameters[data.strand];

    std::vector<BlockTransitionsR7> transitions(num_kmers);
    
    for(uint32_t ki = 0; ki < num_kmers; ++ki) {

        // probability of skipping k_i from k_(i - 1)
        float p_skip = ki > 0 ? calculate_skip_probability_r7(sequence, data, ki - 1, ki) : 0.0f;

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
        BlockTransitionsR7& bt = transitions[ki];

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
class ProfileHMMForwardOutputR7
{
    public:
        ProfileHMMForwardOutputR7(FloatMatrix* p) : p_fm(p), lp_end(-INFINITY) {}
        
        //
        inline void update_4(uint32_t row, uint32_t col, float m, float e, float k, float s, float lp_emission)
        {
            float sum_1 = add_logs(m, e);
            float sum_2 = add_logs(k, s);
            float sum = add_logs(sum_1, sum_2) + lp_emission;
            set(*p_fm, row, col, sum);
        }

        // add in the probability of ending the alignment at row,col
        inline void update_end(float v, uint32_t row, uint32_t col)
        {
            lp_end = add_logs(lp_end, v);
        }

        // get the log probability stored at a particular row/column
        inline float get(uint32_t row, uint32_t col) const
        {
            return ::get(*p_fm, row, col);
        }

        // get the log probability for the end state
        inline float get_end() const
        {
            return lp_end;
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
        ProfileHMMForwardOutputR7(); // not allowed
        FloatMatrix* p_fm;
        float lp_end;
};

// Output writer for the Viterbi Algorithm
class ProfileHMMViterbiOutputR7
{
    public:
        ProfileHMMViterbiOutputR7(FloatMatrix* pf, UInt8Matrix* pb) : p_fm(pf), p_bm(pb), lp_end(-INFINITY) {}
        
        inline void update_4(uint32_t row, uint32_t col, float m, float e, float k, float s, float lp_emission)
        {
            // probability update
            float max = std::max(std::max(m, e), 
                                 std::max(k, s));

            set(*p_fm, row, col, max + lp_emission);

            // backtrack update
            uint8_t from;
            if(max == m)
                from = PSR7_MATCH;
            else if(max == e)
                from = PSR7_EVENT_SPLIT;
            else if(max == k)
                from = PSR7_KMER_SKIP;
            else if(max == s)
                from = PSR7_PRE_SOFT;
            set(*p_bm, row, col, from);
        }
        
        // add in the probability of ending the alignment at row,col
        inline void update_end(float v, uint32_t row, uint32_t col)
        {
            if(v > lp_end) {
                lp_end = v;
                end_row = row;
                end_col = col;
            }
        }

        // get the log probability stored at a particular row/column
        inline float get(uint32_t row, uint32_t col) const
        {
            return ::get(*p_fm, row, col);
        }

        // get the log probability for the end state
        inline float get_end() const
        {
            return lp_end;
        }
        
        // get the row/col that lead to the end state
        inline void get_end_cell(uint32_t& row, uint32_t& col)
        {
            row = end_row;
            col = end_col;
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
        ProfileHMMViterbiOutputR7(); // not allowed

        FloatMatrix* p_fm;
        UInt8Matrix* p_bm;

        float lp_end;
        uint32_t end_row;
        uint32_t end_col;
};

// Allocate a vector with the model probabilities of skipping the first i events
inline std::vector<float> make_pre_flanking_r7(const HMMInputData& data,
                                               const TransitionParameters& parameters,
                                               const uint32_t e_start,
                                               const uint32_t num_events)
{
    std::vector<float> pre_flank(num_events + 1, 0.0f);
    
    // base cases

    // no skipping
    pre_flank[0] = log(1 - parameters.trans_start_to_clip);

    // skipping the first event
    // this includes the transition probability into and out of the skip state
    pre_flank[1] = log(parameters.trans_start_to_clip) + // transition from start to the background state
                   log_probability_background(*data.read, e_start, data.strand) + // emit from background
                   log(1 - parameters.trans_clip_self); // transition to silent pre state

    // skip the remaining events
    for(size_t i = 2; i < pre_flank.size(); ++i) {
        uint32_t event_idx = e_start + (i - 1) * data.event_stride;
        pre_flank[i] = log(parameters.trans_clip_self) +
                       log_probability_background(*data.read, event_idx, data.strand) + // emit from background
                       pre_flank[i - 1]; // this accounts for the transition from the start & to the silent pre
    
    }

    return pre_flank;
}

// Allocate a vector with the model probabilities of skipping the remaining
// events after the alignment of event i
inline std::vector<float> make_post_flanking_r7(const HMMInputData& data,
                                                const TransitionParameters& parameters,
                                                const uint32_t e_start,
                                                const uint32_t num_events)
{
    // post_flank[i] means that the i-th event was the last one
    // aligned and the remainder should be emitted from the background model
    std::vector<float> post_flank(num_events, 0.0f);

    // base case, all events aligned
    post_flank[num_events - 1] = log(1 - parameters.trans_start_to_clip);

    if(num_events > 1) {
        // base case, all events aligned but 1
        {
            uint32_t event_idx = e_start + (num_events - 1) * data.event_stride; // last event
            assert(event_idx == data.event_stop_idx);
            post_flank[num_events - 2] = log(parameters.trans_start_to_clip) + // transition from pre to background state
                                         log_probability_background(*data.read, event_idx, data.strand) + // emit from background
                                         log(1 - parameters.trans_clip_self); // transition to silent pre state
        }

        for(int i = num_events - 3; i >= 0; --i) {
            uint32_t event_idx = e_start + (i + 1) * data.event_stride;
            post_flank[i] = log(parameters.trans_clip_self) +
                            log_probability_background(*data.read, event_idx, data.strand) + // emit from background
                            post_flank[i + 1]; // this accounts for the transition from start, and to silent pre
        }
    }
    return post_flank;
}

// This function fills in a matrix with the result of running the HMM.
// The templated ProfileHMMOutput class allows one to run either Viterbi
// or the Forward algorithm.
template<class ProfileHMMOutput>
inline float profile_hmm_fill_generic_r7(const HMMInputSequence& _sequence,
                                         const HMMInputData& _data,
                                         const uint32_t _e_start,
                                         uint32_t flags,
                                         ProfileHMMOutput& output)
{
    PROFILE_FUNC("profile_hmm_fill_generic")
    HMMInputSequence sequence = _sequence;
    HMMInputData data = _data;
    assert( (data.rc && data.event_stride == -1) || (!data.rc && data.event_stride == 1));

#if HMM_REVERSE_FIX
    if(data.event_stride == -1) {
        sequence.swap();
        uint32_t tmp = data.event_stop_idx;
        data.event_stop_idx = data.event_start_idx;
        data.event_start_idx = tmp;
        data.event_stride = 1;
        data.rc = false;
    }
#endif

    uint32_t e_start = data.event_start_idx;
    
    const TransitionParameters& parameters = data.read->parameters[data.strand];

    // Calculate number of blocks
    // A block of the HMM is a set of PSR7_KMER_SKIP, PSR7_EVENT_SPLIT, PSR7_MATCH
    // events for one kmer
    uint32_t num_blocks = output.get_num_columns() / PSR7_NUM_STATES;
    uint32_t last_event_row_idx = output.get_num_rows() - 1;

    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks
    uint32_t last_kmer_idx = num_kmers - 1;
    
    std::vector<BlockTransitionsR7> transitions = calculate_transitions_r7(num_kmers, sequence, data);
 
    // Precompute kmer ranks
    const uint32_t k = data.pore_model->k;

    // Make sure the HMMInputSequence's alphabet matches the state space of the read
    assert( data.pore_model->states.size() == sequence.get_num_kmer_ranks(k) );

    std::vector<uint32_t> kmer_ranks(num_kmers);
    for(size_t ki = 0; ki < num_kmers; ++ki)
        kmer_ranks[ki] = sequence.get_kmer_rank(ki, k, data.rc);

    size_t num_events = output.get_num_rows() - 1;

    std::vector<float> pre_flank = make_pre_flanking_r7(data, parameters, e_start, num_events);
    std::vector<float> post_flank = make_post_flanking_r7(data, parameters, e_start, num_events);
    
    // The model is currently constrainted to always transition
    // from the terminal/clipped state to the first kmer (and from the
    // last kmer to the terminal/clipping state so these are log(1.0).
    // They are kept as variables as it might be relaxed later.
    float lp_sm, lp_ms;
    lp_sm = lp_ms = 0.0f;

    // Fill in matrix
    for(uint32_t row = 1; row < output.get_num_rows(); row++) {

        // Skip the first block which is the start state, it was initialized above
        // Similarily skip the last block, which is calculated in the terminate() function
        for(uint32_t block = 1; block < num_blocks - 1; block++) {

            // retrieve transitions
            uint32_t kmer_idx = block - 1;
            BlockTransitionsR7& bt = transitions[kmer_idx];

            uint32_t prev_block = block - 1;
            uint32_t prev_block_offset = PSR7_NUM_STATES * prev_block;
            uint32_t curr_block_offset = PSR7_NUM_STATES * block;
            
            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * data.event_stride;
            uint32_t rank = kmer_ranks[kmer_idx];
            float lp_emission_m = log_probability_match_r7(*data.read, *data.pore_model, rank, event_idx, data.strand);
            float lp_emission_e = log_probability_event_insert_r7(*data.read, *data.pore_model, rank, event_idx, data.strand);
            
            // state PSR7_MATCH
            float m_m = bt.lp_mm + output.get(row - 1, prev_block_offset + PSR7_MATCH);
            float m_e = bt.lp_em + output.get(row - 1, prev_block_offset + PSR7_EVENT_SPLIT);
            float m_k = bt.lp_km + output.get(row - 1, prev_block_offset + PSR7_KMER_SKIP);

            // m_s is the probability of going from the start state
            // to this kmer. The start state is (currently) only 
            // allowed to go to the first kmer. If ALLOW_PRE_CLIP
            // is defined, we allow all events before this one to be skipped,
            // with a penalty;
            float m_s = (kmer_idx == 0 &&
                            (event_idx == e_start ||
                             (flags & HAF_ALLOW_PRE_CLIP))) ? lp_sm + pre_flank[row - 1] : -INFINITY;
            
            output.update_4(row, curr_block_offset + PSR7_MATCH, m_m, m_e, m_k, m_s, lp_emission_m);

            // state PSR7_EVENT_SPLIT
            float e_m = bt.lp_me + output.get(row - 1, curr_block_offset + PSR7_MATCH);
            float e_e = bt.lp_ee + output.get(row - 1, curr_block_offset + PSR7_EVENT_SPLIT);
            output.update_4(row, curr_block_offset + PSR7_EVENT_SPLIT, e_m, e_e, -INFINITY, -INFINITY, lp_emission_e);

            // state PSR7_KMER_SKIP
            float k_m = bt.lp_mk + output.get(row, prev_block_offset + PSR7_MATCH);
            float k_k = bt.lp_kk + output.get(row, prev_block_offset + PSR7_KMER_SKIP);
            output.update_4(row, curr_block_offset + PSR7_KMER_SKIP, k_m, -INFINITY, k_k, -INFINITY, 0.0f); // no emission

            // If POST_CLIP is enabled we allow the last kmer to transition directly
            // to the end after any event. Otherwise we only allow it from the 
            // last kmer/event match.
            if(kmer_idx == last_kmer_idx && ( (flags & HAF_ALLOW_POST_CLIP) || row == last_event_row_idx)) {
                float lp1 = lp_ms + output.get(row, curr_block_offset + PSR7_MATCH) + post_flank[row - 1];
                float lp2 = lp_ms + output.get(row, curr_block_offset + PSR7_EVENT_SPLIT) + post_flank[row - 1];
                float lp3 = lp_ms + output.get(row, curr_block_offset + PSR7_KMER_SKIP) + post_flank[row - 1];

                output.update_end(lp1, row, curr_block_offset + PSR7_MATCH);
                output.update_end(lp2, row, curr_block_offset + PSR7_EVENT_SPLIT);
                output.update_end(lp3, row, curr_block_offset + PSR7_KMER_SKIP);
            }

#ifdef DEBUG_LOCAL_ALIGNMENT
            printf("[%d %d] start: %.2lf  pre: %.2lf fm: %.2lf\n", event_idx, kmer_idx, m_s + lp_emission_m, pre_flank[row - 1], output.get(row, curr_block_offset + PSR7_MATCH));
            printf("[%d %d]   end: %.2lf post: %.2lf\n", event_idx, kmer_idx, lp_end, post_flank[row - 1]);
#endif

#ifdef DEBUG_FILL    
            printf("Row %u block %u\n", row, block);
            printf("\tTransitions: p_mx [%.3lf %.3lf %.3lf]\n", bt.lp_mm, bt.lp_me, bt.lp_mk);
            printf("\t             p_ex [%.3lf %.3lf %.3lf]\n", bt.lp_em, bt.lp_ee, 0.0f);
            printf("\t             p_lx [%.3lf %.3lf %.3lf]\n", bt.lp_km, 0.0, bt.lp_kk);

            printf("\tPSR7_MATCH -- Transitions: [%.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mm, bt.lp_em, bt.lp_km, 
                    output.get(row - 1, prev_block_offset + PSR7_MATCH),
                    output.get(row - 1, prev_block_offset + PSR7_EVENT_SPLIT),
                    output.get(row - 1, prev_block_offset + PSR7_KMER_SKIP),
                    0.0f);
            printf("\tPSR7_EVENT_SPLIT -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_me, bt.lp_ee,
                    output.get(row - 1, curr_block_offset + PSR7_MATCH),
                    output.get(row - 1, curr_block_offset + PSR7_EVENT_SPLIT),
                    0.0f);

            printf("\tPSR7_KMER_SKIP -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mk, bt.lp_kk,
                    output.get(row, prev_block_offset + PSR7_MATCH),
                    output.get(row, prev_block_offset + PSR7_KMER_SKIP),
                    0.0f);

            printf("\tEMISSION: %.2lf %.2lf\n", lp_emission_m, lp_emission_e);
#endif
        }
    }
    
    return output.get_end();
}

