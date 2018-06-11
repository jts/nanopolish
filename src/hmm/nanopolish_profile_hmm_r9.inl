//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm_r9 -- Profile Hidden Markov Model
// for R9 data
//

//#define USE_EXTERNAL_PARAMS 1

#define TRANS_CLIP_SELF 0.9
#define TRANS_START_TO_CLIP 0.5

inline std::vector<BlockTransitions> calculate_transitions(uint32_t num_kmers, const HMMInputSequence& sequence, const HMMInputData& data)
{
    std::vector<BlockTransitions> transitions(num_kmers);
    
    double read_events_per_base = data.read->events_per_base[data.strand];

    for(uint32_t ki = 0; ki < num_kmers; ++ki) {

        // probability of skipping k_i from k_(i - 1)
        //float p_stay = 0.4;
        float p_stay = 1 - (1 / read_events_per_base); 
#ifndef USE_EXTERNAL_PARAMS
        float p_skip = 0.0025; 
        float p_bad = 0.001;
        float p_bad_self = p_bad;
        float p_skip_self = 0.3;
#else
        extern float g_p_skip, g_p_skip_self, g_p_bad, g_p_bad_self;
        float p_skip = g_p_skip;
        float p_skip_self = g_p_skip_self;
        float p_bad = g_p_bad;
        float p_bad_self = g_p_bad_self;
#endif
        // transitions from match state in previous block
        float p_mk = p_skip; // probability of not observing an event at all
        float p_mb = p_bad; // probabilty of observing a bad event
        float p_mm_self = p_stay; // probability of observing additional events from this k-mer
        float p_mm_next = 1.0f - p_mm_self - p_mk - p_mb; // normal movement from state to state

        // transitions from event split state in previous block
        float p_bb = p_bad_self;
        float p_bk, p_bm_next, p_bm_self;
        p_bk = p_bm_next = p_bm_self = (1.0f - p_bb) / 3;

        // transitions from kmer skip state in previous block
        float p_kk = p_skip_self;
        float p_km = 1.0f - p_kk;
        // p_kb not needed, equivalent to B->K

        // log-transform and store
        BlockTransitions& bt = transitions[ki];

        bt.lp_mk = log(p_mk);
        bt.lp_mb = log(p_mb);
        bt.lp_mm_self = log(p_mm_self);
        bt.lp_mm_next = log(p_mm_next);

        bt.lp_bb = log(p_bb);
        bt.lp_bk = log(p_bk);
        bt.lp_bm_next = log(p_bm_next);
        bt.lp_bm_self = log(p_bm_self);
        
        bt.lp_kk = log(p_kk);
        bt.lp_km = log(p_km);
    }

    return transitions;
}

// Output writer for the Forward Algorithm
class ProfileHMMForwardOutputR9
{
    public:
        ProfileHMMForwardOutputR9(FloatMatrix* p) : p_fm(p), lp_end(-INFINITY) {}
        
        //
        inline void update_cell(uint32_t row, uint32_t col, const HMMUpdateScores& scores, float lp_emission)
        {
            float sum = scores.x[0];
            for(auto i = 1; i < HMT_NUM_MOVEMENT_TYPES; ++i) {
                sum = add_logs(sum, scores.x[i]);
            }
            sum += lp_emission;
            set(*p_fm, row, col, sum);
        }

        // add in the probability of ending the alignment at row,col
        inline void update_end(float v, uint32_t, uint32_t)
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
        ProfileHMMForwardOutputR9(); // not allowed
        FloatMatrix* p_fm;
        float lp_end;
};

// Output writer for the Viterbi Algorithm
class ProfileHMMViterbiOutputR9
{
    public:
        ProfileHMMViterbiOutputR9(FloatMatrix* pf, UInt8Matrix* pb) : p_fm(pf), p_bm(pb), lp_end(-INFINITY) {}
        
        inline void update_cell(uint32_t row, uint32_t col, const HMMUpdateScores& scores, float lp_emission)
        {
            // probability update
            float max = scores.x[0];
            uint8_t from = 0;
            for(auto i = 1; i < HMT_NUM_MOVEMENT_TYPES; ++i) {
                max = scores.x[i] > max ? scores.x[i] : max;
                from = max == scores.x[i] ? i : from;
            }

            set(*p_fm, row, col, max + lp_emission);
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
        ProfileHMMViterbiOutputR9(); // not allowed

        FloatMatrix* p_fm;
        UInt8Matrix* p_bm;

        float lp_end;
        uint32_t end_row;
        uint32_t end_col;
};

// Allocate a vector with the model probabilities of skipping the first i events
inline std::vector<float> make_pre_flanking(const HMMInputData& data,
                                            const uint32_t e_start,
                                            const uint32_t num_events)
{
    std::vector<float> pre_flank(num_events + 1, 0.0f);
    
    // base cases

    // no skipping
    pre_flank[0] = log(1 - TRANS_START_TO_CLIP);

    // skipping the first event
    // this includes the transition probability into and out of the skip state
    pre_flank[1] = log(TRANS_START_TO_CLIP) + // transition from start to the background state
                   log_probability_background(*data.read, e_start, data.strand) + // emit from background
                   log(1 - TRANS_CLIP_SELF); // transition to silent pre state

    // skip the remaining events
    for(size_t i = 2; i < pre_flank.size(); ++i) {
        uint32_t event_idx = e_start + (i - 1) * data.event_stride;
        pre_flank[i] = log(TRANS_CLIP_SELF) + 
                       log_probability_background(*data.read, event_idx, data.strand) + // emit from background
                       pre_flank[i - 1]; // this accounts for the transition from the start & to the silent pre
    
    }

    return pre_flank;
}

// Allocate a vector with the model probabilities of skipping the remaining
// events after the alignment of event i
inline std::vector<float> make_post_flanking(const HMMInputData& data,
                                             const uint32_t e_start,
                                             const uint32_t num_events)
{
    // post_flank[i] means that the i-th event was the last one
    // aligned and the remainder should be emitted from the background model
    std::vector<float> post_flank(num_events, 0.0f);

    // base case, all events aligned
    post_flank[num_events - 1] = log(1 - TRANS_START_TO_CLIP);

    if(num_events > 1) {
        // base case, all events aligned but 1
        {
            uint32_t event_idx = e_start + (num_events - 1) * data.event_stride; // last event
            assert(event_idx == data.event_stop_idx);
            post_flank[num_events - 2] = log(TRANS_START_TO_CLIP) + // transition from pre to background state
                                         log_probability_background(*data.read, event_idx, data.strand) + // emit from background
                                         log(1 - TRANS_CLIP_SELF); // transition to silent pre state
        }

        for(int i = num_events - 3; i >= 0; --i) {
            uint32_t event_idx = e_start + (i + 1) * data.event_stride;
            post_flank[i] = log(TRANS_CLIP_SELF) +
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
inline float profile_hmm_fill_generic_r9(const HMMInputSequence& _sequence,
                                         const HMMInputData& _data,
                                         const uint32_t,
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
    
    // Calculate number of blocks
    // A block of the HMM is a set of states for one kmer
    uint32_t num_blocks = output.get_num_columns() / PSR9_NUM_STATES; // num_columns is the number of HMM STATES
    uint32_t last_event_row_idx = output.get_num_rows() - 1;

    // Precompute the transition probabilites for each kmer block
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks
    uint32_t last_kmer_idx = num_kmers - 1;
    
    std::vector<BlockTransitions> transitions = calculate_transitions(num_kmers, sequence, data);
 
    // Precompute kmer ranks
    const uint32_t k = data.pore_model->k;

    // Make sure the HMMInputSequence's alphabet matches the state space of the read
    assert( data.pore_model->states.size() == sequence.get_num_kmer_ranks(k) );

    std::vector<uint32_t> kmer_ranks(num_kmers);
    for(size_t ki = 0; ki < num_kmers; ++ki)
        kmer_ranks[ki] = sequence.get_kmer_rank(ki, k, data.rc);

    size_t num_events = output.get_num_rows() - 1;

    std::vector<float> pre_flank = make_pre_flanking(data, e_start, num_events);
    std::vector<float> post_flank = make_post_flanking(data, e_start, num_events);
    
    // The model is currently constrainted to always transition
    // from the terminal/clipped state to the first kmer (and from the
    // last kmer to the terminal/clipping state so these are log(1.0).
    // They are kept as variables as it might be relaxed later.
    float lp_sm, lp_ms;
    lp_sm = lp_ms = 0.0f;

    // the penalty is controlled by the transition probability
    float BAD_EVENT_PENALTY = 0.0f;

    // Fill in matrix
    for(uint32_t row = 1; row < output.get_num_rows(); row++) {

        // Skip the first block which is the start state, it was initialized above
        // Similarily skip the last block, which is calculated in the terminate() function
        for(uint32_t block = 1; block < num_blocks - 1; block++) {

            // retrieve transitions
            uint32_t kmer_idx = block - 1;
            BlockTransitions& bt = transitions[kmer_idx];

            uint32_t prev_block = block - 1;
            uint32_t prev_block_offset = PSR9_NUM_STATES * prev_block;
            uint32_t curr_block_offset = PSR9_NUM_STATES * block;
            
            // Emission probabilities
            uint32_t event_idx = e_start + (row - 1) * data.event_stride;
            uint32_t rank = kmer_ranks[kmer_idx];
            float lp_emission_m = log_probability_match_r9(*data.read, *data.pore_model, rank, event_idx, data.strand);
            float lp_emission_b = BAD_EVENT_PENALTY;
            
            HMMUpdateScores scores;

            // state PSR9_MATCH
            scores.x[HMT_FROM_SAME_M] = bt.lp_mm_self + output.get(row - 1, curr_block_offset + PSR9_MATCH);
            scores.x[HMT_FROM_PREV_M] = bt.lp_mm_next + output.get(row - 1, prev_block_offset + PSR9_MATCH);
            scores.x[HMT_FROM_SAME_B] = bt.lp_bm_self + output.get(row - 1, curr_block_offset + PSR9_BAD_EVENT);
            scores.x[HMT_FROM_PREV_B] = bt.lp_bm_next + output.get(row - 1, prev_block_offset + PSR9_BAD_EVENT);
            scores.x[HMT_FROM_PREV_K] = bt.lp_km + output.get(row - 1, prev_block_offset + PSR9_KMER_SKIP);

            // m_s is the probability of going from the start state
            // to this kmer. The start state is (currently) only 
            // allowed to go to the first kmer. If ALLOW_PRE_CLIP
            // is defined, we allow all events before this one to be skipped,
            // with a penalty;
            scores.x[HMT_FROM_SOFT] = (kmer_idx == 0 &&
                                        (event_idx == e_start ||
                                             (flags & HAF_ALLOW_PRE_CLIP))) ? lp_sm + pre_flank[row - 1] : -INFINITY;
            
            output.update_cell(row, curr_block_offset + PSR9_MATCH, scores, lp_emission_m);

            // state PSR9_BAD_EVENT
            scores.x[HMT_FROM_SAME_M] = bt.lp_mb + output.get(row - 1, curr_block_offset + PSR9_MATCH);
            scores.x[HMT_FROM_PREV_M] = -INFINITY; // not allowed
            scores.x[HMT_FROM_SAME_B] = bt.lp_bb + output.get(row - 1, curr_block_offset + PSR9_BAD_EVENT);
            scores.x[HMT_FROM_PREV_B] = -INFINITY;
            scores.x[HMT_FROM_PREV_K] = -INFINITY;
            scores.x[HMT_FROM_SOFT] = -INFINITY;
            output.update_cell(row, curr_block_offset + PSR9_BAD_EVENT, scores, lp_emission_b);

            // in cu this is where the shared memory sync on prev states would go.
            // state PSR9_KMER_SKIP
            scores.x[HMT_FROM_SAME_M] = -INFINITY;
            scores.x[HMT_FROM_PREV_M] = bt.lp_mk + output.get(row, prev_block_offset + PSR9_MATCH);
            scores.x[HMT_FROM_SAME_B] = -INFINITY;
            scores.x[HMT_FROM_PREV_B] = bt.lp_bk + output.get(row, prev_block_offset + PSR9_BAD_EVENT);
            scores.x[HMT_FROM_PREV_K] = bt.lp_kk + output.get(row, prev_block_offset + PSR9_KMER_SKIP);
            scores.x[HMT_FROM_SOFT] = -INFINITY;
            output.update_cell(row, curr_block_offset + PSR9_KMER_SKIP, scores, 0.0f); // no emission

            // If POST_CLIP is enabled we allow the last kmer to transition directly
            // to the end after any event. Otherwise we only allow it from the 
            // last kmer/event match.
            if(kmer_idx == last_kmer_idx && ( (flags & HAF_ALLOW_POST_CLIP) || row == last_event_row_idx)) {
                float lp1 = lp_ms + output.get(row, curr_block_offset + PSR9_MATCH) + post_flank[row - 1];
                float lp2 = lp_ms + output.get(row, curr_block_offset + PSR9_BAD_EVENT) + post_flank[row - 1];
                float lp3 = lp_ms + output.get(row, curr_block_offset + PSR9_KMER_SKIP) + post_flank[row - 1];

                output.update_end(lp1, row, curr_block_offset + PSR9_MATCH);
                output.update_end(lp2, row, curr_block_offset + PSR9_BAD_EVENT);
                output.update_end(lp3, row, curr_block_offset + PSR9_KMER_SKIP);
            }

#ifdef DEBUG_LOCAL_ALIGNMENT
            printf("[%d %d] start: %.2lf  pre: %.2lf fm: %.2lf\n", event_idx, kmer_idx, m_s + lp_emission_m, pre_flank[row - 1], output.get(row, curr_block_offset + PSR9_MATCH));
            printf("[%d %d]   end: %.2lf post: %.2lf\n", event_idx, kmer_idx, lp_end, post_flank[row - 1]);
#endif

#ifdef DEBUG_FILL    
            printf("Row %u block %u\n", row, block);

            printf("\tPSR9_MATCH -- Transitions: [%.3lf %.3lf %.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf %.2lf %.2lf] out: %.2lf\n", 
                    bt.lp_mm_self, bt.lp_mm_next, bt.lp_bm_self, bt.lp_bm_next, bt.lp_km, 
                    output.get(row - 1, prev_block_offset + PSR9_MATCH),
                    output.get(row - 1, curr_block_offset + PSR9_MATCH),
                    output.get(row - 1, prev_block_offset + PSR9_BAD_EVENT),
                    output.get(row - 1, curr_block_offset + PSR9_BAD_EVENT),
                    output.get(row - 1, prev_block_offset + PSR9_KMER_SKIP),
                    output.get(row, curr_block_offset + PSR9_MATCH));
            printf("\tPSR9_BAD_EVENT -- Transitions: [%.3lf %.3lf] Prev: [%.2lf %.2lf] out: %.2lf\n", 
                    bt.lp_mb, bt.lp_bb,
                    output.get(row - 1, curr_block_offset + PSR9_MATCH),
                    output.get(row - 1, curr_block_offset + PSR9_BAD_EVENT),
                    output.get(row, curr_block_offset + PSR9_BAD_EVENT));

            printf("\tPSR9_KMER_SKIP -- Transitions: [%.3lf %.3lf %.3lf] Prev: [%.2lf %.2lf %.2lf] sum: %.2lf\n", 
                    bt.lp_mk, bt.lp_bk, bt.lp_kk,
                    output.get(row, prev_block_offset + PSR9_MATCH),
                    output.get(row, prev_block_offset + PSR9_BAD_EVENT),
                    output.get(row, prev_block_offset + PSR9_KMER_SKIP),
                    output.get(row, curr_block_offset + PSR9_KMER_SKIP));

            printf("\tEMISSION: %.2lf %.2lf\n", lp_emission_m, lp_emission_b);
#endif
        }
    }
    
    return output.get_end();
}

