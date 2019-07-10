//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_raw_loader - utilities and helpers for loading
// data directly from raw nanopore files without events
//
#include "nanopolish_profile_hmm.h"
#include "nanopolish_haplotype.h"

//#define DEBUG_GENERIC 1

// Structure to keep track of the lower-left position in the band
struct BandOrigin
{
    int event_idx;
    int kmer_idx;
};

const uint8_t SHMM_FROM_D = 0;
const uint8_t SHMM_FROM_U = 1;
const uint8_t SHMM_FROM_L = 2;
const uint8_t SHMM_FROM_INVALID = 3;

class AdaptiveBandedViterbi
{
    public:
        AdaptiveBandedViterbi()
        {
            this->band_scores = NULL;
            this->trace = NULL;
            this->bandwidth = 0;
            this->n_fills = 0;
            this->n_bands = 0;
            this->initialized = false;
        }

        ~AdaptiveBandedViterbi()
        {
            free(this->band_scores);
            this->band_scores = NULL;
            
            free(this->trace);
            this->trace = NULL;
        }

        inline int get_start_state_offset(size_t band_idx) const
        {
            assert(band_idx == 0);
            int offset = get_offset_for_kmer_in_band(band_idx, -1);
            assert(get_offset_for_event_in_band(band_idx, -1) == offset);
            return offset;
        }

        inline int get_end_trim_kmer_state() const
        {
            return this->n_kmers;
        }
        
        inline int get_offset_for_event_in_band(size_t band_idx, int event_idx) const
        {
            return this->band_origins[band_idx].event_idx - event_idx;
        }

        inline int get_offset_for_kmer_in_band(size_t band_idx, int kmer_idx) const
        {
            return kmer_idx - this->band_origins[band_idx].kmer_idx;
        }

        inline int get_event_at_band_offset(size_t band_idx, int offset) const
        {
            return this->band_origins[band_idx].event_idx - offset;
        }

        inline int get_kmer_at_band_offset(size_t band_idx, int offset) const
        {
            return this->band_origins[band_idx].kmer_idx + offset;
        }

        inline int event_kmer_to_band(int event_idx, int kmer_idx) const
        {
            return (event_idx + 1) + (kmer_idx + 1);
        }

        inline bool is_offset_valid(int band_offset) const
        {
            return band_offset >= 0 && band_offset < this->bandwidth;
        }

        inline float get(size_t band_idx, int band_offset) const
        {
            return this->is_offset_valid(band_offset) ? 
                this->band_scores[band_idx * this->bandwidth + band_offset] : -INFINITY;
        }
        
        inline float get_by_event_kmer(int event_idx, int kmer_idx) const
        {
            size_t band_idx = event_kmer_to_band(event_idx, kmer_idx);
            int band_offset = get_offset_for_kmer_in_band(band_idx, kmer_idx);
            return this->get(band_idx, band_offset);
        }
        
        inline uint8_t get_trace(size_t band_idx, int band_offset) const
        {
            return this->is_offset_valid(band_offset) ? 
                this->trace[band_idx * this->bandwidth + band_offset] : 0;
        }
        

        inline void set(size_t band_idx, int band_offset, float value, uint8_t from) {
            size_t idx = band_idx * this->bandwidth + band_offset;
            this->band_scores[idx] = value;
            this->trace[idx] = from;
        }

        inline void set_start_state(float value)
        {
            size_t band_idx = 0;
            int offset = get_start_state_offset(band_idx);
            this->set(band_idx, offset, value, SHMM_FROM_D);
        }

        inline void set_start_trim_for_band(size_t band_idx, float log_trim_penalty)
        {
            int trim_offset = this->get_offset_for_kmer_in_band(band_idx, -1);
            if(!this->is_offset_valid(trim_offset)) {
                return;
            }

            int event_idx = this->get_event_at_band_offset(band_idx, trim_offset);
            if(event_idx >= 0) {
                this->set(band_idx, trim_offset, (event_idx + 1) * log_trim_penalty, SHMM_FROM_U);
            } else {
                this->set(band_idx, trim_offset, -INFINITY, SHMM_FROM_INVALID);
            }
        }

        inline void set3(size_t band_idx, int band_offset, float score_d, float score_u, float score_l)
        {
            float max_score = score_d;
            uint8_t from = SHMM_FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? SHMM_FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? SHMM_FROM_L : from;
#ifdef DEBUG_GENERIC
            int event_idx = this->get_event_at_band_offset(band_idx, band_offset);
            int kmer_idx = this->get_kmer_at_band_offset(band_idx, band_offset);
            fprintf(stderr, "[ada-generic] band: (%zu, %d) ek: (%d %d) set3(%.2f, %.2f, %.2f) from: %d\n", band_idx, band_offset, event_idx, kmer_idx, score_d, score_u, score_l, from);
#endif
            this->set(band_idx, band_offset, max_score, from);
            this->n_fills += 1;
        }

        inline void terminate()
        {
            int terminal_event_idx = this->n_events;
            int terminal_kmer_idx = this->n_kmers;

            float score_d = this->get_by_event_kmer(terminal_event_idx - 1, terminal_kmer_idx - 1);
            float score_u = this->get_by_event_kmer(terminal_event_idx - 1, terminal_kmer_idx);
            float score_l = this->get_by_event_kmer(terminal_event_idx, terminal_kmer_idx - 1);

            int band_idx = this->event_kmer_to_band(terminal_event_idx, terminal_kmer_idx);
            int band_offset = this->get_offset_for_kmer_in_band(band_idx, terminal_kmer_idx);
            this->set3(band_idx, band_offset, score_d, score_u, score_l);
        }

        inline BandOrigin move_band_down(const BandOrigin& curr_origin) const
        {
            return { curr_origin.event_idx + 1, curr_origin.kmer_idx };
        }
        
        inline BandOrigin move_band_right(const BandOrigin& curr_origin) const
        {
            return { curr_origin.event_idx, curr_origin.kmer_idx + 1 };
        }

        void initialize_common(const SquiggleRead& read, const std::string& sequence, size_t k, size_t strand_idx, const AdaBandedParameters& parameters)
        {
            this->initialized = true;
            this->parameters = parameters;
            this->n_events = read.events[strand_idx].size();
            this->n_kmers = sequence.size() - k + 1;
            this->bandwidth = parameters.bandwidth;

            // +2 for the start/end/trim states
            this->n_bands = (n_events + 2) + (n_kmers + 2);
            this->band_scores = (float*)malloc(sizeof(float) * this->n_bands * this->bandwidth);
            if(this->band_scores == NULL){
                fprintf(stderr,"Memory allocation failed at %s\n",__func__);
                exit(1);
            }

            this->trace = (uint8_t*)malloc(sizeof(uint8_t) * this->n_bands * this->bandwidth);
            if(this->trace == NULL){
                fprintf(stderr,"Memory allocation failed at %s\n",__func__);
                exit(1);
            }

            for (size_t i = 0; i < n_bands; i++) {
                for (int j = 0; j < bandwidth; j++) {
                    set(i, j, -INFINITY, 0);
                }
            }

            this->band_origins.resize(n_bands);
            for(size_t i = 0; i < this->band_origins.size(); ++i) {
                this->band_origins[i] = { -1, -1 };
            }
        }

        void initialize_adaptive(const SquiggleRead& read, 
                                 const std::string& sequence, 
                                 size_t k, 
                                 size_t strand_idx, 
                                 const AdaBandedParameters& parameters)
        {
            this->initialize_common(read, sequence, k, strand_idx, parameters);

            // initialize positions of first two bands
            int half_bandwidth = this->bandwidth / 2;
            this->band_origins[0].event_idx = half_bandwidth - 1;
            this->band_origins[0].kmer_idx = -1 - half_bandwidth;
            this->band_origins[1] = move_band_down(this->band_origins[0]);
        }

        void initialize_guided(const SquiggleRead& read, 
                               const Haplotype& haplotype, 
                               const EventAlignmentRecord& event_alignment_record, 
                               size_t k, 
                               size_t strand_idx, 
                               const AdaBandedParameters& parameters)
        {
            this->initialize_common(read, haplotype.get_sequence(), k, strand_idx, parameters);

            const std::vector<AlignedPair> event_to_haplotype_alignment = event_alignment_record.aligned_events;

            // make a map from event to the position in the haplotype it aligns to
            int ref_offset = haplotype.get_reference_position();
            std::vector<int> event_to_kmer(this->n_events, -1);
            for(size_t i = 0; i < event_to_haplotype_alignment.size(); ++i) {
                assert(event_to_haplotype_alignment[i].read_pos < event_to_kmer.size());
                event_to_kmer[event_to_haplotype_alignment[i].read_pos] = event_to_haplotype_alignment[i].ref_pos - ref_offset;
            }

            // First two bands need to be manually initialized
            int half_bandwidth = this->bandwidth / 2;
            this->band_origins[0].event_idx = half_bandwidth - 1;
            this->band_origins[0].kmer_idx = -1 - half_bandwidth;
            this->band_origins[1] = move_band_down(this->band_origins[0]);

            // many events will have an unassigned position, fill them in here using the last seen kmer
            // and set band coordinates
            int prev_kmer_idx = 0;
            for(size_t i = 0; i < event_to_kmer.size(); ++i) {
                if(event_to_kmer[i] == -1) {
                    event_to_kmer[i] = prev_kmer_idx;
                } else {
                    prev_kmer_idx = event_to_kmer[i];
                }

                int event_idx = i;
                int kmer_idx = event_to_kmer[i];
                int band_idx = this->event_kmer_to_band(event_idx, kmer_idx);
                if(band_idx < this->band_origins.size()) {
                    this->band_origins[band_idx].event_idx = half_bandwidth + event_idx;
                    this->band_origins[band_idx].kmer_idx = kmer_idx - half_bandwidth;
                }
                /*
                fprintf(stderr, "lookup[%zu] = %d band_idx: %d origin: %d %d\n", 
                    event_idx, kmer_idx, band_idx, this->band_origins[band_idx].event_idx, this->band_origins[band_idx].kmer_idx);
                */
            }
        }

        int get_num_bands() const { return this->n_bands; }
        int get_num_fills() const { return this->n_fills; }
        bool is_initialized() const { return this->initialized; }

        void determine_band_origin(size_t band_idx)
        {

            // band position already set, do nothing
            if(this->band_origins[band_idx].event_idx >= 0 && this->band_origins[band_idx].kmer_idx >= 0) {
                return;
            }
            // Determine placement of this band according to Suzuki's adaptive algorithm
            // When both ll and ur are out-of-band (ob) we alternate movements
            // otherwise we decide based on scores
            float ll = this->get(band_idx - 1, 0);
            float ur = this->get(band_idx - 1, this->bandwidth - 1);
            bool ll_ob = ll == -INFINITY;
            bool ur_ob = ur == -INFINITY;

            bool right = false;
            if(ll_ob && ur_ob) {
                right = band_idx % 2 == 1;
            } else {
                right = ll < ur; // Suzuki's rule
            }

            if(right) {
                this->band_origins[band_idx] = move_band_right(this->band_origins[band_idx - 1]);
            } else {
                this->band_origins[band_idx] = move_band_down(this->band_origins[band_idx - 1]);
            }
        }

        void get_offset_range_for_band(size_t band_idx, int& min_offset, int& max_offset) const
        {
            // Get the offsets for the first and last event and kmer
            // We restrict the inner loop to only these values
            int kmer_min_offset = this->get_offset_for_kmer_in_band(band_idx, 0);
            int kmer_max_offset = this->get_offset_for_kmer_in_band(band_idx, this->n_kmers);

            int event_min_offset = this->get_offset_for_event_in_band(band_idx, this->n_events - 1);
            int event_max_offset = this->get_offset_for_event_in_band(band_idx, -1);

            min_offset = std::max(kmer_min_offset, event_min_offset);
            min_offset = std::max(min_offset, 0);

            max_offset = std::min(kmer_max_offset, event_max_offset);
            max_offset = std::min(max_offset, (int)this->bandwidth);
        }

        std::vector<AlignedPair> backtrack() const
        {
            // Backtrack to compute alignment
            std::vector<AlignedPair> out;

            float max_score = -INFINITY;
            int curr_event_idx = this->n_events;
            int curr_kmer_idx = this->n_kmers;

#ifdef DEBUG_GENERIC
            fprintf(stderr, "[ada-generic-back] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, this->get_by_event_kmer(curr_event_idx, curr_kmer_idx));
#endif

            while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {

                // emit current alignment
                if(curr_kmer_idx != this->n_kmers) {
                    out.push_back({curr_kmer_idx, curr_event_idx});
                }

#ifdef DEBUG_GENERIC
                fprintf(stderr, "[ada-generic-back] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx);
#endif
                // position in band
                size_t band_idx = this->event_kmer_to_band(curr_event_idx, curr_kmer_idx);
                size_t offset = this->get_offset_for_event_in_band(band_idx, curr_event_idx);
                assert(this->get_offset_for_kmer_in_band(band_idx, curr_kmer_idx) == offset);

                uint8_t from = this->get_trace(band_idx,offset);
                if(from == SHMM_FROM_D) {
                    curr_kmer_idx -= 1;
                    curr_event_idx -= 1;
                } else if(from == SHMM_FROM_U) {
                    curr_event_idx -= 1;
                } else {
                    curr_kmer_idx -= 1;
                }
            }
            std::reverse(out.begin(), out.end());
            return out;
        }

    private:

        float* band_scores;
        uint8_t* trace;
        std::vector<BandOrigin> band_origins;
        AdaBandedParameters parameters;
        size_t n_kmers;
        size_t n_events;
        size_t n_bands;
        size_t n_fills;
        size_t bandwidth;
        bool initialized;
};

template<class GenericBandedHMMResult>
void generic_banded_simple_hmm(SquiggleRead& read,
                               const PoreModel& pore_model,
                               const std::string& sequence,
                               const AdaBandedParameters parameters,
                               GenericBandedHMMResult& hmm_result)
{
    size_t strand_idx = 0;
    size_t k = pore_model.k;
    const Alphabet* alphabet = pore_model.pmalphabet;
    size_t n_events = read.events[strand_idx].size();
    size_t n_kmers = sequence.size() - k + 1;

#ifdef DEBUG_GENERIC
    fprintf(stderr, "[ada] aligning read %s\n", read.read_name.substr(0,6).c_str());
#endif

    // transition penalties
    double events_per_kmer = (double)n_events / n_kmers;
    double p_stay = 1 - (1 / events_per_kmer);
    double lp_skip = log(parameters.p_skip);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(parameters.p_trim);
 
    // Initialize

    // Precompute k-mer ranks
    std::vector<size_t> kmer_ranks(n_kmers);
    for(size_t i = 0; i < n_kmers; ++i) {
        kmer_ranks[i] = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
    }

    assert(hmm_result.is_initialized());

    // initialize first two bands as a special case

    // band 0: score zero in the central cell
    hmm_result.set_start_state(0.0f);

    // band 1: set trim
    hmm_result.set_start_trim_for_band(1, lp_trim);

#ifdef DEBUG_GENERIC
    fprintf(stderr, "[generic] trim-init bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, hmm_result.get(1,first_trim_offset));
#endif

    // fill in remaining bands
    for(int band_idx = 2; band_idx < hmm_result.get_num_bands(); ++band_idx) {

        hmm_result.determine_band_origin(band_idx);

        // update trim state
        hmm_result.set_start_trim_for_band(band_idx, lp_trim);

        // determine the range of offsets in this band we should fill in
        int min_offset, max_offset;
        hmm_result.get_offset_range_for_band(band_idx, min_offset, max_offset);

        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = hmm_result.get_event_at_band_offset(band_idx, offset);
            int kmer_idx = hmm_result.get_kmer_at_band_offset(band_idx, offset);

            size_t kmer_rank = kmer_ranks[kmer_idx];
            
            float up   = hmm_result.get_by_event_kmer(event_idx - 1, kmer_idx);
            float left = hmm_result.get_by_event_kmer(event_idx, kmer_idx - 1);
            float diag = hmm_result.get_by_event_kmer(event_idx - 1, kmer_idx - 1);

            float lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
            float score_d = diag + lp_step + lp_emission;
            float score_u = up + lp_stay + lp_emission;
            float score_l = left + lp_skip;
            hmm_result.set3(band_idx, offset, score_d, score_u, score_l);

#ifdef DEBUG_GENERIC
            fprintf(stderr, "[ada-gen-fill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left);
            fprintf(stderr, "[ada-gen-fill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
            fprintf(stderr, "[ada-gen-fill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d rank: %zu emit: %.2lf\n", 
                band_idx, offset, event_idx, kmer_idx, hmm_result.get(band_idx, offset), hmm_result.get_trace(band_idx, offset), kmer_rank, lp_emission);
#endif
        }

        // if there is an end trim state in this band, set it here
        int end_trim_kmer_state = hmm_result.get_end_trim_kmer_state();
        assert(end_trim_kmer_state == n_kmers);
        int offset = hmm_result.get_offset_for_kmer_in_band(band_idx, end_trim_kmer_state);
        if(hmm_result.is_offset_valid(offset)) {
            int event_idx = hmm_result.get_event_at_band_offset(band_idx, offset);
            float score_d = hmm_result.get_by_event_kmer(event_idx - 1, n_kmers - 1) + lp_trim;
            float score_u = hmm_result.get_by_event_kmer(event_idx - 1, end_trim_kmer_state) + lp_trim;
            float score_l = hmm_result.get_by_event_kmer(event_idx, n_kmers - 1) + lp_trim;
            hmm_result.set3(band_idx, offset, score_d, score_u, score_l);
        }
    }

    hmm_result.terminate();

/*
    // Debug, print some of the score matrix
    for(int col = 0; col <= 10; ++col) {
        for(int row = 0; row < 100; ++row) {
            int kmer_idx = col - 1;
            int event_idx = row - 1;
            int band_idx = hmm_result.event_kmer_to_band(event_idx, kmer_idx);
            int offset = hmm_result.get_offset_for_kmer_in_band(band_idx, kmer_idx);
            assert(offset == hmm_result.get_offset_for_event_in_band(band_idx, event_idx));
            assert(event_idx == hmm_result.get_event_at_band_offset(band_idx, offset));
            fprintf(stderr, "[ada-gen-fill] ei: %d ki: %d bi: %d o: %d s: %.2f\n", event_idx, kmer_idx, band_idx, offset, hmm_result.get(band_idx, offset));
        }
    }
*/
}
