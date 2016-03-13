//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_transition_parameters -- parameters for khmm model
//
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "nanopolish_transition_parameters.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_model_names.h"

TransitionParameters::TransitionParameters()
{
    // initialize training data
    TransitionTrainingData& td = training_data;
    td.n_matches = 0;
    td.n_merges = 0;
    td.n_skips = 0;

    //
    allocate_matrix(td.state_transitions, 3, 3);
    for(unsigned i = 0; i < td.state_transitions.n_rows; ++i) {
        for(unsigned j = 0; j < td.state_transitions.n_cols; ++j) {
            set(td.state_transitions, i, j, 0);
        }
    }

    //
    // Initialize transition parameters
    //

    // these are fixed across all models
    skip_bin_width = 0.5;
    skip_probabilities.resize(30);
 
    trans_start_to_clip = 0.5f;
    trans_clip_self = 0.90f;
}

//
TransitionParameters::~TransitionParameters()
{
    free_matrix(training_data.state_transitions);
}

void TransitionParameters::initialize(const std::string& model_name)
{
    is_initialized = true;
    ModelMetadata model_data = get_model_metadata_from_name(model_name);

    if(model_data.kit == KV_SQK005) {
        initialize_sqkmap005();
    } else if(model_data.kit == KV_SQK006) {
        if(model_data.strand_idx == T_IDX) {
            initialize_sqkmap006_template();
        } else {
            initialize_sqkmap006_complement();
        }
    } else {
        printf("Warning: unknown model: %s\n", model_name.c_str());
        initialize_sqkmap005();
    }
}

void TransitionParameters::initialize_sqkmap005()
{
    assert(!skip_probabilities.empty());
    trans_m_to_e_not_k = 0.15f;
    trans_e_to_e = 0.33f;

    skip_probabilities[0] = 0.51268137;
    skip_probabilities[1] = 0.47243219;
    skip_probabilities[2] = 0.42888741;
    skip_probabilities[3] = 0.34932588;
    skip_probabilities[4] = 0.27427068;
    skip_probabilities[5] = 0.22297225;
    skip_probabilities[6] = 0.17585147;
    skip_probabilities[7] = 0.14705882;
    skip_probabilities[8] = 0.12183525;
    skip_probabilities[9] = 0.11344997;
    skip_probabilities[10] = 0.10069393;
    skip_probabilities[11] = 0.09153005;
    skip_probabilities[12] = 0.08765206;
    skip_probabilities[13] = 0.08491435;
    skip_probabilities[14] = 0.08272553;
    skip_probabilities[15] = 0.07747396;
    skip_probabilities[16] = 0.08439116;
    skip_probabilities[17] = 0.07819045;
    skip_probabilities[18] = 0.07337461;
    skip_probabilities[19] = 0.07020490;
    skip_probabilities[20] = 0.06869961;
    skip_probabilities[21] = 0.06576609;
    skip_probabilities[22] = 0.06923376;
    skip_probabilities[23] = 0.06239092;
    skip_probabilities[24] = 0.06586513;
    skip_probabilities[25] = 0.07372986;
    skip_probabilities[26] = 0.07050360;
    skip_probabilities[27] = 0.07228916;
    skip_probabilities[28] = 0.05855856;
    skip_probabilities[29] = 0.06842737;

}

void TransitionParameters::initialize_sqkmap006_template()
{
    assert(!skip_probabilities.empty());
    trans_m_to_e_not_k = 0.17f;
    trans_e_to_e = 0.55f;

    skip_probabilities[0] = 0.487;
    skip_probabilities[1] = 0.412;
    skip_probabilities[2] = 0.311;
    skip_probabilities[3] = 0.229;
    skip_probabilities[4] = 0.174;
    skip_probabilities[5] = 0.134;
    skip_probabilities[6] = 0.115;
    skip_probabilities[7] = 0.103;
    skip_probabilities[8] = 0.096;
    skip_probabilities[9] = 0.092;
    skip_probabilities[10] = 0.088;
    skip_probabilities[11] = 0.087;
    skip_probabilities[12] = 0.084;
    skip_probabilities[13] = 0.085;
    skip_probabilities[14] = 0.083;
    skip_probabilities[15] = 0.082;
    skip_probabilities[16] = 0.085;
    skip_probabilities[17] = 0.083;
    skip_probabilities[18] = 0.084;
    skip_probabilities[19] = 0.082;
    skip_probabilities[20] = 0.080;
    skip_probabilities[21] = 0.085;
    skip_probabilities[22] = 0.088;
    skip_probabilities[23] = 0.086;
    skip_probabilities[24] = 0.087;
    skip_probabilities[25] = 0.089;
    skip_probabilities[26] = 0.085;
    skip_probabilities[27] = 0.090;
    skip_probabilities[28] = 0.087;
    skip_probabilities[29] = 0.096;
}

void TransitionParameters::initialize_sqkmap006_complement()
{
    assert(!skip_probabilities.empty());
    trans_m_to_e_not_k = 0.14f;
    trans_e_to_e = 0.49f;

    skip_probabilities[0] = 0.531;
    skip_probabilities[1] = 0.478;
    skip_probabilities[2] = 0.405;
    skip_probabilities[3] = 0.327;
    skip_probabilities[4] = 0.257;
    skip_probabilities[5] = 0.207;
    skip_probabilities[6] = 0.172;
    skip_probabilities[7] = 0.154;
    skip_probabilities[8] = 0.138;
    skip_probabilities[9] = 0.132;
    skip_probabilities[10] = 0.127;
    skip_probabilities[11] = 0.123;
    skip_probabilities[12] = 0.117;
    skip_probabilities[13] = 0.115;
    skip_probabilities[14] = 0.113;
    skip_probabilities[15] = 0.113;
    skip_probabilities[16] = 0.115;
    skip_probabilities[17] = 0.109;
    skip_probabilities[18] = 0.109;
    skip_probabilities[19] = 0.107;
    skip_probabilities[20] = 0.104;
    skip_probabilities[21] = 0.105;
    skip_probabilities[22] = 0.108;
    skip_probabilities[23] = 0.106;
    skip_probabilities[24] = 0.111;
    skip_probabilities[25] = 0.114;
    skip_probabilities[26] = 0.118;
    skip_probabilities[27] = 0.119;
    skip_probabilities[28] = 0.110;
    skip_probabilities[29] = 0.119;
}

// 
double TransitionParameters::get_skip_probability(double k_level1, double k_level2) const
{
    assert(is_initialized);
    size_t bin = get_skip_bin(k_level1, k_level2);
    assert(bin < skip_probabilities.size());
    return skip_probabilities[bin];
}

//
int statechar2index(char s)
{
    switch(s) {
        case 'M': return 0;
        case 'E': return 1;
        case 'K': return 2;
    }
    assert(false);
    return 0;
}

//
void TransitionParameters::add_transition_observation(char state_from, char state_to)
{
    int f_idx = statechar2index(state_from);
    int t_idx = statechar2index(state_to);

    int count = get(training_data.state_transitions, f_idx, t_idx);
    set(training_data.state_transitions, f_idx, t_idx, count + 1);
}

void TransitionParameters::add_training_from_alignment(const HMMInputSequence& sequence,
                                                       const HMMInputData& data,
                                                       const std::vector<HMMAlignmentState>& alignment,
                                                       size_t ignore_edge_length)
{
    // do nothing if the alignment is too short
    if(alignment.size() <= ignore_edge_length) {
        return;
    }

    const PoreModel& pm = data.read->pore_model[data.strand];
    const uint32_t k = pm.k;

    size_t n_kmers = sequence.length() - k + 1;
#ifdef PRINT_TRAINING_MESSAGES
    uint32_t strand_idx = 0;
#endif
    char prev_s = 'M';

    for(size_t pi = 0; pi < alignment.size(); ++pi) {

        uint32_t ei = alignment[pi].event_idx;
        uint32_t ki = alignment[pi].kmer_idx;
        char s = alignment[pi].state;
    
        // Record transition observations
        // We do not record observations for merge states as there was no kmer transitions
        // We also do not record observations for the beginning of the matches as the
        // alignment may be poor due to edge effects
        if(pi > ignore_edge_length && pi < alignment.size() - ignore_edge_length) {
 
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

                uint32_t rank_1 = sequence.get_kmer_rank(transition_kmer_from, k, data.rc);
                uint32_t rank_2 = sequence.get_kmer_rank(transition_kmer_to, k, data.rc);
            
                GaussianParameters level_1 = pm.get_scaled_parameters(rank_1);
                GaussianParameters level_2 = pm.get_scaled_parameters(rank_2);
            
#ifdef PRINT_TRAINING_MESSAGES
                printf("TRAIN_SKIP\t%d\t%.3lf\t%.3lf\t%c\n", strand_idx, level_1.mean, level_2.mean, s);
#endif
                KmerTransitionObservation to = { level_1.mean, level_2.mean, s };
                training_data.kmer_transitions.push_back(to);
            }

            // State-to-state transition
            add_transition_observation(prev_s, s);

            // emission
            //float level = data.read->get_drift_corrected_level(ei, data.strand);
            //float sd = data.read->events[data.strand][ei].stdv;
            //float duration = data.read->get_duration(ei, data.strand);
            if(ki >= n_kmers)
                printf("%zu %d %d %zu %.2lf %c\n", pi, ei, ki, n_kmers, alignment[pi].l_fm, s);
            
            assert(ki < n_kmers);
            //uint32_t rank = sequence.get_kmer_rank(ki, k, data.rc);
        
            //GaussianParameters model = pm.get_scaled_parameters(rank);
            //float norm_level = (level - model.mean) / model.stdv;

            prev_s = s;
        }

        // summary
        training_data.n_matches += (s == 'M');
        training_data.n_merges += (s == 'E');
        training_data.n_skips += (s == 'K');
    }
}

void TransitionParameters::train()
{
    TransitionTrainingData& td = training_data;

    //
    // Profile HMM transitions
    //

    size_t sum_m_not_k = get(td.state_transitions, statechar2index('M'), statechar2index('M')) + 
                         get(td.state_transitions, statechar2index('M'), statechar2index('E'));

    size_t me = get(td.state_transitions, statechar2index('M'), statechar2index('E'));
    double p_me_not_k = (double)me / sum_m_not_k;

    size_t sum_e = 0;
    for(unsigned j = 0; j < td.state_transitions.n_cols; ++j) {
        sum_e += get(td.state_transitions, statechar2index('E'), j);
    }
    
    size_t ee = get(td.state_transitions, statechar2index('E'), statechar2index('E'));
    double p_ee = (double)ee / sum_e;

#ifdef SHOW_TRAINING_RESULT
    fprintf(stderr, "TRANSITIONS\n");
    fprintf(stderr, "M->E|not_k: %lf\n", p_me_not_k);
    fprintf(stderr, "E->E: %lf\n", p_ee);
    for(int i = 0; i < td.state_transitions.n_rows; ++i) {
        fprintf(stderr, "\t%c: ", "MEK"[i]);
        for(int j = 0; j < td.state_transitions.n_cols; ++j) {
            fprintf(stderr, "%d ", get(td.state_transitions, i, j));
        }
        fprintf(stderr, "\n");
    }
#endif

    if(sum_e == 0 || sum_m_not_k == 0) {
        // insufficient data to train, use defaults
        return;
    }

    trans_m_to_e_not_k = p_me_not_k;
    trans_e_to_e = p_ee;

    //
    // Signal-dependent skip probability
    //

    // Initialize observations with pseudocounts from the current model
    size_t num_bins = skip_probabilities.size();
    uint32_t pseudocount = 100;
    std::vector<double> total_observations(num_bins, 0.0f);
    std::vector<double> skip_observations(num_bins, 0.0f);

    for(size_t bin = 0; bin < num_bins; bin++) {
        skip_observations[bin] = skip_probabilities[bin] * pseudocount;
        total_observations[bin] = pseudocount;
    }

    for(size_t oi = 0; oi < td.kmer_transitions.size(); ++oi) {
        const KmerTransitionObservation& to = td.kmer_transitions[oi];
        bool is_skip = to.state == 'K';
        size_t bin = get_skip_bin(to.level_1, to.level_2);

        skip_observations[bin] += is_skip;
        total_observations[bin] += 1;
    }

    // Update probabilities
    for(size_t bin = 0; bin < num_bins; bin++) {
        skip_probabilities[bin] = skip_observations[bin] / total_observations[bin];
#ifdef SHOW_TRAINING_RESULT
        fprintf(stderr, "SKIPLEARN -- %zu %.3lf %.3lf %.3lf\n", bin, skip_observations[bin], total_observations[bin], skip_probabilities[bin]);
#endif
    }
}
