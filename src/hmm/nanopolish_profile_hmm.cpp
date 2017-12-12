//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_profile_hmm -- Profile Hidden Markov Model
//
#include <algorithm>
#include "nanopolish_profile_hmm.h"
#include "nanopolish_profile_hmm_r9.h"
#include "nanopolish_profile_hmm_r7.h"

// convenience function to run the HMM over multiple inputs and sum the result
float profile_hmm_score(const HMMInputSequence& sequence, const std::vector<HMMInputData>& data, const uint32_t flags)
{
    float score = 0.0f;
    for(size_t i = 0; i < data.size(); ++i) {
        score += profile_hmm_score(sequence, data[i], flags);
    }
    return score;
}

float profile_hmm_score(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags)
{
    if(data.read->pore_type == PT_R9) {
        return profile_hmm_score_r9(sequence, data, flags);
    } else {
        return profile_hmm_score_r7(sequence, data, flags);
    }
}

float profile_hmm_score_set(const std::vector<HMMInputSequence>& sequences, const HMMInputData& data, const uint32_t flags)
{
    assert(!sequences.empty());
    assert(std::string(sequences[0].get_alphabet()->get_name()) == "nucleotide");
    assert(std::string(data.pore_model->pmalphabet->get_name()) == "nucleotide");
    
    HMMInputData alt_data = data;
    size_t num_models = sequences.size();
    double num_model_penalty = log(num_models);

    // Score the first sequence using the nucleotide alphabet
    double score = profile_hmm_score(sequences[0], data, flags) - num_model_penalty;

    // Score the remainder using the alternative alphabets/pore models
    for(size_t seq_idx = 1; seq_idx < sequences.size(); ++seq_idx) {
        alt_data.pore_model = alt_data.read->get_model(alt_data.strand, sequences[seq_idx].get_alphabet()->get_name());
        assert(alt_data.pore_model != NULL);

        // Score the methylated sequence
        double alt_score = profile_hmm_score(sequences[seq_idx], alt_data, flags) - num_model_penalty;
        score = add_logs(score, alt_score);
    }

    return score;
}

std::vector<HMMAlignmentState> profile_hmm_align(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags)
{
    if(data.read->pore_type == PT_R9) {
        return profile_hmm_align_r9(sequence, data, flags);
    } else {
        return profile_hmm_align_r7(sequence, data, flags);
    }
}
