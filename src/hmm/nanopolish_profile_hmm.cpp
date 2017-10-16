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
    if(data.read->pore_model_metadata[data.strand].is_r9()) {
        return profile_hmm_score_r9(sequence, data, flags);
    } else {
        return profile_hmm_score_r7(sequence, data, flags);
    }
}

std::vector<HMMAlignmentState> profile_hmm_align(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags)
{
    if(data.read->pore_model_metadata[data.strand].is_r9()) {
        return profile_hmm_align_r9(sequence, data, flags);
    } else {
        return profile_hmm_align_r7(sequence, data, flags);
    }
}
