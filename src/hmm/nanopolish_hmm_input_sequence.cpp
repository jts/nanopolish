//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_hmm_input_sequence -- a nucleotide sequence
// that is input into the hidden Markov model
//
#include "nanopolish_hmm_input_sequence.h"

HMMInputSequence::HMMInputSequence(const std::string& seq) : m_seq(seq)
{

}
