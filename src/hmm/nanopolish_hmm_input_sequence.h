//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_hmm_input_sequence -- a nucleotide sequence
// that is input into the hidden Markov model
//
#ifndef NANOPOLISH_HMM_INPUT_SEQUENCE
#define NANOPOLISH_HMM_INPUT_SEQUENCE

#include <string>
#include "nanopolish_common.h"
#include "nanopolish_alphabet.h"

// This class is a general wrapper around a string
// that allows different alphabets to be used.
// 
class HMMInputSequence
{
    public:
    
        // constructors
        HMMInputSequence(const std::string& seq);

        //
        size_t length() const { return m_seq.length(); }

        // returns the i-th kmer of the sequence
        inline std::string get_kmer(uint32_t i, uint32_t k) const
        {
            return m_seq.substr(i, k);
        }

        // get the lexicographic rank of the i-th kmer
        // if the do_rc flag is set, return the rank of
        // reverse-complemented version of the ki-th kmer
        // NOT the ki-th kmer of the reverse-complemented sequence
        inline uint32_t get_kmer_rank(uint32_t i, bool do_rc) const
        {
            const char* p = m_seq.c_str() + i;
            return ! do_rc ?  DNAAlphabet::kmer_rank(p, K) : DNAAlphabet::rc_kmer_rank(p, K);
        }

    private:
        HMMInputSequence(); // not allowed

        std::string m_seq;
};

#endif
