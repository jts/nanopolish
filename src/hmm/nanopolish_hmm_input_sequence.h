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

//
// This class is a general wrapper around a string
// that allows different alphabets to be used.
// 
class HMMInputSequence
{
    public:
    
        // constructors
        HMMInputSequence(const std::string& seq) : 
                             m_seq(seq),
                             m_alphabet(&gDNAAlphabet)
        {
            m_rc_seq = m_alphabet->reverse_complement(seq);
        }

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
            return ! do_rc ? _kmer_rank(i) : _rc_kmer_rank(i);
        }

    private:

        inline uint32_t _kmer_rank(uint32_t i) const
        {
            return m_alphabet->kmer_rank(m_seq.c_str() + i, K);
        }

        inline uint32_t _rc_kmer_rank(uint32_t i) const
        {
            return m_alphabet->kmer_rank(m_rc_seq.c_str() + (length() - i - K), K);
        }

        HMMInputSequence(); // not allowed
        Alphabet* m_alphabet;

        std::string m_seq;
        std::string m_rc_seq;
};

#endif
