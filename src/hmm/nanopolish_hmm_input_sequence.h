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
                             m_alphabet(&gDNAAlphabet),
                             m_seq(seq)
        {
            m_rc_seq = m_alphabet->reverse_complement(seq);
        }
        
        HMMInputSequence(const std::string& fwd,
                         const Alphabet* alphabet) : 
                             m_seq(fwd),
                             m_alphabet(alphabet)
        {
            m_rc_seq = m_alphabet->reverse_complement(m_seq);
        }

        HMMInputSequence(const std::string& fwd,
                         const std::string& rc,
                         const Alphabet* alphabet) : 
                             m_alphabet(alphabet),
                             m_seq(fwd),
                             m_rc_seq(rc)
        {

        }

        //
        const std::string& get_sequence() const { return m_seq; }

        //
        const Alphabet* get_alphabet() const { return m_alphabet; }

        //
        size_t length() const { return m_seq.length(); }

        // swap sequence and its reverse complement
        void swap() { m_seq.swap(m_rc_seq); }

        // returns the i-th kmer of the sequence
        inline std::string get_kmer(uint32_t i, uint32_t k, bool do_rc) const
        {
            return ! do_rc ? m_seq.substr(i, k) : 
                             m_rc_seq.substr(m_rc_seq.length() - i - k, k);
        }
        
        // get the number of kmer ranks supported by the alphabet for this sequence
        size_t get_num_kmer_ranks(size_t k) const { return m_alphabet->get_num_strings(k); }

        // get the lexicographic rank of the i-th kmer
        // if the do_rc flag is set, return the rank of
        // reverse-complemented version of the ki-th kmer
        // NOT the ki-th kmer of the reverse-complemented sequence
        inline uint32_t get_kmer_rank(uint32_t i, uint32_t k, bool do_rc) const
        {
            return ! do_rc ? _kmer_rank(i, k) : _rc_kmer_rank(i, k);
        }

    private:

        inline uint32_t _kmer_rank(uint32_t i, uint32_t k) const
        {
            return m_alphabet->kmer_rank(m_seq.c_str() + i, k);
        }

        inline uint32_t _rc_kmer_rank(uint32_t i, uint32_t k) const
        {
            return m_alphabet->kmer_rank(m_rc_seq.c_str() + (length() - i - k), k);
        }

        HMMInputSequence(); // not allowed
        const Alphabet* m_alphabet;

        std::string m_seq;
        std::string m_rc_seq;
};

#endif
