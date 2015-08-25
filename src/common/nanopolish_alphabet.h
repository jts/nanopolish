//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alphabet -- support for multiple alphabets
//
#ifndef NANOPOLISH_ALPHABET_H
#define NANOPOLISH_ALPHABET_H

#include <inttypes.h>

// A table to map { A, C, G, T } => { 0, 1, 2, 3 }
extern const uint8_t dna_base_rank[];

struct DNABaseMap
{
    static const uint8_t rank[256];
    static const char* base;
    static const char* complement;
    static const uint32_t size;
};

// DNABaseMap with methyl-cytosine
struct MethylCytosineBaseMap
{
    static const uint8_t rank[256];
    static const char* base;
    static const char* complement;
    static const uint32_t size;
};

template<class BaseMap>
class Alphabet
{
    public:
        inline static uint8_t rank(char b) { return BaseMap::rank[b]; }
        inline static char base(uint8_t rank) { return BaseMap::base[rank]; }
        inline static char complement(char b) { return BaseMap::complement[rank(b)]; }
        inline static uint32_t size() { return BaseMap::size; }

        // return the lexicographic rank of the kmer amongst all strings of 
        // length k for this alphabet
        inline static uint32_t kmer_rank(const char* str, uint32_t k)
        {
            uint32_t p = 1;
            uint32_t r = 0;

            // from last base to first
            for(uint32_t i = 0; i < k; ++i) {
                r += rank(str[k - i - 1]) * p;
                p *= BaseMap::size;
            }
            return r;
        }
        
        inline static uint32_t rc_kmer_rank(const char* str, uint32_t k)
        {
            uint32_t p = 1;
            uint32_t r = 0;

            // from first base to last
            for(uint32_t i = 0; i < k; ++i) {
                r += rank(complement(str[i])) * p;
                p *= BaseMap::size;
            }
            return r;
        }

        // Increment the input string to be the next sequence in lexicographic order
        inline static void lexicographic_next(std::string& str)
        {
            int carry = 1;
            int i = str.size() - 1;
            do {
                uint32_t r = rank(str[i]) + carry;
                str[i] = base(r % BaseMap::size);
                carry = r / BaseMap::size;
                i -= 1;
            } while(carry > 0 && i >= 0);
        }

        // Reverse-complement a string
        inline static std::string reverse_complement(const std::string& seq)
        {
            std::string out(seq.length(), 'A');
            size_t last_pos = seq.length() - 1;
            for(int i = last_pos; i >= 0; --i) {
                out[last_pos - i] = complement(seq[i]);
            }
            return out;
        }

};

typedef Alphabet<DNABaseMap> DNAAlphabet;
typedef Alphabet<MethylCytosineBaseMap> MethylCytosineAlphabet;

#endif
