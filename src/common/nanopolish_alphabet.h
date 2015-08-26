//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alphabet -- support for multiple alphabets
//
#ifndef NANOPOLISH_ALPHABET_H
#define NANOPOLISH_ALPHABET_H

#include <string>
#include <inttypes.h>
#include <assert.h>
#include "nanopolish_iupac.h"

// A table to map { A, C, G, T } => { 0, 1, 2, 3 }
extern const uint8_t dna_base_rank[];

// Abstract base class for alphabets
class Alphabet
{
    public:
        // basic functions
        virtual uint8_t rank(char b) const = 0;
        virtual char base(uint8_t r) const = 0;
        virtual uint32_t size() const = 0;

        // return the lexicographic rank of the kmer amongst all strings of 
        // length k for this alphabet
        inline uint32_t kmer_rank(const char* str, uint32_t k) const
        {
            uint32_t p = 1;
            uint32_t r = 0;

            // from last base to first
            for(uint32_t i = 0; i < k; ++i) {
                r += rank(str[k - i - 1]) * p;
                p *= size();
            }
            return r;
        }
        
        // Increment the input string to be the next sequence in lexicographic order
        inline void lexicographic_next(std::string& str) const
        {
            int carry = 1;
            int i = str.size() - 1;
            do {
                uint32_t r = rank(str[i]) + carry;
                str[i] = base(r % size());
                carry = r / size();
                i -= 1;
            } while(carry > 0 && i >= 0);
        }

        // reverse complement a string over this alphabet
        virtual std::string reverse_complement(const std::string& seq) const = 0;

        // remove ambiguous nucleotides from the string
        virtual std::string disambiguate(const std::string& seq) const = 0; 
};

struct DNAAlphabet : public Alphabet
{
    static const uint8_t _rank[256];
    static const char* _base;
    static const char* _complement;
    static const uint32_t _size;

    virtual uint8_t rank(char b) const { return _rank[b]; }
    virtual char base(uint8_t r) const { return _base[r]; }
    virtual uint32_t size() const { return _size; }

    virtual std::string reverse_complement(const std::string& seq) const
    {
        std::string out(seq.length(), 'A');
        size_t last_pos = seq.length() - 1;
        for(int i = last_pos; i >= 0; --i) {
            out[last_pos - i] = _complement[_rank[seq[i]]];
        }
        return out;
    }

    // return a new copy of the string with ambiguous characters changed
    virtual std::string disambiguate(const std::string& str) const
    {
        std::string out(str);
        for(size_t i = 0; i < str.length(); ++i) {
            assert(IUPAC::isValid(str[i]));
            out[i] = IUPAC::getPossibleSymbols(str[i])[0];
        }
        return out;
    }

};

// DNABaseMap with methyl-cytosine
struct MethylCpGAlphabet : public Alphabet
{
    static const uint8_t _rank[256];
    static const char* _base;
    static const char* _complement;
    static const uint32_t _size;

    virtual uint8_t rank(char b) const { return _rank[b]; }
    virtual char base(uint8_t r) const { return _base[r]; }
    virtual uint32_t size() const { return _size; }

    virtual std::string reverse_complement(const std::string& seq) const
    {
        std::string out(seq.length(), 'A');
        size_t i = 0; // input
        int j = seq.length() - 1; // output
        while(i < seq.length()) {
            if(seq[i] == 'M') {
                
                out[j--] = 'G';
                i += 1;

                // CpG methylation model requires M to be followed by G
                // (if there is space)
                if(j >= 0) {
                    assert(i < seq.length());
                    assert(seq[i] == 'G');
                    out[j--] = 'M';
                    ++i;
                }
            } else {
                out[j--] = DNAAlphabet::_complement[DNAAlphabet::_rank[seq[i++]]];
            }
        }
        return out;
    }

    // return a new copy of the string with ambiguous characters changed
    virtual std::string disambiguate(const std::string& str) const
    {
        std::string out(str);
        for(size_t i = 0; i < str.length(); ++i) {
            if(str[i] == 'M' && i != str.length() - 1 && str[i + 1] == 'G') {
                // CpG site, assume its methylated not an ambiguity symbol
                out[i] = 'M';
            } else {
                assert(IUPAC::isValid(str[i]));
                out[i] = IUPAC::getPossibleSymbols(str[i])[0];
            }
        }
        return out;
    }

};

// Global alphabet objects that can be re-used
extern DNAAlphabet gDNAAlphabet;
extern MethylCpGAlphabet gMCpGAlphabet;

#endif
