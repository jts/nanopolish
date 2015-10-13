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
#include <cstring>
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
        virtual char complement(char b) const = 0;
        virtual uint32_t size() const = 0;

        // support for methylated bases with recognition sequences
        virtual size_t num_recognition_sites() const = 0;
        virtual size_t recognition_length() const = 0;
        virtual const char* get_recognition_site(size_t i) const = 0;
        virtual const char* get_recognition_site_methylated(size_t i) const = 0;
        virtual const char* get_recognition_site_methylated_complement(size_t i) const = 0;

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

        // returns the number of unique strings of length l for this alphabet
        inline size_t get_num_strings(size_t l) const
        {
            size_t s = size();
            size_t n = 1;
            for(size_t i = 0; i < l; ++i) {
                n *= s;
            }
            return n;
        }

        virtual std::string reverse_complement(const std::string& str) const
        {
            std::string out(str.length(), 'A');
            size_t i = 0; // input
            int j = str.length() - 1; // output
            while(i < str.length()) {
                int recognition_index = -1;
                size_t match_length = 0;

                // Does this location (partially) match a methylated recognition site?
                for(size_t j = 0; j < num_recognition_sites(); ++j) {

                    // Require the recognition site to be completely matched
                    size_t cl = std::min((size_t)recognition_length(), str.length() - i);
                    if(str.compare(i, cl, get_recognition_site_methylated(j), cl) == 0) {
                        
                        // Matches a recognition site
                        recognition_index = j;
                        match_length = cl;
                        break;
                    }
                }

                // If this subsequence matched a methylated recognition site,
                // copy the complement of the site to the output
                if(recognition_index != -1) {
                    for(size_t k = 0; k < match_length; ++k) {
                        out[j--] = get_recognition_site_methylated_complement(recognition_index)[k];
                        i += 1;
                    }
                } else {
                    // complement a single base
                    assert(str[i] != 'M');
                    out[j--] = complement(str[i++]);
                }
            }
            return out;
        }

        // return a new copy of the string with IUPAC ambiguity characters changed
        virtual std::string disambiguate(const std::string& str) const
        {
            std::string out(str);
            size_t i = 0;
            while(i < out.length()) {
                size_t stride = 1;
                bool is_recognition_site = false;

                // Does this location (partially) match a methylated recognition site?
                for(size_t j = 0; j < num_recognition_sites(); ++j) {

                    // Require the recognition site to be completely matched
                    size_t cl = std::min((size_t)recognition_length(), out.length() - i);
                    if(out.compare(i, cl, get_recognition_site_methylated(j), cl) == 0) {
                        stride = cl; // skip to end of match
                        is_recognition_site = true;
                        break;
                    }
                }
                
                // disambiguate if not a recognition site
                if(!is_recognition_site) {
                    assert(IUPAC::isValid(out[i]));
                    out[i] = IUPAC::getPossibleSymbols(out[i])[0];
                    stride = 1;
                }

                i += stride;
            }
            return out;
        }

        // If the alphabet supports methylated bases, convert str
        // to a methylated string using the recognition sites
        virtual std::string methylate(const std::string& str) const
        {
            std::string out(str);
            size_t i = 0;
            while(i < out.length()) {
                size_t stride = 1;

                // Does this location match a recognition site?
                for(size_t j = 0; j < num_recognition_sites(); ++j) {

                    // Require the recognition site to be completely matched
                    if(str.compare(i, recognition_length(), get_recognition_site(j)) == 0) {

                        // Replace by the methylated version
                        out.replace(i, recognition_length(), get_recognition_site_methylated(j));
                        stride = recognition_length(); // skip to end of match
                        break;
                    }
                }

                i += stride;
            }
            return out;
        }

        // Remove methylated bases according to the recognition site
        std::string unmethylate(const std::string& str) const
        {
            std::string out(str);
            size_t i = 0;
            while(i < out.length()) {
                size_t stride = 1;

                // Does this location (partially) match a methylated recognition site?
                for(size_t j = 0; j < num_recognition_sites(); ++j) {

                    // Require the recognition site to be completely matched
                    size_t cl = std::min((size_t)recognition_length(), out.length() - i);
                    if(out.compare(i, cl, get_recognition_site_methylated(j), cl) == 0) {

                        // Replace by the unmethylated version
                        out.replace(i, cl, get_recognition_site(j), cl);
                        stride = cl; // skip to end of match
                        break;
                    }
                }

                i += stride;
            }
            return out;
        }

        // does this alphabet contain all of the nucleotides in bases?
        virtual inline bool contains_all(const char *bases) const = 0;
};

struct DNAAlphabet : public Alphabet
{
    // members
    static const uint8_t _rank[256];
    static const char* _base;
    static const char* _complement;
    static const uint32_t _size;

    // functions
    virtual uint8_t rank(char b) const { return _rank[b]; }
    virtual char base(uint8_t r) const { return _base[r]; }
    virtual char complement(char b) const { return _complement[_rank[b]]; }
    virtual uint32_t size() const { return _size; }

    // no methylation in this alphabet
    virtual size_t num_recognition_sites() const { return 0; }
    virtual size_t recognition_length() const { return 0; }
    virtual const char* get_recognition_site(size_t i) const { return NULL; }
    virtual const char* get_recognition_site_methylated(size_t i) const { return NULL; }
    virtual const char* get_recognition_site_methylated_complement(size_t i) const { 
        return NULL;
    }

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

// DNABaseMap with methyl-cytosine
struct MethylCpGAlphabet : public Alphabet
{
    static const uint8_t _rank[256];
    static const char* _base;
    static const char* _complement;
    static const uint32_t _size;

    // methylation support
    static const uint32_t _num_recognition_sites;
    static const uint32_t _recognition_length;
    static const char* _recognition_sites[];
    static const char* _recognition_sites_methylated[];
    static const char* _recognition_sites_methylated_complement[];

    virtual uint8_t rank(char b) const { return _rank[b]; }
    virtual char base(uint8_t r) const { return _base[r]; }
    virtual char complement(char b) const { return _complement[_rank[b]]; }
    virtual uint32_t size() const { return _size; }

    virtual size_t num_recognition_sites() const { return _num_recognition_sites; }
    virtual size_t recognition_length() const { return _recognition_length; }
    virtual const char* get_recognition_site(size_t i) const { return _recognition_sites[i]; }
    virtual const char* get_recognition_site_methylated(size_t i) const { return _recognition_sites_methylated[i]; }
    virtual const char* get_recognition_site_methylated_complement(size_t i) const { 
        return _recognition_sites_methylated_complement[i]; 
    }

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const 
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

// Global alphabet objects that can be re-used
extern DNAAlphabet gDNAAlphabet;
extern MethylCpGAlphabet gMCpGAlphabet;

const Alphabet *best_alphabet(const char *bases);
#endif
