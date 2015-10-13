//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alphabet -- support for multiple alphabets
//
#include <cassert>
#include <vector>
#include "nanopolish_alphabet.h"

//
// DNAAlphabet
// 

const uint8_t DNAAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const char* DNAAlphabet::_base = "ACGT";
const char* DNAAlphabet::_complement = "TGCA";
const uint32_t DNAAlphabet::_size = 4;

//
// methyl-cytosine in CG context
//
const uint8_t MethylCpGAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylCpGAlphabet::_base = "ACGMT";
const char* MethylCpGAlphabet::_complement = "TGCGA";
const uint32_t MethylCpGAlphabet::_size = 5;

const uint32_t MethylCpGAlphabet::_num_recognition_sites = 1;
const uint32_t MethylCpGAlphabet::_recognition_length = 2;
const char* MethylCpGAlphabet::_recognition_sites[] = { "CG" };
const char* MethylCpGAlphabet::_recognition_sites_methylated[] = { "MG" };
const char* MethylCpGAlphabet::_recognition_sites_methylated_complement[] = { "GM" };

//
// Dam methylation: methyl-adenine in GATC context
//
const uint8_t MethylDamAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylDamAlphabet::_base = "ACGMT";
const char* MethylDamAlphabet::_complement = "TGCTA";
const uint32_t MethylDamAlphabet::_size = 5;

const uint32_t MethylDamAlphabet::_num_recognition_sites = 1;
const uint32_t MethylDamAlphabet::_recognition_length = 4;
const char* MethylDamAlphabet::_recognition_sites[] = { "GATC" };
const char* MethylDamAlphabet::_recognition_sites_methylated[] = { "GMTC" };
const char* MethylDamAlphabet::_recognition_sites_methylated_complement[] = { "CTMG" };

//
// Dcm methylation: methyl-cytosine in CCAGG and CCTGG context
//
const uint8_t MethylDcmAlphabet::_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,3,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

const char* MethylDcmAlphabet::_base = "ACGMT";
const char* MethylDcmAlphabet::_complement = "TGCGA";
const uint32_t MethylDcmAlphabet::_size = 5;

const uint32_t MethylDcmAlphabet::_num_recognition_sites = 2;
const uint32_t MethylDcmAlphabet::_recognition_length = 5;
const char* MethylDcmAlphabet::_recognition_sites[] = { "CCAGG", "CCTGG" };
const char* MethylDcmAlphabet::_recognition_sites_methylated[] = { "CMAGG", "CMTGG" };
const char* MethylDcmAlphabet::_recognition_sites_methylated_complement[] = { "GGTMC", "GGAMC" };

// Global objects
DNAAlphabet gDNAAlphabet;
MethylCpGAlphabet gMCpGAlphabet;
MethylDamAlphabet gMethylDamAlphabet;

// Select the alphabet that best matches bases
const Alphabet *best_alphabet(const char *bases)
{
    std::vector<const Alphabet *>list = {&gDNAAlphabet, &gMCpGAlphabet};

    for (auto alphabet: list)
        if (alphabet->contains_all(bases))
            return alphabet;

    return nullptr;                
}

