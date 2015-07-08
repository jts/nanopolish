//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_klcs -- function to compute the longest
// common subsequence of k-mers for two strings
//
#ifndef NANOPOLISH_KLCS_H
#define NANOPOLISH_KLCS_H

#include <stdint.h>
#include <vector>
#include <string>
#include "nanopolish_matrix.h"

// The indices of a k-mer match in a pair of sequences
struct kLCSPair
{
    uint32_t i;
    uint32_t j;
};
typedef std::vector<kLCSPair> kLCSResult;

kLCSResult kLCS(const std::string& a, const std::string& b, const int k);

#endif
