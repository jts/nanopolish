//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_iupac -- handle iupac ambiguity codes
//
#ifndef NANOPOLISH_IUPAC_H
#define NANOPOLISH_IUPAC_H

#include <string>

// IUPAC ambiguity alphabet
namespace IUPAC
{
    // Returns true if c is [ACGT]
    bool isUnambiguous(char c);

    // Returns true if c is a valid ambiguity code
    bool isAmbiguous(char c);

    // Returns true if c is a valid symbol in this alphabet
    bool isValid(char c);

    // Returns a string defining the possible unambiguous bases for each symbol
    // in the alphabet
    std::string getPossibleSymbols(char c);
}

#endif
