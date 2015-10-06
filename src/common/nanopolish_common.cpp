//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_common -- Data structures and definitions
// shared across files
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

// Split a string into parts based on the delimiter
std::vector<std::string> split(std::string in, char delimiter)
{
    std::vector<std::string> out;
    size_t lastPos = 0;
    size_t pos = in.find_first_of(delimiter);

    while(pos != std::string::npos)
    {
        out.push_back(in.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
        pos = in.find_first_of(delimiter, lastPos);
    }
    out.push_back(in.substr(lastPos));
    return out;
}
