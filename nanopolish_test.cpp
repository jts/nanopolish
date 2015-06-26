//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_test -- test driver program
//
#define CATCH_CONFIG_MAIN
#include <stdio.h>
#include <string>
#include "logsum.h"
#include "catch.hpp"

// This code needs to be run before any of the program logic
// It sets up pre-computed values and caches
void initialize()
{
    p7_FLogsumInit();
}

int sum(int a, int b)
{
    return a + b;
}

TEST_CASE( "Basic sums", "[sum]" ) {
    REQUIRE( sum(0, 1) == 1 ); // pass
    REQUIRE( sum(10, 11) == 20 ); // fail
    REQUIRE( sum(2, 3) == 5 ); // pass
}
