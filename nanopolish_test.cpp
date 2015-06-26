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
#include "nanopolish_common.h"
#include "nanopolish_emissions.h"

// This code needs to be run before any of the program logic
// It sets up pre-computed values and caches
void initialize()
{
    p7_FLogsumInit();
}

TEST_CASE( "string functions", "[string_functions]" ) {

    // kmer rank
    REQUIRE( kmer_rank("AAAAA", 5) == 0 );
    REQUIRE( kmer_rank("GATGA", 5) == 568 );
    REQUIRE( kmer_rank("TTTTT", 5) == 1023 );
    REQUIRE( kmer_rank("GATGA", 5) == rc_kmer_rank("TCATC", 5 ) );

    // lexicographic increment
    std::string str = "AAAAA";
    lexicographic_next(str);
    REQUIRE( str == "AAAAC" );

    str = "AAAAT";
    lexicographic_next(str);
    REQUIRE( str == "AAACA" );

    // complement, reverse complement
    REQUIRE( complement('A') == 'T' );
    REQUIRE( reverse_complement("GATGA") == "TCATC" );

}

TEST_CASE( "math", "[math]") {
    GaussianParameters params;
    params.mean = 4;
    params.stdv = 2;
    params.log_stdv = log(params.stdv);

    REQUIRE( normal_pdf(2.25, params) == Approx(0.1360275) );
    REQUIRE( log_normal_pdf(2.25, params) == Approx(log(normal_pdf(2.25, params))) );
}
