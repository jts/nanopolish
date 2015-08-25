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
#include "nanopolish_alphabet.h"
#include "nanopolish_emissions.h"
#include "nanopolish_profile_hmm.h"

// This code needs to be run before any of the program logic
// It sets up pre-computed values and caches
void initialize()
{
    p7_FLogsumInit();
}

TEST_CASE( "alphabet", "[alphabet]" ) {

    // DNA alphabet
    REQUIRE( DNAAlphabet::rank('A') == 0 );
    REQUIRE( DNAAlphabet::rank('C') == 1 );
    REQUIRE( DNAAlphabet::rank('G') == 2 );
    REQUIRE( DNAAlphabet::rank('T') == 3 );

    REQUIRE( DNAAlphabet::base(0) == 'A' );
    REQUIRE( DNAAlphabet::base(1) == 'C' );
    REQUIRE( DNAAlphabet::base(2) == 'G' );
    REQUIRE( DNAAlphabet::base(3) == 'T' );

    // MethylCytosine alphabet
    REQUIRE( MethylCytosineAlphabet::rank('A') == 0 );
    REQUIRE( MethylCytosineAlphabet::rank('C') == 1 );
    REQUIRE( MethylCytosineAlphabet::rank('G') == 2 );
    REQUIRE( MethylCytosineAlphabet::rank('M') == 3 );
    REQUIRE( MethylCytosineAlphabet::rank('T') == 4 );

    REQUIRE( MethylCytosineAlphabet::base(0) == 'A' );
    REQUIRE( MethylCytosineAlphabet::base(1) == 'C' );
    REQUIRE( MethylCytosineAlphabet::base(2) == 'G' );
    REQUIRE( MethylCytosineAlphabet::base(3) == 'M' );
    REQUIRE( MethylCytosineAlphabet::base(4) == 'T' );

    // Collectively test lexicographic_next, kmer_rank, rc_kmer_rank
    // and reverse_complement
    uint8_t k = 3;
    uint32_t num_strings = pow((double)MethylCytosineAlphabet::size(), (double)k);
    
    std::string kmer(k, 'A');
    for(size_t i = 0; i < num_strings - 1; ++i) {

        // check that the rank(rc(str)) is equal to rc_rank(str)
        std::string rc = MethylCytosineAlphabet::reverse_complement(kmer);
        REQUIRE( MethylCytosineAlphabet::kmer_rank(rc.c_str(), k) == 
                 MethylCytosineAlphabet::rc_kmer_rank(kmer.c_str(), k) );

        // check lexicographic next
        std::string next = kmer;
        MethylCytosineAlphabet::lexicographic_next(next);
        REQUIRE( next > kmer );
        int rank_diff = MethylCytosineAlphabet::kmer_rank(next.c_str(), k) - 
                        MethylCytosineAlphabet::kmer_rank(kmer.c_str(), k);
        REQUIRE( rank_diff == 1);
        kmer = next;
    }
    REQUIRE(kmer == "TTT");
}

TEST_CASE( "string functions", "[string_functions]" ) {

    // kmer rank
    REQUIRE( DNAAlphabet::kmer_rank("AAAAA", 5) == 0 );
    REQUIRE( DNAAlphabet::kmer_rank("GATGA", 5) == 568 );
    REQUIRE( DNAAlphabet::kmer_rank("TTTTT", 5) == 1023 );
    REQUIRE( DNAAlphabet::kmer_rank("GATGA", 5) == DNAAlphabet::rc_kmer_rank("TCATC", 5 ) );

    // lexicographic increment
    std::string str = "AAAAA";
    DNAAlphabet::lexicographic_next(str);
    REQUIRE( str == "AAAAC" );

    str = "AAAAT";
    DNAAlphabet::lexicographic_next(str);
    REQUIRE( str == "AAACA" );

    // complement, reverse complement
    REQUIRE( DNAAlphabet::complement('A') == 'T' );
    REQUIRE( DNAAlphabet::reverse_complement("GATGA") == "TCATC" );
}

TEST_CASE( "math", "[math]") {
    GaussianParameters params;
    params.mean = 4;
    params.stdv = 2;
    params.log_stdv = log(params.stdv);

    REQUIRE( normal_pdf(2.25, params) == Approx(0.1360275) );
    REQUIRE( log_normal_pdf(2.25, params) == Approx(log(normal_pdf(2.25, params))) );
}

std::string event_alignment_to_string(const std::vector<AlignmentState>& alignment)
{
    std::string out;
    for(size_t i = 0; i < alignment.size(); ++i) {
        out.append(1, alignment[i].state);
    }
    return out;
}

TEST_CASE( "hmm", "[hmm]") {

    // read the FAST5
    SquiggleRead sr("test_read", "test/data/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch8_file30_strand.fast5");
    sr.transform();

    // The reference sequence to align to:
    std::string ref_subseq = "ATCAGTAAAATAACGTAGAGCGGTAACCTTGCCATAAAGGTCGAGTTTA"
                             "TTACCATCCTTGTTATAGACTTCGGCAGCGTGTGCTACGTTCGCAGCT";

    // Generate a HMMInputData structure to tell the HMM
    // which part of the read to align
    HMMInputData input[2];
    
    // template strand
    input[0].read = &sr;
    input[0].event_start_idx = 3;
    input[0].event_stop_idx = 88;
    input[0].event_stride = 1;
    input[0].rc = false;
    input[0].strand = 0;
    
    // complement strand
    input[1].read = &sr;
    input[1].event_start_idx = 6788;
    input[1].event_stop_idx = 6697;
    input[1].event_stride = -1;
    input[1].rc = true;
    input[1].strand = 1;

    // expected output
    std::string expected_alignment[2];
    expected_alignment[0] = 
        "MMMMMEMKMKMMMMMMMKMMMKMMMKMMMMMMMMMKKMMEEEMMMMMMKMMMM" 
        "MMMKMMMMMKMKMKMEMKKMKMKKMMMMMMEMMMMKMKMEEMMMMKMEEEEEM";

    expected_alignment[1] = 
        "MMKMMMKMEEMMKMKMKMEMMMKMMMKMEMMMKMMMKMMMMMMMMMKKMEMMMM"
        "EMMMMMMMMKMKKMMMMMMMEMMMMMKMMMMMKMEMMMMMKMMMMMEEEEEEEEM";

    double expected_viterbi_last_state[2] = { -237.7690734863, -266.1903076172 };
    double expected_forward[2] = { -221.1331481934, -262.7046203613 };

    for(int si = 0; si <= 1; ++si) {

        // viterbi align
        std::vector<AlignmentState> event_alignment = profile_hmm_align(ref_subseq, input[si]);
        std::string ea_str = event_alignment_to_string(event_alignment);
    
        // check
        REQUIRE( ea_str == expected_alignment[si]);
        REQUIRE( event_alignment.back().l_fm == Approx(expected_viterbi_last_state[si]));

        // forward algorithm
        double lp = profile_hmm_score(ref_subseq, input[si]);
        REQUIRE(lp == Approx(expected_forward[si]));
    }
}
