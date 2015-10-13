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
    DNAAlphabet dna_alphabet;
    MethylCpGAlphabet mc_alphabet;

    REQUIRE( dna_alphabet.rank('A') == 0 );
    REQUIRE( dna_alphabet.rank('C') == 1 );
    REQUIRE( dna_alphabet.rank('G') == 2 );
    REQUIRE( dna_alphabet.rank('T') == 3 );

    REQUIRE( dna_alphabet.base(0) == 'A' );
    REQUIRE( dna_alphabet.base(1) == 'C' );
    REQUIRE( dna_alphabet.base(2) == 'G' );
    REQUIRE( dna_alphabet.base(3) == 'T' );

    // MethylCpG alphabet
    REQUIRE( mc_alphabet.rank('A') == 0 );
    REQUIRE( mc_alphabet.rank('C') == 1 );
    REQUIRE( mc_alphabet.rank('G') == 2 );
    REQUIRE( mc_alphabet.rank('M') == 3 );
    REQUIRE( mc_alphabet.rank('T') == 4 );

    REQUIRE( mc_alphabet.base(0) == 'A' );
    REQUIRE( mc_alphabet.base(1) == 'C' );
    REQUIRE( mc_alphabet.base(2) == 'G' );
    REQUIRE( mc_alphabet.base(3) == 'M' );
    REQUIRE( mc_alphabet.base(4) == 'T' );

    // Collectively test lexicographic_next and kmer_rank 
    uint8_t k = 3;
    uint32_t num_strings = pow((double)mc_alphabet.size(), (double)k);
    
    std::string kmer(k, 'A');
    for(size_t i = 0; i < num_strings - 1; ++i) {

        // check lexicographic next
        std::string next = kmer;
        mc_alphabet.lexicographic_next(next);
        REQUIRE( next > kmer );
        int rank_diff = mc_alphabet.kmer_rank(next.c_str(), k) - 
                        mc_alphabet.kmer_rank(kmer.c_str(), k);
        REQUIRE( rank_diff == 1);
        kmer = next;
    }
    REQUIRE(kmer == "TTT");

    // Test the methylate function in the CpG alphabet
    REQUIRE( mc_alphabet.methylate("C") == "C");
    REQUIRE( mc_alphabet.methylate("G") == "G");
    REQUIRE( mc_alphabet.methylate("CG") == "MG");
    REQUIRE( mc_alphabet.methylate("GC") == "GC");
    REQUIRE( mc_alphabet.methylate("CGCG") == "MGMG");
    REQUIRE( mc_alphabet.methylate("AAGCGT") == "AAGMGT");
    REQUIRE( mc_alphabet.methylate("CGGCGT") == "MGGMGT");
    REQUIRE( mc_alphabet.methylate("CGCGC") == "MGMGC");

    // Test unmethylate
    REQUIRE( mc_alphabet.unmethylate("C") == "C");
    REQUIRE( mc_alphabet.unmethylate("CG") == "CG");
    REQUIRE( mc_alphabet.unmethylate("M") == "C");
    REQUIRE( mc_alphabet.unmethylate("MG") == "CG");
    REQUIRE( mc_alphabet.unmethylate("MT") == "MT");

    // Test disambiguate
    REQUIRE( mc_alphabet.disambiguate("") == "");
    REQUIRE( mc_alphabet.disambiguate("M") == "M");
    REQUIRE( mc_alphabet.disambiguate("MT") == "AT");
    REQUIRE( mc_alphabet.disambiguate("MG") == "MG");
    REQUIRE( mc_alphabet.disambiguate("AMG") == "AMG");
    REQUIRE( mc_alphabet.disambiguate("CAM") == "CAM");

    // Test reverse complement
    REQUIRE( mc_alphabet.reverse_complement("M") == "G");
    REQUIRE( mc_alphabet.reverse_complement("C") == "G");
    REQUIRE( mc_alphabet.reverse_complement("MG") == "MG");
    REQUIRE( mc_alphabet.reverse_complement("AM") == "GT");
    REQUIRE( mc_alphabet.reverse_complement("AMG") == "MGT");
    REQUIRE( mc_alphabet.reverse_complement("AAAMG") == "MGTTT");
    REQUIRE( mc_alphabet.reverse_complement("MGMG") == "MGMG");
    REQUIRE( mc_alphabet.reverse_complement("MGAMG") == "MGTMG");
    REQUIRE( mc_alphabet.reverse_complement("GTACATG") == dna_alphabet.reverse_complement("GTACATG"));
}

TEST_CASE( "string functions", "[string_functions]" ) {
    DNAAlphabet dna_alphabet;

    // kmer rank
    REQUIRE( dna_alphabet.kmer_rank("AAAAA", 5) == 0 );
    REQUIRE( dna_alphabet.kmer_rank("GATGA", 5) == 568 );
    REQUIRE( dna_alphabet.kmer_rank("TTTTT", 5) == 1023 );

    // lexicographic increment
    std::string str = "AAAAA";
    dna_alphabet.lexicographic_next(str);
    REQUIRE( str == "AAAAC" );

    str = "AAAAT";
    dna_alphabet.lexicographic_next(str);
    REQUIRE( str == "AAACA" );

    // complement, reverse complement
    REQUIRE( dna_alphabet.reverse_complement("GATGA") == "TCATC" );
}

TEST_CASE( "math", "[math]") {
    GaussianParameters params;
    params.mean = 4;
    params.stdv = 2;
    params.log_stdv = log(params.stdv);

    REQUIRE( normal_pdf(2.25, params) == Approx(0.1360275) );
    REQUIRE( log_normal_pdf(2.25, params) == Approx(log(normal_pdf(2.25, params))) );
}

std::string event_alignment_to_string(const std::vector<HMMAlignmentState>& alignment)
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

    double expected_viterbi_last_state[2] = { -237.7690734863, -266.2348022461 };
    double expected_forward[2] = { -221.1331481934, -262.7491455078 };

    for(int si = 0; si <= 1; ++si) {

        // viterbi align
        std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(ref_subseq, input[si]);
        std::string ea_str = event_alignment_to_string(event_alignment);
    
        // check
        REQUIRE( ea_str == expected_alignment[si]);
        REQUIRE( event_alignment.back().l_fm == Approx(expected_viterbi_last_state[si]));

        // forward algorithm
        double lp = profile_hmm_score(ref_subseq, input[si]);
        REQUIRE(lp == Approx(expected_forward[si]));
    }
}
