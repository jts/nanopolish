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
#include <array>
#include <vector>
#include <random>

#include "logsum.h"
#include "catch.hpp"
#include "nanopolish_common.h"
#include "nanopolish_alphabet.h"
#include "nanopolish_emissions.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_variant_db.h"
#include "training_core.hpp"
#include "invgauss.hpp"
#include "logger.hpp"

TEST_CASE( "alphabet", "[alphabet]" ) {

    // DNA alphabet
    DNAAlphabet dna_alphabet;
    MethylCpGAlphabet mc_alphabet;
    MethylDamAlphabet dam_alphabet;
    MethylDcmAlphabet dcm_alphabet;

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

    // unmethylate
    REQUIRE( mc_alphabet.unmethylate("C") == "C");
    REQUIRE( mc_alphabet.unmethylate("CG") == "CG");
    REQUIRE( mc_alphabet.unmethylate("M") == "C");
    REQUIRE( mc_alphabet.unmethylate("MG") == "CG");
    REQUIRE( mc_alphabet.unmethylate("MT") == "MT");

    // disambiguate
    REQUIRE( mc_alphabet.disambiguate("") == "");
    REQUIRE( mc_alphabet.disambiguate("M") == "M");
    REQUIRE( mc_alphabet.disambiguate("MT") == "AT");
    REQUIRE( mc_alphabet.disambiguate("MG") == "MG");
    REQUIRE( mc_alphabet.disambiguate("AMG") == "AMG");
    REQUIRE( mc_alphabet.disambiguate("CAM") == "CAM");

    // reverse complement
    REQUIRE( mc_alphabet.reverse_complement("M") == "G");
    REQUIRE( mc_alphabet.reverse_complement("C") == "G");
    REQUIRE( mc_alphabet.reverse_complement("G") == "C");
    REQUIRE( mc_alphabet.reverse_complement("MG") == "MG");
    REQUIRE( mc_alphabet.reverse_complement("CG") == "CG");
    REQUIRE( mc_alphabet.reverse_complement("AM") == "GT");
    REQUIRE( mc_alphabet.reverse_complement("AMG") == "MGT");
    REQUIRE( mc_alphabet.reverse_complement("AAAMG") == "MGTTT");
    REQUIRE( mc_alphabet.reverse_complement("MGMG") == "MGMG");
    REQUIRE( mc_alphabet.reverse_complement("MGAMG") == "MGTMG");
    REQUIRE( mc_alphabet.reverse_complement("GTACATG") == dna_alphabet.reverse_complement("GTACATG"));

    // Dam methylation tests
    
    // methylate
    REQUIRE( dam_alphabet.methylate("") == "");
    REQUIRE( dam_alphabet.methylate("G") == "G");
    REQUIRE( dam_alphabet.methylate("GA") == "GA");
    REQUIRE( dam_alphabet.methylate("GAT") == "GAT");
    REQUIRE( dam_alphabet.methylate("GATC") == "GMTC");
    REQUIRE( dam_alphabet.methylate("GATCG") == "GMTCG");
    REQUIRE( dam_alphabet.methylate("GATCGA") == "GMTCGA");
    REQUIRE( dam_alphabet.methylate("GATCGAT") == "GMTCGAT");
    REQUIRE( dam_alphabet.methylate("GATCGATC") == "GMTCGMTC");
    REQUIRE( dam_alphabet.methylate("GMTCGATC") == "GMTCGMTC");
    REQUIRE( dam_alphabet.methylate("GMTCGMTC") == "GMTCGMTC");

    // unmethylate
    REQUIRE( dam_alphabet.unmethylate("M") == "A");
    REQUIRE( dam_alphabet.unmethylate("MT") == "AT");
    REQUIRE( dam_alphabet.unmethylate("MTC") == "ATC");
    REQUIRE( dam_alphabet.unmethylate("GM") == "GA");
    REQUIRE( dam_alphabet.unmethylate("GMT") == "GAT");
    REQUIRE( dam_alphabet.unmethylate("GMTC") == "GATC");
    REQUIRE( dam_alphabet.unmethylate("GMTCG") == "GATCG");
    REQUIRE( dam_alphabet.unmethylate("GMTCGM") == "GATCGA");
    REQUIRE( dam_alphabet.unmethylate("GMTCGMTC") == "GATCGATC");
    REQUIRE( dam_alphabet.unmethylate("GMTCGMT") == "GATCGAT");
    REQUIRE( dam_alphabet.unmethylate("GMTCGM") == "GATCGA");
    REQUIRE( dam_alphabet.unmethylate("MA") == "MA");
    REQUIRE( dam_alphabet.unmethylate("MT") == "AT");
    REQUIRE( dam_alphabet.unmethylate("GM") == "GA");
    REQUIRE( dam_alphabet.unmethylate("CM") == "CM");

    // disambiguate
    REQUIRE( dam_alphabet.disambiguate("") == "");
    REQUIRE( dam_alphabet.disambiguate("GMTC") == "GMTC");
    REQUIRE( dam_alphabet.disambiguate("M") == "M");
    REQUIRE( dam_alphabet.disambiguate("MT") == "MT");
    REQUIRE( dam_alphabet.disambiguate("MTC") == "MTC");
    REQUIRE( dam_alphabet.disambiguate("GM") == "GM");
    REQUIRE( dam_alphabet.disambiguate("GMT") == "GMT");
    REQUIRE( dam_alphabet.disambiguate("GMA") == "GAA");

    // reverse complement
    REQUIRE( dam_alphabet.reverse_complement("") == "");
    REQUIRE( dam_alphabet.reverse_complement("M") == "T");
    REQUIRE( dam_alphabet.reverse_complement("G") == "C");
    REQUIRE( dam_alphabet.reverse_complement("GM") == "TC");
    REQUIRE( dam_alphabet.reverse_complement("GMT") == "MTC");
    REQUIRE( dam_alphabet.reverse_complement("GMTC") == "GMTC");
    REQUIRE( dam_alphabet.reverse_complement("MTC") == "GMT");
    REQUIRE( dam_alphabet.reverse_complement("TC") == "GA");
    REQUIRE( dam_alphabet.reverse_complement("C") == "G");
    REQUIRE( dam_alphabet.reverse_complement("GATC") == "GATC");
    REQUIRE( dam_alphabet.reverse_complement("ATC") == "GAT");
    REQUIRE( dam_alphabet.reverse_complement("TC") == "GA");
    REQUIRE( dam_alphabet.reverse_complement("GAT") == "ATC");

    //
    // Dcm methylation tests
    //

    // methylate
    REQUIRE( dcm_alphabet.methylate("") == "");
    REQUIRE( dcm_alphabet.methylate("C") == "C");
    REQUIRE( dcm_alphabet.methylate("CC") == "CC");

    // first recognition site
    REQUIRE( dcm_alphabet.methylate("CCA") == "CCA");
    REQUIRE( dcm_alphabet.methylate("CCAG") == "CCAG");
    REQUIRE( dcm_alphabet.methylate("CCAGG") == "CMAGG");
    REQUIRE( dcm_alphabet.methylate("CAGG") == "CAGG");
    REQUIRE( dcm_alphabet.methylate("AGG") == "AGG");
    
    // second recognition site
    REQUIRE( dcm_alphabet.methylate("CCT") == "CCT");
    REQUIRE( dcm_alphabet.methylate("CCTG") == "CCTG");
    REQUIRE( dcm_alphabet.methylate("CCTGG") == "CMTGG");
    REQUIRE( dcm_alphabet.methylate("CTGG") == "CTGG");
    REQUIRE( dcm_alphabet.methylate("TGG") == "TGG");

    // both recognition sites
    REQUIRE( dcm_alphabet.methylate("CCAGGCCTGG") == "CMAGGCMTGG");
    REQUIRE( dcm_alphabet.methylate("CCAGGCCTG") == "CMAGGCCTG");

    // unmethylate
    REQUIRE( dcm_alphabet.unmethylate("M") == "C");
    REQUIRE( dcm_alphabet.unmethylate("MA") == "CA");
    REQUIRE( dcm_alphabet.unmethylate("MT") == "CT");
    REQUIRE( dcm_alphabet.unmethylate("MAG") == "CAG");
    REQUIRE( dcm_alphabet.unmethylate("MTG") == "CTG");
    REQUIRE( dcm_alphabet.unmethylate("MAGG") == "CAGG");
    REQUIRE( dcm_alphabet.unmethylate("MTGG") == "CTGG");
    
    REQUIRE( dcm_alphabet.unmethylate("CM") == "CC");
    REQUIRE( dcm_alphabet.unmethylate("GM") == "GM");
    REQUIRE( dcm_alphabet.unmethylate("MC") == "MC");

    // disambiguate
    REQUIRE( dcm_alphabet.disambiguate("") == "");
    REQUIRE( dcm_alphabet.disambiguate("M") == "M");
    REQUIRE( dcm_alphabet.disambiguate("CM") == "CM");
    REQUIRE( dcm_alphabet.disambiguate("GM") == "GA");
    REQUIRE( dcm_alphabet.disambiguate("MA") == "MA");
    REQUIRE( dcm_alphabet.disambiguate("MT") == "MT");
    REQUIRE( dcm_alphabet.disambiguate("MC") == "AC");

    // reverse complement
    REQUIRE( dcm_alphabet.reverse_complement("") == "");
    REQUIRE( dcm_alphabet.reverse_complement("M") == "G");
    REQUIRE( dcm_alphabet.reverse_complement("MT") == "AG");
    REQUIRE( dcm_alphabet.reverse_complement("MTG") == "MAG");
    REQUIRE( dcm_alphabet.reverse_complement("MTGG") == "CMAG");
    
    REQUIRE( dcm_alphabet.reverse_complement("MA") == "TG");
    REQUIRE( dcm_alphabet.reverse_complement("MAG") == "MTG");
    REQUIRE( dcm_alphabet.reverse_complement("MAGG") == "CMTG");
    
    REQUIRE( dcm_alphabet.reverse_complement("CM") == "GG");
    REQUIRE( dcm_alphabet.reverse_complement("CCAGG") == "CCTGG");
    REQUIRE( dcm_alphabet.reverse_complement("CCTGG") == "CCAGG");
    REQUIRE( dcm_alphabet.reverse_complement("CMAGG") == "CMTGG");
    REQUIRE( dcm_alphabet.reverse_complement("CMTGG") == "CMAGG");
    
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

    // suffix functions
    REQUIRE( ends_with("abcd", "cd") );
    REQUIRE( ! ends_with("abcd", "bc") );
    REQUIRE( ! ends_with("abcd", "e") );
    REQUIRE( ends_with("abcd", "d") );
    REQUIRE( ends_with("abcd", "") );
}

TEST_CASE( "math", "[math]") {
    GaussianParameters params;
    params.mean = 4;
    params.stdv = 2;
    params.log_stdv = log(params.stdv);

    REQUIRE( normal_pdf(2.25, params) == Approx(0.1360275) );
    REQUIRE( log_normal_pdf(2.25, params) == Approx(log(normal_pdf(2.25, params))) );
}

TEST_CASE( "scalings", "[scalings]") {

    SquiggleRead test_read;
    size_t strand = 0;
    test_read.base_model[0] = PoreModelSet::get_model("r9.4_450bps", "nucleotide", "template", 6);
    const PoreModel* pore_model = test_read.base_model[0];
    test_read.scalings[strand].set4(10.0f, 1.2, 0.5, 1.3);
    
    assert(pore_model != NULL);

    // Generate events from the pore model
    const SquiggleScalings& scalings = test_read.scalings[strand];
    size_t n_events = 100;
    double duration = 10.0 / 4000.0f;
    size_t rank = 100;
    std::default_random_engine generator;

    std::vector<float> lp_truth;
    std::vector<float> z_truth;

    for(size_t i = 0; i < n_events; ++i) {

        PoreModelStateParams params = pore_model->states[rank];
        SquiggleEvent event;
        event.stdv = 1.0f; //unused
        event.start_time = (i * duration);
        event.duration = duration;
        
        // Generate level from N(a + bu + ct, (do)^2)
        float g_mean = scalings.shift + scalings.scale * params.level_mean + event.start_time * scalings.drift;
        float g_stdv = scalings.var * params.level_stdv;

        std::normal_distribution<float> distribution(g_mean, g_stdv);
        event.mean = distribution(generator);
        test_read.events[0].push_back(event);

        GaussianParameters gp = { g_mean, g_stdv };
        lp_truth.push_back(log_normal_pdf(event.mean, gp));
        z_truth.push_back((event.mean - gp.mean) / gp.stdv);
    }

    // Calculate log-probability of observing an event by scaling the event
    for(size_t i = 0; i < n_events; ++i) {
        float lp = log_probability_match_r9(test_read, *pore_model, rank, i, strand);
        float z = z_score(test_read, *pore_model, rank, i, strand);
        REQUIRE( lp == Approx(lp_truth[i]) );
        REQUIRE( z == Approx(z_truth[i]) );
    }
}

size_t factorial(size_t n)
{
    if(n == 0 || n == 1) {
        return 1;
    }

    size_t out = 1;
    for(size_t i = 2; i <= n; ++i) {
        out *= i;
    }
    return out;
}

void test_combinations(size_t n, size_t k, CombinationOption option, std::vector<std::string> expected)
{
    Combinations c(n, k, option);
    size_t idx = 0;
    while(!c.done()) {
        std::string cs = c.get_as_string();
        REQUIRE(cs == expected[idx]);
        idx++;
        c.next();
    }
    REQUIRE(idx == expected.size());
}

TEST_CASE( "combinations", "[combinations]") {
    test_combinations(1, 1, CO_WITHOUT_REPLACEMENT, {"0"});
    test_combinations(2, 1, CO_WITHOUT_REPLACEMENT, { "0", "1" });
    test_combinations(2, 2, CO_WITHOUT_REPLACEMENT, { "0,1" });
    test_combinations(3, 2, CO_WITHOUT_REPLACEMENT, { "0,1", "0,2", "1,2" });
    test_combinations(4, 4, CO_WITHOUT_REPLACEMENT, { "0,1,2,3", });
    REQUIRE(factorial(4) == 24);

    // Count we see the right number of combinations
    size_t n = 10;
    size_t k = 4;
    size_t n_comb = factorial(n) / (factorial(k) * factorial(n - k));
    size_t count = 0;
    Combinations c(n, k);
    while(!c.done()) {
        count++;
        c.next();
    }
    REQUIRE(count == n_comb);

    // With replacement
    test_combinations(1, 1, CO_WITH_REPLACEMENT, {"0"});
    test_combinations(2, 1, CO_WITH_REPLACEMENT, { "0", "1" });
    test_combinations(2, 2, CO_WITH_REPLACEMENT, { "0,0", "0,1", "1,1" });
    test_combinations(3, 2, CO_WITH_REPLACEMENT, { "0,0", "0,1", "0,2", "1,1", "1,2", "2,2"});
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

    // load the FAST5
    std::string read_id = "01234567-0123-0123-0123-0123456789ab:2D_000:2d";
    std::string fast5_path = "test/data/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch8_file30_strand.fast5";

    ReadDB read_db;
    read_db.add_signal_path(read_id, fast5_path);

    SquiggleRead sr(read_id, read_db);

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
    input[0].pore_model = sr.get_model(0, "nucleotide");

    // complement strand
    input[1].read = &sr;
    input[1].event_start_idx = 6788;
    input[1].event_stop_idx = 6697;
    input[1].event_stride = -1;
    input[1].rc = true;
    input[1].strand = 1;
    input[1].pore_model = sr.get_model(1, "nucleotide");
    
    // expected output
    std::string expected_alignment[2];
    expected_alignment[0] = 
        "MMMMMEMKMKMMMMMMMKMMMKMMMKMMMMMMMMMKKMMEEEMMMMMMKMMMM" 
        "MMMKMMMMMKMKMKMEMKKMKMKKMMMMMMEMMMMKMKMEEMMMMKMEEEEEM";

    expected_alignment[1] = 
        "MMKMMMKMEEMMKMKMKMEMMMKMMMKMEMMMKMMMKMMMMMMMMMKKMEMMM"
        "EMMMMMMMMMKMKKMMMMMMMEMMMMMKMMMMMKMEMMMMMKMMMMMEEEEEEEEM";

    double expected_viterbi_last_state[2] = { -237.7808380127, -267.9027709961 };
    double expected_forward[2] = { -216.053604126, -254.5881347656 };

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

std::vector< StateTrainingData >
generate_training_data(const ParamMixture& mixture, size_t n_data,
                       const std::array< float, 2 >& scaled_read_var_rg = { .5f, 1.5f },
                       const std::array< float, 2 >& read_scale_sd_rg = { .5f, 1.5f },
                       const std::array< float, 2 >& read_var_sd_rg = { .5f, 1.5f })
{
    // check parameter sizes
    size_t n_components = mixture.log_weights.size();
    assert(mixture.params.size() == n_components);
    assert(scaled_read_var_rg[0] < scaled_read_var_rg[1]);
    assert(read_scale_sd_rg[0] < read_scale_sd_rg[1]);
    assert(read_var_sd_rg[0] < read_var_sd_rg[1]);
    // set random seed
    //seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    // catch takes care of managing the random seed
    std::mt19937 rg(std::rand());
    typedef std::discrete_distribution< size_t > discrete_dist;
    typedef std::uniform_real_distribution< float > uniform_dist;
    typedef std::normal_distribution< float > normal_dist;
    typedef inverse_gaussian_distribution< float > inverse_gaussian_dist;
    // generate data
    std::vector< float > weights(n_components);
    for (size_t j = 0; j < n_components; ++j)
    {
        weights[j] = std::exp(mixture.log_weights[j]);
        LOG("gen_data", debug) << "weights " << j << " "
                               << std::fixed << std::setprecision(2) << weights[j] << " ("
                               << mixture.log_weights[j] << ")" << std::endl;
    }
    std::vector< float > level_mean_sum(n_components, 0.0);
    std::vector< float > sd_mean_sum(n_components, 0.0);
    std::vector< size_t > population_size(n_components, 0);
    std::vector< StateTrainingData > data(n_data);
    for (size_t i = 0; i < n_data; ++i)
    {
        // draw population
        size_t j = discrete_dist(weights.begin(), weights.end())(rg);
        ++population_size[j];
        assert(j < n_components);
        // draw scaled_read_var
        data[i].scaled_read_var = uniform_dist(scaled_read_var_rg[0], scaled_read_var_rg[1])(rg);
        data[i].log_scaled_read_var = std::log(data[i].scaled_read_var);
        // draw read_scale_sd
        data[i].read_scale_sd = uniform_dist(read_scale_sd_rg[0], read_scale_sd_rg[1])(rg);
        data[i].log_read_scale_sd = std::log(data[i].read_scale_sd);
        // draw read_var_sd
        data[i].read_var_sd = uniform_dist(read_var_sd_rg[0], read_var_sd_rg[1])(rg);
        data[i].log_read_var_sd = std::log(data[i].read_var_sd);
        // scale the state
        auto scaled_params = mixture.params[j];
        scaled_params.level_stdv *= data[i].scaled_read_var;
        scaled_params.level_log_stdv += data[i].log_scaled_read_var;
        scaled_params.sd_lambda *= data[i].read_var_sd / data[i].read_scale_sd;
        scaled_params.sd_log_lambda += data[i].log_read_var_sd - data[i].log_read_scale_sd;
        // draw level_mean & level_stdv
        data[i].level_mean = normal_dist(scaled_params.level_mean, scaled_params.level_stdv)(rg);
        data[i].log_level_mean = std::log(data[i].level_mean);
        data[i].level_stdv = inverse_gaussian_dist(scaled_params.sd_mean, scaled_params.sd_lambda)(rg);
        data[i].log_level_stdv = std::log(data[i].level_stdv);
        level_mean_sum[j] += data[i].level_mean;
        sd_mean_sum[j] += data[i].level_stdv;
        LOG("gen_data", debug1)
            << "data " << i << " " << j << " "
            << data[i].level_mean << " "
            << data[i].level_stdv << " "
            << data[i].scaled_read_var << " "
            << data[i].read_scale_sd << " "
            << data[i].read_var_sd << std::endl;
    }
    for (size_t j = 0; j < n_components; ++j)
    {
        LOG("gen_data", debug)
            << "population " << j << " "
            << std::fixed << std::setprecision(3) << level_mean_sum[j] / population_size[j] << " "
            << sd_mean_sum[j] / population_size[j] << std::endl;
    }
    return data;
}

TEST_CASE("training", "[training]")
{
    const unsigned n_data = 1000;
    const float um_rate = .1;
    PoreModelStateParams um_params;
    um_params.level_mean = 65.0;
    um_params.level_stdv = 1.0;
    um_params.sd_mean = 0.8;
    um_params.sd_lambda = 7.0;
    um_params.update_sd_stdv();
    um_params.update_logs();
    //Logger::set_default_level(level_wrapper::debug);

    // first, we test gaussian training only
    SECTION("gaussian")
    {
        float delta_um_rate = .05;
        float delta_level_mean = 10.0;
        ParamMixture gen_mixture;
        gen_mixture.log_weights.push_back(std::log(um_rate + delta_um_rate));
        gen_mixture.log_weights.push_back(std::log(1 - (um_rate + delta_um_rate)));
        gen_mixture.params.push_back(um_params);
        gen_mixture.params.push_back(um_params);
        gen_mixture.params[1].level_mean += delta_level_mean;
        auto data = generate_training_data(gen_mixture, n_data);
        ParamMixture in_mixture;
        in_mixture.log_weights.push_back(std::log(um_rate));
        in_mixture.log_weights.push_back(std::log(1 - um_rate));
        in_mixture.params.push_back(um_params);
        in_mixture.params.push_back(um_params);
        // encourage the second component to capture points not well fit by the first
        in_mixture.params[1].level_stdv += 1.0;
        in_mixture.params[1].update_logs();
        auto out_mixture = train_gaussian_mixture(data, in_mixture);
        CHECK( std::exp(out_mixture.log_weights[0]) == Approx( um_rate + delta_um_rate ).epsilon(.05) );
        CHECK( out_mixture.params[0].level_mean == Approx( um_params.level_mean ).epsilon(.05) );
        CHECK( out_mixture.params[1].level_mean == Approx( um_params.level_mean + delta_level_mean ).epsilon(.05) );
    }
}
