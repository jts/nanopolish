//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_polya_estimator.cpp -- estimate the length
// of poly-A tails for each read, this time at the sample level.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <iterator>
#include <limits> // for -INFTY
#include "htslib/faidx.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_read_db.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_polya_estimator.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_emissions.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

using namespace std::placeholders;

//
// Getopt
//
#define SUBPROGRAM "polya"

static const char *POLYA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *POLYA_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Estimate the length of the poly-A tail on direct RNA reads\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:w:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "reads",            required_argument, NULL, 'r' },
    { "bam",              required_argument, NULL, 'b' },
    { "genome",           required_argument, NULL, 'g' },
    { "window",           required_argument, NULL, 'w' },
    { "threads",          required_argument, NULL, 't' },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_polya_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << POLYA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << POLYA_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(argc - optind > 0) {
        opt::region = argv[optind++];
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::num_threads <= 0) {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
        die = true;
    }

    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << POLYA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

// Basic HMM struct with fixed transition/emission parameters;
// all of the below is relative to a **scaled & shifted** set of events.
struct PolyAHMM {
    // State labels: [S]TART => 0 ; [L]EADER => 1 ; [A]DAPTER => 2; [P]OLYA => 3; [C]LIFF => 4; [T]RANSCRIPT => 5
    // N.B.: `state transitions` is used to compute log probabilities, as viterbi decoding is done in log-space.
    // state transition probabilities (S->L->A->[P<->C]->T):
    float state_transitions[6][6] = { {0.10f, 0.90f, 0.00f, 0.00f, 0.00f, 0.00f},     // S -> S (10%), S -> L (90%)
                                      {0.00f, 0.90f, 0.10f, 0.00f, 0.00f, 0.00f},     // L -> A (10%), L -> L (90%)
                                      {0.00f, 0.00f, 0.95f, 0.05f, 0.00f, 0.00f},     // A -> P (05%), A -> A (95%)
                                      {0.00f, 0.00f, 0.00f, 0.89f, 0.01f, 0.10f},     // P -> P (89%), P -> C (01%), P -> T (10%)
                                      {0.00f, 0.00f, 0.00f, 0.99f, 0.01f, 0.00f},     // C -> P (99%), C -> C (01%)
                                      {0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f} };   // T -> T (100%)
    float start_probs[6] = { 0.50f, 0.50f, 0.00f, 0.00f, 0.00f, 0.00f }; // 50/50 chance of starting on L or S

    // emission parameters, from empirical MLE on manually-flagged reads:
    // START, LEADER, and POLYA have Gaussian emissions;
    // ADAPTER, TRANSCRIPT have Gaussian mixture emissions;
    // CLIFF has Uniform emissions.
    GaussianParameters s_emission = {70.2737f, 3.7743f};
    GaussianParameters l_emission = {110.973f, 5.237f};
    GaussianParameters a0_emission = {79.347f, 8.3702f};
    GaussianParameters a1_emission = {63.3126f, 2.7464f};
    float a0_coeff = 0.874f;
    float a1_coeff = 0.126f;
    GaussianParameters p_emission = {108.883f, 3.257f};
    float c_begin = 70.0f;
    float c_end = 140.0f;
    float c_log_prob = -4.2485f; // natural log of [1/(140-70)]
    GaussianParameters t0_emission = {79.679f, 6.966f};
    GaussianParameters t1_emission = {105.784f, 16.022f};
    float t0_coeff = 0.346f;
    float t1_coeff = 0.654f;

    // Compute log-probabilities in the constructor:
    float log_state_transitions[6][6];
    float log_start_probs[6];
    PolyAHMM() {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                if (state_transitions[i][j] > 0.00f) {
                    log_state_transitions[i][j] = std::log(state_transitions[i][j]);
                } else {
                    log_state_transitions[i][j] = -INFINITY;
                }
            }

            if (start_probs[i] > 0.00f) {
                log_start_probs[i] = std::log(start_probs[i]);
            } else {
                log_start_probs[i] = -INFINITY;
            }
        } 
    }
};

// Get the log-probability of seeing `x` given we're in state `state` of an HMM
// N.B.: we scale the emission parameters (events are **not** scaled).
float emit_log_proba(const float x, const PolyAHMM hmm, const uint8_t state,
                     const float scale, const float shift, const float var)
{
    // sometimes samples can exceed reasonable bounds due to mechanical issues;
    // in that case, we should clamp it to 100:
    float xx;
    if (x > 200.0f || x < 40.0f) {
        xx = 100.0f;
    } else {
        xx = x;
    }

    // compute on a case-by-case basis to handle heterogeneous probability distributions
    float log_probs;
    if (state == 0) {
        // START state:
        GaussianParameters emission = hmm.s_emission;
        emission.mean = shift + scale*emission.mean;
        emission.stdv = var * emission.stdv;
        emission.log_stdv = std::log(emission.stdv);
        log_probs = log_normal_pdf(xx, emission);
    }
    if (state == 1) {
        // LEADER state:
        GaussianParameters emission = hmm.l_emission;
        emission.mean = shift + scale*emission.mean;
        emission.stdv = var * emission.stdv;
        emission.log_stdv = std::log(emission.stdv);
        log_probs = log_normal_pdf(xx, emission);
    }
    if (state == 2) {
        // ADAPTER state: compute log of gaussian mixture probability
        GaussianParameters emission0 = hmm.a0_emission;
        emission0.mean = shift + scale*emission0.mean;
        emission0.stdv = var * emission0.stdv;
        emission0.log_stdv = std::log(emission0.stdv);
        GaussianParameters emission1 = hmm.a1_emission;
        emission1.mean = shift + scale*emission1.mean;
        emission1.stdv = var * emission1.stdv;
        emission1.log_stdv = std::log(emission1.stdv);
        float a0_coeff = hmm.a0_coeff;
        float a1_coeff = hmm.a1_coeff;
        float mixture_proba = (a0_coeff * normal_pdf(xx, emission0)) + (a1_coeff * normal_pdf(xx, emission1));
        log_probs = std::log(mixture_proba);
    }
    if (state == 3) {
        // POLYA state:
        GaussianParameters emission = hmm.p_emission;
        emission.mean = shift + scale*emission.mean;
        emission.stdv = var * emission.stdv;
        emission.log_stdv = std::log(emission.stdv);
        log_probs = log_normal_pdf(xx, emission);
    }
    if (state == 4) {
        // CLIFF state: middle-out uniform distribution
        if ((xx > hmm.c_begin) && (xx <  hmm.c_end)) {
            log_probs = hmm.c_log_prob;
        } else {
            log_probs = -INFINITY;
        }
    }
    if (state == 5) {
        // TRANSCRIPT state: compute log of gaussian mixture probability
        GaussianParameters emission0 = hmm.t0_emission;
        emission0.mean = shift + scale*emission0.mean;
        emission0.stdv = var * emission0.stdv;
        emission0.log_stdv = std::log(emission0.stdv);
        GaussianParameters emission1 = hmm.t1_emission;
        emission1.mean = shift + scale*emission1.mean;
        emission1.stdv = var * emission1.stdv;
        emission1.log_stdv = std::log(emission1.stdv);
        float t0_coeff = hmm.t0_coeff;
        float t1_coeff = hmm.t1_coeff;
        float mixture_proba = (t0_coeff * normal_pdf(xx, emission0)) + (t1_coeff * normal_pdf(xx, emission1));
        log_probs = std::log(mixture_proba);
    }
    return log_probs;
}

// Run viterbi algorithm given sample sequence and an HMM;
// returns a struct ViterbiOutputs composed of viterbi probs
// and a vector of integers from {0,1,2,3,4,5} == {S,L,A,P,C,T}.
// N.B.: this algorithm takes place in log-space for numerical stability; when we
// write `scores`, this refers to log(P).
struct ViterbiOutputs {
    std::vector<float> scores;
    std::vector<uint8_t> labels;
};

// Viterbi decoding, in the 3'->5' direction:
ViterbiOutputs polya_viterbi(const SquiggleRead& sr, const PolyAHMM hmm)
{
    // get scale/shift/var values, number of events:
    float scale = static_cast<float>(sr.scalings[0].scale);
    float shift = static_cast<float>(sr.scalings[0].shift);
    float var = static_cast<float>(sr.scalings[0].var);
    size_t num_samples = sr.samples.size();

    // create/initialize viterbi scores and backpointers:
    std::vector<float> init_scores(6, -std::numeric_limits<float>::infinity()); // log(0.0) == -INFTY
    std::vector<uint8_t> init_bptrs(6, 6); // 6 == INIT symbol
    std::vector< std::vector<float> > viterbi_scores(num_samples, init_scores);
    std::vector< std::vector<uint8_t> > viterbi_bptrs(num_samples, init_bptrs);

    // forward viterbi pass; fill up backpointers:
    // weight initially distributed between S and L:
    viterbi_scores[0][0] = hmm.log_start_probs[0] + emit_log_proba(sr.samples[num_samples-1], hmm, 0, scale, shift, var);
    viterbi_scores[0][1] = hmm.log_start_probs[1] + emit_log_proba(sr.samples[num_samples-1], hmm, 1, scale, shift, var);
    for (size_t i = 1; i < num_samples; ++i) {
        // `t` moves from 3'->5' on the vector of samples, in opposite direction of `i`:
        size_t t = (num_samples-1-i);
        // get individual incoming state scores:
        float s_to_s = viterbi_scores.at(i-1)[0] + hmm.log_state_transitions[0][0];
        float s_to_l = viterbi_scores.at(i-1)[0] + hmm.log_state_transitions[0][1];
        float l_to_l = viterbi_scores.at(i-1)[1] + hmm.log_state_transitions[1][1];
        float l_to_a = viterbi_scores.at(i-1)[1] + hmm.log_state_transitions[1][2];
        float a_to_a = viterbi_scores.at(i-1)[2] + hmm.log_state_transitions[2][2];
        float a_to_p = viterbi_scores.at(i-1)[2] + hmm.log_state_transitions[2][3];
        float p_to_p = viterbi_scores.at(i-1)[3] + hmm.log_state_transitions[3][3];
        float p_to_c = viterbi_scores.at(i-1)[3] + hmm.log_state_transitions[3][4];
        float p_to_t = viterbi_scores.at(i-1)[3] + hmm.log_state_transitions[3][5];
        float c_to_c = viterbi_scores.at(i-1)[4] + hmm.log_state_transitions[4][4];
        float c_to_p = viterbi_scores.at(i-1)[4] + hmm.log_state_transitions[4][3];
        float t_to_t = viterbi_scores.at(i-1)[5] + hmm.log_state_transitions[5][5];

        // update the viterbi scores for each state at this timestep:
        viterbi_scores.at(i)[0] = s_to_s + emit_log_proba(sr.samples[t], hmm, 0, scale, shift, var);
        viterbi_scores.at(i)[1] = std::max(l_to_l, s_to_l) + emit_log_proba(sr.samples[t], hmm, 1, scale, shift, var);
        viterbi_scores.at(i)[2] = std::max(a_to_a, l_to_a) + emit_log_proba(sr.samples[t], hmm, 2, scale, shift, var);
        viterbi_scores.at(i)[3] = std::max(p_to_p, std::max(a_to_p, c_to_p)) + emit_log_proba(sr.samples[t], hmm, 3, scale, shift, var);
        viterbi_scores.at(i)[4] = std::max(c_to_c, p_to_c) + emit_log_proba(sr.samples[t], hmm, 4, scale, shift, var);
        viterbi_scores.at(i)[5] = std::max(p_to_t, t_to_t) + emit_log_proba(sr.samples[t], hmm, 5, scale, shift, var);

        // backpointers:
        // START: S can only come from S
        viterbi_bptrs.at(i)[0] = 0;
        // LEADER: L->L or S->L
        if (s_to_l < l_to_l) {
            viterbi_bptrs.at(i)[1] = 1;
        } else {
            viterbi_bptrs.at(i)[1] = 0;
        }
        // ADAPTER:
        if (l_to_a < a_to_a) {
            viterbi_bptrs.at(i)[2] = 2;
        } else {
            viterbi_bptrs.at(i)[2] = 1;
        }
        // POLYA:
        if ((a_to_p < p_to_p) && (c_to_p < p_to_p)) {
            viterbi_bptrs.at(i)[3] = 3;
        } else if ((p_to_p < a_to_p) && (c_to_p < a_to_p)) {
            viterbi_bptrs.at(i)[3] = 2;
        } else {
            viterbi_bptrs.at(i)[3] = 4;
        }
        // CLIFF:
        if (p_to_c < c_to_c) {
            viterbi_bptrs.at(i)[4] = 4;
        } else {
            viterbi_bptrs.at(i)[4] = 3;
        }
        // TRANSCRIPT:
        if (p_to_t < t_to_t) {
            viterbi_bptrs.at(i)[5] = 5;
        } else {
            viterbi_bptrs.at(i)[5] = 3;
        }

    }

    // backwards viterbi pass:
    // allocate `regions` vector of same dimensions as sample sequence;
    // clamp final state to 'T' ~ transcript:
    std::vector<uint8_t> regions(num_samples, 0);
    std::vector<float> scores(num_samples, 0);
    regions[num_samples-1] = 5;
    scores[num_samples-1] = viterbi_scores.at(num_samples-1)[5];
    // loop backwards and keep appending best states:
    for (size_t j=(num_samples-2); j > 0; --j) {
        regions[j] = viterbi_bptrs.at(j)[regions.at(j+1)];
        scores[j] = viterbi_scores.at(j)[regions.at(j+1)];
    }

    // put into struct and return:
    ViterbiOutputs output_vectors = { scores, regions };
    return output_vectors;
}

// HMM construct/viterbi wrapper for segmentation of a squiggle read
ViterbiOutputs segment_read_into_regions(const SquiggleRead& sr)
{
    PolyAHMM polya_hmm;
    ViterbiOutputs output_vectors = polya_viterbi(sr, polya_hmm);
    return output_vectors;
}

// helper fn that gets the final indices of each of the first three regions
struct RegionIxs {
    size_t start;   // final index of S; might not exist if skipped over
    size_t leader;  // final index of L, as indicated by 3'->5' viterbi
    size_t adapter; // final index of A, as indicated by 3'->5' viterbi
    size_t polya;   // final index of P/C, as indicated by 3'->5' viterbi
    size_t cliffs;  // number of observed 'CLIFF' samples
};
RegionIxs get_region_indices(const std::vector<uint8_t> region_labels)
{
    // initial values for indices should preserve expected order:
    RegionIxs ixs = { 0, 1, 2, 3, 0 };

    // loop through sequence and collect values:
    for (std::vector<uint8_t>::size_type i = 0; i < region_labels.size(); ++i) {
        // call end of START:
        if (region_labels[i] == 0 && region_labels[i+1] == 1) {
            ixs.start = static_cast<size_t>(i);
        }
        // call end of leader:
        if (region_labels[i] == 1 && region_labels[i+1] == 2) {
            ixs.leader = static_cast<size_t>(i);
        }
        // call end of adapter:
        if (region_labels[i] == 2 && region_labels[i+1] == 3) {
            ixs.adapter = static_cast<size_t>(i);
        }
        // call end of polya:
        if (region_labels[i] == 3 && region_labels[i+1] == 5) {
            ixs.polya = static_cast<size_t>(i);
        }
        // increment cliff counter:
        if (region_labels[i] == 4) {
            ixs.cliffs++;
        }
    }

    // set sensible (easy to QC-filter) default values if not all four detected:
    // S-end is always detected (min value == 0)
    if (ixs.leader == 1 || ixs.adapter == 2 || ixs.polya == 3) {
        ixs.leader = region_labels.size() - 3;
        ixs.adapter = region_labels.size() - 2;
        ixs.polya = region_labels.size() - 1;
    }
    
    return ixs;
}

// Write Poly(A) region segmentation data to TSV
void estimate_polya_for_single_read_hmm(const ReadDB& read_db,
                                        const faidx_t* fai,
                                        FILE* out_fp,
                                        const bam_hdr_t* hdr,
                                        const bam1_t* record,
                                        size_t read_idx,
                                        int region_start,
                                        int region_end)
{
    //----- load a squiggle read
    std::string read_name = bam_get_qname(record);

    //----- get length of suffix of the read that was softclipped:
    size_t n_cigar = record->core.n_cigar;
    uint32_t prefix_cigar = bam_get_cigar(record)[0];
    uint32_t suffix_cigar = bam_get_cigar(record)[n_cigar - 1];

    uint32_t prefix_clip = bam_cigar_oplen(prefix_cigar);
    uint32_t suffix_clip = bam_cigar_oplen(suffix_cigar);

    SquiggleRead sr(read_name, read_db, SRF_LOAD_RAW_SAMPLES);

    //----- print clipping data if `verbose > 2` set:
    if (opt::verbose > 2) {
        fprintf(stderr, "[polya] read: %s length: %zu prefix clip: %zu suffix clip %zu\n",
                read_name.c_str(), sr.read_sequence.length(), prefix_clip, suffix_clip);
    }
    std::string sequenced_transcript = sr.read_sequence;

    //----- QC: skip this read if long skip at end or if most of transcript wasnt aligned
    if (suffix_clip > 200 || (double)(prefix_clip + suffix_clip) / sequenced_transcript.length() > 0.2) {
        return;
    }

    //----- QC: skip if no events:
    if (sr.events[0].empty()) {
        return;
    }

    //----- perform HMM-based regional segmentation and get regions:
    ViterbiOutputs v_out = segment_read_into_regions(sr);
    RegionIxs region_indices = get_region_indices(v_out.labels);

    //----- compute output values:
    // start and end times (sample indices) of the poly(A) tail, in original 3'->5' time-direction:
    // (n.b.: everything in 5'->3' order due to inversion in SquiggleRead constructor, but our
    // `region_indices` struct has everything in 3'->5' order)
    double num_samples = sr.samples.size();
    double polya_sample_start = region_indices.adapter + 1;
    double polya_sample_end = region_indices.polya;
    double adapter_sample_start = region_indices.leader;
    double leader_sample_start = region_indices.start;
    // calculate duration of poly(A) region (in seconds)
    double duration = (region_indices.polya - (region_indices.adapter + 1)) / sr.sample_rate;
    // calculate read duration (length of transcript, in seconds) and read rate:
    double read_duration = (num_samples - polya_sample_end) / sr.sample_rate;
    double read_rate = (sequenced_transcript.length() - suffix_clip) / read_duration;
    // length of the poly(A) tail, in nucleotides:
    double polya_length = duration * read_rate;
    // number of cliffs observed:
    size_t num_cliffs = region_indices.cliffs;
    // estimated adapter length, in nucleotides:
    double adapter_duration = (region_indices.adapter - (region_indices.leader - 1)) / sr.sample_rate;
    double adapter_length = adapter_duration * read_rate;

    //----- QC: check for NOREGION and faulty ADAPTER:
    // `adapter_qc_tol` is the upper-tolerance for number of estimated adapter nucleotides; discovered empirically
    double num_leader_samples = adapter_sample_start;
    double num_adapter_samples = polya_sample_start - adapter_sample_start;
    double num_polya_samples = polya_sample_end - polya_sample_start;
    double adapter_qc_tol = 300.0f;
    std::string qc_tag;
    if (num_leader_samples < 200.0 || num_adapter_samples < 200.0 || num_polya_samples < 200.0) {
        qc_tag = "NOREGION";
    } else if (adapter_length > adapter_qc_tol) {
        qc_tag = "ADAPTER";
    } else {
        qc_tag = "PASS";
    }

    //----- print to TSV:
    #pragma omp critical
    fprintf(out_fp, "polya-annotation\t%s\t%zu\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%zu\t%s\n",
            read_name.c_str(), record->core.pos, polya_sample_start, polya_sample_end,
            read_rate, polya_length, num_cliffs, qc_tag.c_str());

    //----- if `verbose == 1`, print out the full read segmentation:
    #pragma omp critical
    if (opt::verbose == 1) {
        fprintf(out_fp, "polya-segmentation\t%s\t%zu\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%.2lf\n",
                read_name.c_str(), record->core.pos,
                leader_sample_start, adapter_sample_start, polya_sample_start, polya_sample_end,
                read_rate, polya_length, adapter_length);
    }
    //----- if `verbose >= 2`, print the samples (picoAmps) of the read,
    // up to the first 1000 samples of transcript region:
    #pragma omp critical
    if (opt::verbose >= 2) {
        // copy 5'->3'-oriented samples from squiggleread and reverse back to 3'->5':
        std::vector<float> samples(sr.samples);
        std::reverse(samples.begin(), samples.end());
        const PolyAHMM hmm;
        std::string ref_name(hdr->target_name[record->core.tid]);
        for (size_t i = 0; i < std::min(static_cast<size_t>(polya_sample_end)+1000, samples.size()); ++i) {
            std::string tag;
            if (i < leader_sample_start) {
                tag = "START";
            } else if (i < adapter_sample_start) {
                tag = "LEADER";
            } else if (i < polya_sample_start) {
                tag =  "ADAPTER";
            } else if (i < polya_sample_end) {
                tag = "POLYA";
            } else {
                tag = "TRANSCRIPT";
            }

            double s = samples[i];
            double scaled_s = (s - sr.scalings[0].shift) / sr.scalings[0].scale;
            double s_proba_0 = emit_log_proba(s, hmm, 0, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            double s_proba_1 = emit_log_proba(s, hmm, 1, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            double s_proba_2 = emit_log_proba(s, hmm, 2, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            double s_proba_3 = emit_log_proba(s, hmm, 3, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            double s_proba_4 = emit_log_proba(s, hmm, 4, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            double s_proba_5 = emit_log_proba(s, hmm, 5, sr.scalings[0].scale,
                                              sr.scalings[0].shift, sr.scalings[0].var);
            fprintf(out_fp, "polya-samples\t%s\t%s\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                    read_name.substr(0,6).c_str(), ref_name.c_str(), i, s, scaled_s,
                    s_proba_0, s_proba_1, s_proba_2, s_proba_3, s_proba_4, s_proba_5,
                    tag.c_str());
        }
    }

}

//
int polya_main(int argc, char** argv)
{
    parse_polya_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(estimate_polya_for_single_read_hmm, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    // free allocated values:
    fai_destroy(fai);

    return EXIT_SUCCESS;
}
