//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_polya_estimator.cpp -- estimate the length
// of poly-A tails for each read
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
#include <limits> // for -INFTY

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
    // state transition probabilities:
    float state_transitions[4][4] = { {0.99f, 0.01f, 0.00f, 0.00f},
                                      {0.00f, 0.99f, 0.01f, 0.00f},
                                      {0.00f, 0.00f, 0.99f, 0.01f},
                                      {0.00f, 0.00f, 0.00f, 1.00f} };
    // log-prob transitions:
    float log_state_transitions[4][4] = { {-0.0101f, -4.6052f, -INFINITY, -INFINITY},
                                          {-INFINITY, -0.0101f, -4.6052f, -INFINITY},
                                          {-INFINITY, -INFINITY, -0.0101f, -4.6052f},
                                          {-INFINITY, -INFINITY, -INFINITY,  0.0000f} };
    // default emission parameters, from empirical MLE:
    GaussianParameters l_emission = {110.276f, 2.815f};
    GaussianParameters a_emission = {79.898f, 3.0883f};
    GaussianParameters p_emission = {111.352f, 2.085f}; 
    GaussianParameters t_emission = {99.016f, 4.339f};
    // for reference: hand-tuned emission parameters for older kmer model:
    // LEADER {100.0f, 2.5f};
    // ADAPTER {75.0f, 3.0f};
    // POLYA {99.0f, 2.5f};
    // TRANSCRIPT {93.0f, 5.0f};
};

// Get the log-probability of seeing `x` given we're in state `state` of an HMM
// N.B.: we scale the emission parameters (events are **not** scaled).
float emit_log_proba(const float x, const PolyAHMM hmm, const uint8_t state,
                     const float scale, const float shift, const float var)
{
    GaussianParameters emission;
    if (state == 0) {
        emission = hmm.l_emission;
    }
    if (state == 1) {
        emission = hmm.a_emission;
    }
    if (state == 2) {
        emission = hmm.p_emission;
    }
    if (state == 3) {
        emission = hmm.t_emission;
    }
    emission.mean = shift + scale*emission.mean;
    emission.stdv = var * emission.stdv;
    emission.log_stdv = std::log(emission.stdv);
    float log_probs = log_normal_pdf(x, emission);
    return log_probs;
}

// Run viterbi algorithm given event sequence and an HMM;
// returns a struct ViterbiOutputs composed of viterbi probs
// and a vector of integers from {0,1,2,3} == {L,A,P,T}.
// N.B.: this algorithm takes place in log-space for numerical stability; when we
// write `scores`, this refers to log(P).
struct ViterbiOutputs {
    std::vector<float> scores;
    std::vector<uint8_t> labels;
};
ViterbiOutputs polya_viterbi(const SquiggleRead& sr, const PolyAHMM hmm)
{
    // get scale/shift/var values, number of events:
    float scale = static_cast<float>(sr.scalings[0].scale);
    float shift = static_cast<float>(sr.scalings[0].shift);
    float var = static_cast<float>(sr.scalings[0].var);
    size_t num_events = sr.events[0].size();

    // create/initialize viterbi scores and backpointers:
    std::vector<float> init_scores(4, -std::numeric_limits<float>::infinity()); // log(0.0) == -INFTY
    std::vector<uint8_t> init_bptrs(4, 4); // 4 == INIT symbol
    std::vector< std::vector<float> > viterbi_scores(num_events, init_scores);
    std::vector< std::vector<uint8_t> > viterbi_bptrs(num_events, init_bptrs);

    // forward viterbi pass; fill up backpointers:
    // weight initially on L:
    viterbi_scores[0][0] = emit_log_proba(sr.events[0][num_events-1].mean, hmm, 0, scale, shift, var);
    for (size_t i=1; i < num_events; ++i) {
        // `t` moves from 3'->5' on the vector of events, in opposite direction of `i`:
        size_t t = (num_events-1-i);
        // get scores and update timestep:
        float l_to_l = viterbi_scores.at(i-1)[0] + hmm.log_state_transitions[0][0];
        float l_to_a = viterbi_scores.at(i-1)[0] + hmm.log_state_transitions[0][1];
        float a_to_a = viterbi_scores.at(i-1)[1] + hmm.log_state_transitions[1][1];
        float a_to_p = viterbi_scores.at(i-1)[1] + hmm.log_state_transitions[1][2];
        float p_to_p = viterbi_scores.at(i-1)[2] + hmm.log_state_transitions[2][2];
        float p_to_t = viterbi_scores.at(i-1)[2] + hmm.log_state_transitions[2][3];
        float t_to_t = viterbi_scores.at(i-1)[3] + hmm.log_state_transitions[3][3];

        viterbi_scores.at(i)[0] = l_to_l + emit_log_proba(sr.events[0][t].mean, hmm, 0, scale, shift, var);
        viterbi_scores.at(i)[1] = std::max(l_to_a, a_to_a) + emit_log_proba(sr.events[0][t].mean, hmm, 1, scale, shift, var);
        viterbi_scores.at(i)[2] = std::max(a_to_p, p_to_p) + emit_log_proba(sr.events[0][t].mean, hmm, 2, scale, shift, var);
        viterbi_scores.at(i)[3] = std::max(p_to_t, t_to_t) + emit_log_proba(sr.events[0][t].mean, hmm, 3, scale, shift, var);

        // backpointers:
        uint8_t l_bptr = 0; // L can only come from L
        uint8_t a_bptr = (l_to_a < a_to_a) ? 1 : 0; // A->A : L->A
        uint8_t p_bptr = (a_to_p < p_to_p) ? 2 : 1; // P->P : A->P
        uint8_t t_bptr = (p_to_t < t_to_t) ? 3 : 2; // T->T : P->T
        viterbi_bptrs.at(i)[0] = l_bptr;
        viterbi_bptrs.at(i)[1] = a_bptr;
        viterbi_bptrs.at(i)[2] = p_bptr;
        viterbi_bptrs.at(i)[3] = t_bptr;
    }

    // backwards viterbi pass:
    // allocate `regions` vector of same dimensions as event sequence;
    // clamp final state to 'T' ~ transcript:
    std::vector<uint8_t> regions(num_events, 0);
    std::vector<float> scores(num_events, 0);
    regions[num_events-1] = 3;
    scores[num_events-1] = viterbi_scores.at(num_events-1)[3];
    // loop backwards and keep appending best states:
    for (size_t j=(num_events-2); j > 0; --j) {
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
    size_t leader;
    size_t adapter;
    size_t polya;
};
RegionIxs get_region_indices(const std::vector<uint8_t> regions)
{
    RegionIxs ixs = { 0, 0, 0 };

    // loop through sequence and collect values:
    for (std::vector<uint8_t>::size_type i = 0; i < regions.size(); ++i) {
        // call end of leader:
        if (regions[i] == 0 && regions[i+1] == 1) {
            ixs.leader = static_cast<size_t>(i);
        }
        // call end of adapter:
        if (regions[i] == 1 && regions[i+1] == 2) {
            ixs.adapter = static_cast<size_t>(i);
        }
        // call end of polya:
        if (regions[i] == 2 && regions[i+1] == 3) {
            ixs.polya = static_cast<size_t>(i);
        }
    }

    // set sensible default values (1 event past previous region) if not all three detected:
    // L-end is always detected (min value == 0)
    if (ixs.adapter == 0) {
        // A-end undetected if all events after L are all A's; set final 3 events as A,P,T:
        ixs.adapter = regions.size() - 3;
        ixs.polya = regions.size() - 2;
    } else if (ixs.polya == 0) {
        // P-end undetected if all events after A are all P's; set final events to P,T:
        ixs.polya = regions.size() - 2;
    }

    return ixs;
}

// helper fn that checks a segmented read for an incomplete poly(A) region; returns
// `1` if the poly(A) region is incomplete and 0 if it passes.
// an INCOMPLETE poly(A) call is a situation where a spuriously high/low signal lasting
// approx. 100 samples throws off the HMM and prevents it from seeing the rest of the
// region.
uint8_t qc_incomplete_polya(const SquiggleRead& sr, double polya_start, double polya_end)
{
    // qc parameters, found empirically:
    double multiplier = 2.5;
    double polya_stdv_cutoff = 4.0;
    size_t clip_samples = 50;
    size_t skip_samples = 100;
    size_t downstream_block_size = 100;

    // get number of samples; we need this to compensate for the 5'->3' ordering of `sr.samples`
    size_t num_samples = sr.samples.size();

    // start and stop indices for called poly(A) region, taking clipping and reversal into account:
    size_t polya_start_ix = num_samples - (size_t)polya_end + clip_samples;
    size_t polya_stop_ix = num_samples - (size_t)polya_start - clip_samples;

    // compute STDV of inner region of poly(A) samples:
    double polya_mean_acc = 0.0;
    for (size_t ix = polya_start_ix; ix < polya_stop_ix; ++ix) {
	polya_mean_acc += sr.samples[ix];
    }
    double polya_mean = polya_mean_acc / (polya_stop_ix - polya_start_ix - 2*clip_samples);
    double polya_stdv_acc = 0.0;
    for (size_t ix = polya_start_ix; ix < polya_stop_ix; ++ix) {
	polya_stdv_acc += (sr.samples[ix] - polya_mean) * (sr.samples[ix] - polya_mean);
    }
    double polya_stdv = sqrt(polya_stdv_acc / (polya_stop_ix - polya_start_ix - 2*clip_samples));

    // start and stop indices for downstream region, taking 5'->3' orientation into account:
    size_t downstream_start_ix = num_samples - (size_t)polya_end - skip_samples - downstream_block_size;
    size_t downstream_stop_ix = num_samples -(size_t)polya_end - skip_samples;

    // compute STDV of downstream region:
    double downstream_mean_acc = 0.0;
    for (size_t ix = downstream_start_ix; ix < downstream_stop_ix; ++ix) {
	downstream_mean_acc += sr.samples[ix];
    }
    double downstream_mean = downstream_mean_acc / downstream_block_size;
    double downstream_stdv_acc = 0.0;
    for (size_t ix = downstream_start_ix; ix < downstream_stop_ix; ++ix) {
	downstream_stdv_acc += (sr.samples[ix] - downstream_mean) * (sr.samples[ix] - downstream_mean);
    }
    double downstream_stdv = sqrt(downstream_stdv_acc / downstream_block_size);

    // if the downstream region has a STDV that's more stable than expected, we have an incomplete poly(A):
    uint8_t qc_incomplete_flag;
    if (downstream_stdv < multiplier * polya_stdv && polya_stdv < polya_stdv_cutoff) {
        qc_incomplete_flag = 1;
    } else {
        qc_incomplete_flag = 0;
    }

    //std::cout << "POLYA STDV: " << polya_stdv << "| DOWNSTREAM STDV: " << downstream_stdv << std::endl;

    return qc_incomplete_flag;
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
    ViterbiOutputs viterbi_outs = segment_read_into_regions(sr);
    RegionIxs region_indices = get_region_indices(viterbi_outs.labels);

    //----- compute output values:
    // start and end times (sample indices) of the poly(A) tail, in original 3'->5' time-direction:
    // (n.b.: everything in 5'->3' order due to inversion in SquiggleRead constructor, but our
    // `region_indices` struct has everything in 3'->5' order)
    int STRAND = 0;
    size_t num_events = sr.events[0].size();
    double total_num_samples = sr.samples.size();
    double polya_sample_start = sr.events[STRAND][num_events - region_indices.adapter - 1].start_time * sr.sample_rate;
    double polya_sample_end = sr.events[STRAND][num_events - region_indices.polya].start_time * sr.sample_rate;
    double adapter_sample_start = sr.events[STRAND][num_events - region_indices.leader].start_time * sr.sample_rate;
    // calculate duration of poly(A) region (in seconds)
    double duration = sr.events[STRAND][num_events - region_indices.polya].start_time
        - sr.events[STRAND][num_events - region_indices.adapter - 1].start_time;
    // calculate read duration (length of transcript, in seconds) and read rate:
    double read_duration = (total_num_samples - polya_sample_end) / sr.sample_rate;
    double read_rate = (sequenced_transcript.length() - suffix_clip) / read_duration;
    // length of the poly(A) tail, in nucleotides:
    double polya_length = duration * read_rate;

    //----- QC: check for NOREGION and INCOMPLETE poly(A) calls (with NOREGION taking precedence):
    double num_leader_samples = adapter_sample_start;
    double num_adapter_samples = polya_sample_start - adapter_sample_start;
    double num_polya_samples = polya_sample_end - polya_sample_start;
    std::string qc_tag;
    if (num_leader_samples < 200.0 || num_adapter_samples < 200.0 || num_polya_samples < 200.0) {
        qc_tag = "NOREGION";
    } else if (qc_incomplete_polya(sr,polya_sample_start,polya_sample_end) == 1) {
        qc_tag = "INCOMPLETE";
    } else {
        qc_tag = "PASS";
    }

    //----- print to TSV:
    fprintf(out_fp, "polya-annotation\t%s\t%zu\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%s\n",
            read_name.c_str(), record->core.pos, polya_sample_start, polya_sample_end,
            read_rate, polya_length, qc_tag.c_str());

    //----- if `verbose == 1`, print out the full read segmentation:
    if (opt::verbose == 1) {
        fprintf(out_fp, "polya-segmentation\t%s\t%zu\t%.1lf\t%.1lf\t%.1lf\t%.2lf\t%.2lf\n",
                read_name.c_str(), record->core.pos,
                adapter_sample_start, polya_sample_start, polya_sample_end,
                read_rate, polya_length);
    }
    //----- if `verbose >= 2`, print the samples (picoAmps) of the read,
    // up to the first 1000 samples of transcript region:
    if (opt::verbose >= 2) {
        // copy 5'->3'-oriented samples from squiggleread and reverse back to 3'->5':
        std::vector<float> samples(sr.samples);
        std::reverse(samples.begin(), samples.end());
        const PolyAHMM hmm;
	std::string ref_name(hdr->target_name[record->core.tid]);
        for (size_t i = 0; i < std::min(static_cast<size_t>(polya_sample_end)+1000, samples.size()); ++i) {
            std::string tag;
            if (i < adapter_sample_start) {
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
            fprintf(out_fp, "polya-samples\t%s\t%s\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                    read_name.substr(0,6).c_str(), ref_name.c_str(), i, s, scaled_s,
                    s_proba_0, s_proba_1, s_proba_2, s_proba_3,
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
