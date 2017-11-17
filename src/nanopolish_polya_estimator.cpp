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

//==================================================================================================
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
    // default emission parameters, from eyeballing a bunch of plots:
    GaussianParameters l_emission = {100.0f, 2.5f};  // leader;
    GaussianParameters a_emission = {75.0f, 3.0f};   // adapter;
    GaussianParameters p_emission = {99.0f, 2.5f};   // polya;
    GaussianParameters t_emission = {93.0f, 5.0f};  // transcript;
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
    float scale = static_cast<float>(sr.pore_model[0].scale);
    float shift = static_cast<float>(sr.pore_model[0].shift);
    float var = static_cast<float>(sr.pore_model[0].var);
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

    //----- print TSV line:
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


    // print to TSV:
    fprintf(out_fp, "polya-annotation\t%s\t%zu\t%.1lf\t%.1lf\t%.2lf\t%.2lf\n",
            read_name.c_str(), record->core.pos, polya_sample_start, polya_sample_end,
            read_rate, polya_length);

    //----- if `verbose >= 1`, print the samples (picoAmps) of the read,
    // up to the first 1000 samples of transcript region:
    if (opt::verbose >= 1) {
        // copy 5'->3'-oriented samples from squiggleread and reverse back to 3'->5':
        std::vector<float> samples(sr.samples);
        std::reverse(samples.begin(), samples.end());
        const PolyAHMM hmm;
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
            double scaled_s = (s - sr.pore_model[0].shift) / sr.pore_model[0].scale;
            double s_proba_0 = emit_log_proba(s, hmm, 0, sr.pore_model[0].scale,
                                              sr.pore_model[0].shift, sr.pore_model[0].var);
            double s_proba_1 = emit_log_proba(s, hmm, 1, sr.pore_model[0].scale,
                                              sr.pore_model[0].shift, sr.pore_model[0].var);
            double s_proba_2 = emit_log_proba(s, hmm, 2, sr.pore_model[0].scale,
                                              sr.pore_model[0].shift, sr.pore_model[0].var);
            double s_proba_3 = emit_log_proba(s, hmm, 3, sr.pore_model[0].scale,
                                              sr.pore_model[0].shift, sr.pore_model[0].var);
            fprintf(out_fp, "polya-samples\t%s\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                    read_name.substr(0,6).c_str(), i, s, scaled_s,
                    s_proba_0, s_proba_1, s_proba_2, s_proba_3,
                    tag.c_str());
        }
    }

}

//==================================================================================================

double bhattacharyya_coefficient(const PoreModelStateParams& a,
                                 const PoreModelStateParams& b)
{
    double var_a = pow(a.level_stdv, 2.0);
    double var_b = pow(b.level_stdv, 2.0);
    double term_s = ( (var_a / var_b) + (var_b / var_a) + 2);
    double term_m = (pow(a.level_mean - b.level_mean, 2.0) / (var_a + var_b));
    return 0.25 * log( 0.25 * term_s) + 0.25 * term_m;
}

// Partition the input sequence into segments that have similar current levels based on the pore model
// The similarity between current levels for adjacent k-mers is determined by the Bhattacharyya coefficient
// with a tunable mininum distance
std::vector<size_t> partition_sequence(const std::string& sequence,
                                       double min_distance,
                                       const PoreModel& pore_model)
{
    size_t k = pore_model.k;
    size_t num_kmers = sequence.length() - k + 1;
    std::vector<size_t> partitions(num_kmers);
    size_t curr_partition = 0;
    partitions[0] = curr_partition;
    for(size_t i = 1; i < num_kmers; ++i) {
        uint32_t prev_rank = pore_model.pmalphabet->kmer_rank(sequence.c_str() + i - 1, k);
        uint32_t curr_rank = pore_model.pmalphabet->kmer_rank(sequence.c_str() + i, k);
        const PoreModelStateParams& prev_params = pore_model.states[prev_rank];
        const PoreModelStateParams& curr_params = pore_model.states[curr_rank];
        double bc = bhattacharyya_coefficient(curr_params, prev_params);
        /*
           fprintf(stderr, "a: (%.1lf, %.2lf) b: (%.1lf, %.2lf) bc: %.3lf\n", curr_params.level_mean, curr_params.level_stdv,
           prev_params.level_mean, prev_params.level_stdv,
           bc);
         */
        if(bc < min_distance) {
            partitions[i] = curr_partition;
        } else {
            partitions[i] = ++curr_partition;
        }
    }
    return partitions;
}

//
void estimate_polya_for_single_read(const ReadDB& read_db,
                                    const faidx_t* fai,
                                    FILE* out_fp,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);

    // Get the length of the suffix of the read that has been softclipped
    size_t n_cigar = record->core.n_cigar;
    uint32_t prefix_cigar = bam_get_cigar(record)[0];
    uint32_t suffix_cigar = bam_get_cigar(record)[n_cigar - 1];

    uint32_t prefix_clip = bam_cigar_oplen(prefix_cigar);
    uint32_t suffix_clip = bam_cigar_oplen(suffix_cigar);

    // load read
    SquiggleRead sr(read_name, read_db, SRF_LOAD_RAW_SAMPLES);

    if(opt::verbose > 2) {
        fprintf(stderr, "[polya] read: %s length: %zu prefix clip: %zu suffix clip: %zu\n", read_name.c_str(), sr.read_sequence.length(), prefix_clip, suffix_clip);
    }

    std::string sequenced_transcript = sr.read_sequence;

    // QC - skip reads with long skips at the end or where most of the transcript was not aligned
    if(suffix_clip > 200 || (double)(prefix_clip + suffix_clip) / sequenced_transcript.length() > 0.2) {
        return;
    }

    // Align events to transcript
    std::vector<AlignedPair> alignment = adaptive_banded_simple_event_align(sr, sequenced_transcript);

    // Partition the read into blocks with a similar current value
    std::vector<size_t> partitions = partition_sequence(sequenced_transcript, 0.5, sr.pore_model[0]);
    size_t num_partitions = partitions.back() + 1;

    std::vector<double> start_time_by_partition(num_partitions, INFINITY);
    std::vector<double> events_by_partition(num_partitions, 0);
    std::vector<double> sum_z_by_partition(num_partitions, 0);

    size_t k = sr.pore_model[0].k;

    for(size_t i = 0; i < alignment.size(); ++i) {
        size_t event_idx = alignment[i].read_pos;
        size_t kmer_idx = alignment[i].ref_pos;
        size_t partition = partitions[kmer_idx];
        double time = sr.events[0][event_idx].start_time;
        if(time < start_time_by_partition[partition]) {
            start_time_by_partition[partition] = time;
        }

        std::string kmer = sequenced_transcript.substr(kmer_idx, k);
        size_t polya_rank = sr.pore_model[0].pmalphabet->kmer_rank("AAAAAA", k);
        double z = z_score(sr, polya_rank, event_idx, 0);
        sum_z_by_partition[partition] += z;
        events_by_partition[partition] += 1;
        if(opt::verbose > 2) {
            fprintf(stderr, "[polya] %zu %zu %zu %s partition: %zu time: %.0lf z_score: %.2lf\n",
                i, kmer_idx, event_idx, kmer.c_str(), partition, start_time_by_partition[partition] * sr.sample_rate, z);
        }
    }

    // Select the partition that indicates the poly-A tail using the following conditions:
    // 1) it is within 10bp of the end of the transcripts alignment
    // 2) it is not the last partition in the read (which is part of the non-polyA adapter)
    // 3) at least 500 samples in length
    // 4) it has the minimum average absolute z-score across all partitions meeting 1-3
    double min_duration = 500;
    size_t min_distance_to_alignment_end = 10;

    size_t first_valid_base = sequenced_transcript.size() - suffix_clip - min_distance_to_alignment_end;
    size_t first_valid_partition = partitions[first_valid_base];
    size_t best_partition = 0;
    double best_z = INFINITY;
    for(size_t i = first_valid_partition; i < num_partitions - 1; ++i) {
        double duration_in_samples = (start_time_by_partition[i - 1] - start_time_by_partition[i]) * sr.sample_rate;
        double avg_z = sum_z_by_partition[i] / events_by_partition[i];

        if(opt::verbose > 1) {
            fprintf(stderr, "checking partition %zu duration: %.0f z: %.2f\n", i, duration_in_samples, avg_z);
        }
        if(duration_in_samples > min_duration && abs(avg_z) < best_z) {
            best_z = abs(avg_z);
            best_partition = i;
        }
    }

    if(best_partition == 0) {
        fprintf(stderr, "Could not detect poly-A tail for %s\n", read_name.c_str());
        return;
    }

    // write inferred poly-A position in samples
    double sample_start = start_time_by_partition[best_partition] * sr.sample_rate;
    double sample_end = start_time_by_partition[best_partition - 1] * sr.sample_rate;
    double duration = start_time_by_partition[best_partition - 1] - start_time_by_partition[best_partition];

    // get the length of transcript part of the read, in seconds
    double read_duration = start_time_by_partition[0] - start_time_by_partition[best_partition - 1];
    double read_rate = (sequenced_transcript.length() - suffix_clip) / read_duration;

    double polya_length = duration * read_rate;
    fprintf(out_fp, "polya-annotation\t%s\t%zu\t%zu\t%.1lf\t%.1lf\t%.2lf\t%.2lf\n", read_name.c_str(), record->core.pos, best_partition, sample_start, sample_end, read_rate, polya_length);

    if(opt::verbose >= 1) {
        // reverse samples
        const std::vector<float> samples = sr.samples;

        for(size_t i = 0; i < 2 * sample_end && samples.size(); ++i) {
            std::string tag;
            if(i < sample_start) {
                tag = "PRETAIL";
            } else if(i < sample_end) {
                tag = "POLYA";
            } else {
                tag = "POSTTAIL";
            }

            double s = samples[i];
            double scaled_s = (s - sr.pore_model[0].shift) / sr.pore_model[0].scale;
            fprintf(out_fp, "polya-samples\t%s\t%zu\t%f\t%f\t%s\n", read_name.substr(0,6).c_str(), i, s, scaled_s, tag.c_str());
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
    //auto f = std::bind(estimate_polya_for_single_read, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    auto f = std::bind(estimate_polya_for_single_read_hmm, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    // free allocated values:
    fai_destroy(fai);

    return EXIT_SUCCESS;
}
