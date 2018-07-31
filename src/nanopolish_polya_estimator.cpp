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
"  -w, --window=STR                     only compute the poly-A lengths for reads in window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 1D ONT direct RNA reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome assembly for the reads is in FILE\n"
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
            case 'w': arg >> opt::region; break;
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

// ================================================================================
// Segmentation Hidden Markov Model
//   Define an HMM class `SegmentationHMM` with all relevant functions necessary for
//   segmentation of a squiggle into a series of regions.
// * struct ViterbiOutputs: contains log-prob scores and inferred state sequence
//   from a run of the viterbi algorithm.
// * struct Segmentation: contains ending sample indices for each region of a
//   squiggle's segmentation.
// * SegmentationHMM: class defining a hidden markov model for segmentation of
//   a squiggle. Contains the following members:
//   - state_transitions
//   - start_probs
//   - gaussian parameters defining emission distributions
//   - log-probabilities
//   + viterbi
//   + segment_squiggle
//   + log_probas
// ================================================================================
// struct ViterbiOutputs composed of viterbi probs
// and a vector of integers from {0,1,2,3,4,5} == {S,L,A,P,C,T}.
struct ViterbiOutputs {
    std::vector<float> scores;
    std::vector<uint8_t> labels;
};

// Segmentation struct holds endpoints of distinct regions from a segmented squiggle:
struct Segmentation {
    size_t start;   // final index of S; might not exist if skipped over
    size_t leader;  // final index of L, as indicated by 3'->5' viterbi
    size_t adapter; // final index of A, as indicated by 3'->5' viterbi
    size_t polya;   // final index of P/C, as indicated by 3'->5' viterbi
    size_t cliffs;  // number of observed 'CLIFF' samples
};

// Basic HMM struct with fixed parameters and viterbi/segmentation methods.
// (N.B.: all of the below is relative to a **scaled & shifted** set of events.)
class SegmentationHMM {
private:
    // ----- state space parameters:
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

    // ----- emission parameters:
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

    // log-probabilities are computed in the constructor:
    float log_state_transitions[6][6];
    float log_start_probs[6];

    // ----- inlined computation of emission log-probabilities:
    // Get the log-probability of seeing `x` given we're in state `state` of the HMM
    // N.B.: we scale the emission parameters (events are **not** scaled).
    inline float emit_log_proba(const float x, const uint8_t state) const
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
            log_probs = log_normal_pdf(xx, this->s_emission);
        }
        if (state == 1) {
            // LEADER state:
            log_probs = log_normal_pdf(xx, this->l_emission);
        }
        if (state == 2) {
            // ADAPTER state: compute log of gaussian mixture probability
            float mixture_proba = (this->a0_coeff*normal_pdf(xx,this->a0_emission)) + \
                (this->a1_coeff*normal_pdf(xx, this->a1_emission));
            log_probs = std::log(mixture_proba);
        }
        if (state == 3) {
            // POLYA state:
            log_probs = log_normal_pdf(xx, this->p_emission);
        }
        if (state == 4) {
            // CLIFF state: middle-out uniform distribution
            if ((xx > this->c_begin) && (xx <  this->c_end)) {
                log_probs = this->c_log_prob;
            } else {
                log_probs = -INFINITY;
            }
        }
        if (state == 5) {
            // TRANSCRIPT state: compute log of gaussian mixture probability
            float mixture_proba = (this->t0_coeff*normal_pdf(xx, this->t0_emission)) + \
                (this->t1_coeff*normal_pdf(xx, this->t1_emission));
            log_probs = std::log(mixture_proba);
        }
        return log_probs;
    }

public:
    // ----- constructor: compute logs of params & scale/shift
    SegmentationHMM(float scale, float shift, float var)
    {
        // - - - initialize log-probabilities:
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                if (this->state_transitions[i][j] > 0.00f) {
                    this->log_state_transitions[i][j] = std::log(this->state_transitions[i][j]);
                } else {
                    this->log_state_transitions[i][j] = -INFINITY;
                }
            }

            if (this->start_probs[i] > 0.00f) {
                this->log_start_probs[i] = std::log(this->start_probs[i]);
            } else {
                this->log_start_probs[i] = -INFINITY;
            }
        }
        // - - - update all gaussian parameters by scaling/shifting:
        // START emissions:
        this->s_emission.mean = shift + scale*(this->s_emission.mean);
        this->s_emission.stdv = var * this->s_emission.stdv;
        this->s_emission.log_stdv = std::log(this->s_emission.stdv);
        // LEADER emissions:
        this->l_emission.mean = shift + scale*(this->l_emission.mean);
        this->l_emission.stdv = var * this->l_emission.stdv;
        this->l_emission.log_stdv = std::log(this->l_emission.stdv);
	// ADAPTER emissions:
        this->a0_emission.mean = shift + scale*(this->a0_emission.mean);
        this->a0_emission.stdv = var * this->a0_emission.stdv;
        this->a0_emission.log_stdv = std::log(this->a0_emission.stdv);
        this->a1_emission.mean = shift + scale*(this->a1_emission.mean);
        this->a1_emission.stdv = var * this->a1_emission.stdv;
        this->a1_emission.log_stdv = std::log(this->a1_emission.stdv);
        // POLYA emissions:
        this->p_emission.mean = shift + scale*(this->p_emission.mean);
        this->p_emission.stdv = var * this->p_emission.stdv;
        this->p_emission.log_stdv = std::log(this->p_emission.stdv);
        // TRANSCRIPT emissions:
        this->t0_emission.mean = shift + scale*(this->t0_emission.mean);
        this->t0_emission.stdv = var * this->t0_emission.stdv;
        this->t0_emission.log_stdv = std::log(this->t0_emission.stdv);
        this->t1_emission.mean = shift + scale*(this->t1_emission.mean);
        this->t1_emission.stdv = var * this->t1_emission.stdv;
        this->t1_emission.log_stdv = std::log(this->t1_emission.stdv);
    }
    // ----- destructor: nothing to clean up
    ~SegmentationHMM() { }

    // ----- for a given sample value and shift/scale parameters, return log-probs for each state:
    std::vector<float> log_probas(const float x) const
    {
	std::vector<float> log_proba(6);
        for (uint8_t k = 0; k < 6; ++k) {
            log_proba[k] = this->emit_log_proba(x, k);
        }
        return log_proba;
    }

    // ----- viterbi-decoding of a squiggle into region labels:
    // N.B.1: viterbi decoding happens in the 3'->5' direction.
    // N.B.2: this algorithm takes place in log-space for numerical stability;
    // the `scores` variable refers to log-prob scores.
    ViterbiOutputs viterbi(const SquiggleRead& sr) const
    {
        // count of raw samples:
        size_t num_samples = sr.samples.size();

        // create/initialize viterbi scores and backpointers:
        std::vector<float> init_scores(6, -std::numeric_limits<float>::infinity()); // log(0.0) == -INFTY
        std::vector<uint8_t> init_bptrs(6, 6); // 6 == INIT symbol
        std::vector< std::vector<float> > viterbi_scores(num_samples, init_scores);
        std::vector< std::vector<uint8_t> > viterbi_bptrs(num_samples, init_bptrs);

        // forward viterbi pass; fill up backpointers:
        // weight initially distributed between S and L:
        viterbi_scores[0][0] = this->log_start_probs[0] + this->emit_log_proba(sr.samples[num_samples-1], 0);
        viterbi_scores[0][1] = this->log_start_probs[1] + this->emit_log_proba(sr.samples[num_samples-1], 1);
        for (size_t i = 1; i < num_samples; ++i) {
            // `t` moves from 3'->5' on the vector of samples, in opposite direction of `i`:
            size_t t = (num_samples-1-i);
            // get individual incoming state scores:
            float s_to_s = viterbi_scores.at(i-1)[0] + this->log_state_transitions[0][0];
            float s_to_l = viterbi_scores.at(i-1)[0] + this->log_state_transitions[0][1];
            float l_to_l = viterbi_scores.at(i-1)[1] + this->log_state_transitions[1][1];
            float l_to_a = viterbi_scores.at(i-1)[1] + this->log_state_transitions[1][2];
            float a_to_a = viterbi_scores.at(i-1)[2] + this->log_state_transitions[2][2];
            float a_to_p = viterbi_scores.at(i-1)[2] + this->log_state_transitions[2][3];
            float p_to_p = viterbi_scores.at(i-1)[3] + this->log_state_transitions[3][3];
            float p_to_c = viterbi_scores.at(i-1)[3] + this->log_state_transitions[3][4];
            float p_to_t = viterbi_scores.at(i-1)[3] + this->log_state_transitions[3][5];
            float c_to_c = viterbi_scores.at(i-1)[4] + this->log_state_transitions[4][4];
            float c_to_p = viterbi_scores.at(i-1)[4] + this->log_state_transitions[4][3];
            float t_to_t = viterbi_scores.at(i-1)[5] + this->log_state_transitions[5][5];

            // update the viterbi scores for each state at this timestep:
            viterbi_scores.at(i)[0] = s_to_s + this->emit_log_proba(sr.samples[t], 0);
            viterbi_scores.at(i)[1] = std::max(l_to_l, s_to_l) + this->emit_log_proba(sr.samples[t], 1);
            viterbi_scores.at(i)[2] = std::max(a_to_a, l_to_a) + this->emit_log_proba(sr.samples[t], 2);
            viterbi_scores.at(i)[3] = std::max(p_to_p, std::max(a_to_p, c_to_p)) + this->emit_log_proba(sr.samples[t], 3);
            viterbi_scores.at(i)[4] = std::max(c_to_c, p_to_c) + this->emit_log_proba(sr.samples[t], 4);
            viterbi_scores.at(i)[5] = std::max(p_to_t, t_to_t) + this->emit_log_proba(sr.samples[t], 5);

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

        // format as ViterbiOutputs struct and return:
        ViterbiOutputs output_vectors = { scores, regions };
        return output_vectors;
    }

    // ----- parse a squiggle's viterbi labels into a regional segmentation:
    Segmentation segment_squiggle(const SquiggleRead& sr) const
    {
	ViterbiOutputs viterbi_outs = this->viterbi(sr);

        // compute final sample indices of each region:
        std::vector<uint8_t>& region_labels = viterbi_outs.labels;

        // initial values for indices should preserve expected order:
        Segmentation ixs = { 0, 1, 2, 3, 0 };

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

        // set sensible (easy to QC-filter) default values if not all four detected;
        // S-end is always detected (min value == 0)
        if (ixs.leader == 1 || ixs.adapter == 2 || ixs.polya == 3) {
            ixs.leader = region_labels.size() - 3;
            ixs.adapter = region_labels.size() - 2;
            ixs.polya = region_labels.size() - 1;
        }
        return ixs;
    }
};

// ================================================================================
// Estimate the duration profile for a single read.
//   Estimate the underlying read rate.
// * estimate_eventalign_duration_profile : compute median read rate via collapsed-
//     duration event-alignment.
// * estimate_unaligned_duration_profile : compute median read rate via collapsed-
//     durations, without event-alignment.
// ================================================================================
// Compute a read-rate based on event-alignment, collapsed by consecutive 5mer identity
// (N.B.: deprecated; using non-eventaligned durations seems to work just as well
// while being faster to run.)
double estimate_eventalign_duration_profile(SquiggleRead& sr,
                                            const faidx_t* fai,
                                            const bam_hdr_t* hdr,
                                            const bam1_t* record,
                                            const size_t read_idx)
{
    EventAlignmentParameters params;
    params.sr = &sr;
    params.fai = fai;
    params.hdr = hdr;
    params.record = record;
    params.strand_idx = 0;
    params.read_idx = read_idx;

    std::vector<EventAlignment> alignment_output = align_read_to_ref(params);

    // collect durations, collapsing by k-mer
    std::vector<double> durations_per_kmer;

    size_t prev_ref_position = -1;
    for(const auto& ea : alignment_output) {
        float event_duration = sr.get_duration(ea.event_idx, ea.strand_idx);
        size_t ref_position = ea.ref_position;
        if(ref_position == prev_ref_position) {
            assert(!durations_per_kmer.empty());
            durations_per_kmer.back() += event_duration;
        } else {
            durations_per_kmer.push_back(event_duration);
            prev_ref_position = ref_position;
        }
    }
    std::sort(durations_per_kmer.begin(), durations_per_kmer.end());
    double median_duration = durations_per_kmer[durations_per_kmer.size() / 2];

    // this is our estimator of read rate, currently we use the median duration
    // per k-mer as its more robust to outliers caused by stalls
    double read_rate = 1.0 / median_duration;

    return read_rate;
}

// compute a read-rate based on kmer-to-event mapping, collapsed by consecutive 5mer identity:
double estimate_unaligned_duration_profile(const SquiggleRead& sr,
                                           const faidx_t* fai,
                                           const bam_hdr_t* hdr,
                                           const bam1_t* record,
                                           const size_t read_idx,
                                           const size_t strand_idx)
{
    // get kmer stats:
    size_t basecalled_k = sr.get_base_model(strand_idx)->k;
    size_t num_kmers = sr.read_sequence.length() - basecalled_k + 1;

    // collect durations, collapsing by k-mer:
    std::vector<double> durations_per_kmer(num_kmers);
    for (size_t i = 0; i < sr.base_to_event_map.size(); ++i) {
        size_t start_idx = sr.base_to_event_map[i].indices[strand_idx].start;
        size_t end_idx = sr.base_to_event_map[i].indices[strand_idx].stop;
        // no events for this k-mer
	if (start_idx == -1) {
            continue;
        }
        assert(start_idx <= end_idx);
        for (size_t j = start_idx; j <= end_idx; ++j) {
            durations_per_kmer[i] += sr.get_duration(j, strand_idx);
        }
    }
    
    std::sort(durations_per_kmer.begin(), durations_per_kmer.end());
    assert(durations_per_kmer.size() > 0);
    double median_duration = durations_per_kmer[durations_per_kmer.size() / 2];

    // this is our estimator of read rate, currently we use the median duration
    // per k-mer as its more robust to outliers caused by stalls
    double read_rate = 1.0 / median_duration;

    return read_rate;
}

// fetch the raw event durations for a given read:
std::vector<double> fetch_event_durations(const SquiggleRead& sr,
                                          const faidx_t* fai,
                                          const bam_hdr_t* hdr,
                                          const bam1_t* record,
                                          const size_t read_idx,
                                          const size_t strand_idx)
{
    // get kmer stats:
    size_t basecalled_k = sr.get_base_model(strand_idx)->k;
    size_t num_kmers = sr.read_sequence.length() - basecalled_k + 1;

    // collect durations, collapsing by k-mer:
    std::vector<double> durations_per_kmer(num_kmers);
    for (size_t i = 0; i < sr.base_to_event_map.size(); ++i) {
        size_t start_idx = sr.base_to_event_map[i].indices[strand_idx].start;
        size_t end_idx = sr.base_to_event_map[i].indices[strand_idx].stop;
        // no events for this k-mer
        if (start_idx == -1) {
            continue;
        }
        assert(start_idx <= end_idx);
        for (size_t j = start_idx; j <= end_idx; ++j) {
            durations_per_kmer[i] += sr.get_duration(j, strand_idx);
        }
    }
    assert(durations_per_kmer.size() > 0);

    return durations_per_kmer;
}

// ================================================================================
// Poly-A Tail Length Estimation
//   Estimate the number of nucleotides in the poly-A region.
// * estimate_polya_length : return an estimate of the read rate for this read.
// ================================================================================
// Compute an estimate of the number of nucleotides in the poly-A tail
double estimate_polya_length(const SquiggleRead& sr, const Segmentation& region_indices, const double read_rate)
{
    // start and end times (sample indices) of the poly(A) tail, in original 3'->5' time-direction:
    // (n.b.: everything in 5'->3' order due to inversion in SquiggleRead constructor, but our
    // `region_indices` struct has everything in 3'->5' order)
    double num_samples = sr.samples.size();
    double polya_sample_start = region_indices.adapter + 1;
    double polya_sample_end = region_indices.polya;
    double adapter_sample_start = region_indices.leader;
    double leader_sample_start = region_indices.start;

    // calculate duration of poly(A) region (in seconds)
    double polya_duration = (region_indices.polya - (region_indices.adapter + 1)) / sr.sample_rate;

    // Empirically determined offset to handle modal bias of the estimator:
    double estimation_error_offset = -5;

    // length of the poly(A) tail, in nucleotides:
    double polya_length = polya_duration * read_rate + estimation_error_offset;

    // ensure estimated length is non-negative:
    polya_length = std::max(0.0, polya_length);

    return polya_length;
}

// ================================================================================
// QC Functions: pre-segmentation, post-segmentation, post-estimation
//
// * pre_segmentation_qc: return a bool (true ~> FAIL) performing basic pre-seg QC.
// * post_segmentation_qc: check the segmentation results for failures.
// * post_estimation_qc: sanity check for estimates.
// ================================================================================
// Some basic sanity-check QC before segmentation; this returns true if QC-FAIL and
// false if QC-PASS.
bool pre_segmentation_qc(uint32_t suffix_clip, uint32_t prefix_clip, double transcript_length, const SquiggleRead& sr)
{
    // skip this read if long skip at end or if most of transcript wasnt aligned:
    if (suffix_clip > 200 || (double)(prefix_clip + suffix_clip) / transcript_length > 0.2) {
        return true;
    }
    // skip if no events:
    if (sr.events[0].empty()) {
        return true;
    }
    // return false (i.e. dont skip) if none of the above fail:
    return false;
}

// QC pass after constructing a segmentation; returns a QC flag represented as a string,
// either "PASS" or "NOREGION".
// These tests indicate that something went wrong in the segmentation algorithm.
std::string post_segmentation_qc(const Segmentation& region_indices, const SquiggleRead& sr)
{
    // fetch sizes of LEADER, ADAPTER, POLYA regions:
    double num_leader_samples = region_indices.leader;
    double num_adapter_samples = (region_indices.adapter+1) - region_indices.leader;
    double num_polya_samples = region_indices.polya - (region_indices.adapter+1);

    // check for NOREGION:
    std::string qc_tag;
    if (num_leader_samples < 200.0 || num_adapter_samples < 200.0 || num_polya_samples < 200.0) {
        qc_tag = "NOREGION";
    } else {
        qc_tag = "PASS";
    }
    return qc_tag;
}

// QC pass after performing estimation; returns a QC flag represented as a string.
// Currently returns either "PASS" or "ADAPTER".
// These tests indicate that something went wrong in the estimation algorithm.
std::string post_estimation_qc(const Segmentation& region_indices, const SquiggleRead& sr, double read_rate, double polya_length)
{
    // `adapter_qc_tol` is the (empirically-discovered) upper-tolerance for number of estimated adapter nucleotides:
    double adapter_qc_tol = 300.0f;
    // estimated adapter length, in nucleotides:
    double adapter_duration = (region_indices.adapter - (region_indices.leader - 1)) / sr.sample_rate;
    double adapter_length = adapter_duration * read_rate;

    std::string qc_tag;
    if (adapter_length > adapter_qc_tol) {
        qc_tag = "ADAPTER";
    } else {
        qc_tag = "PASS";
    }
    return qc_tag;
}

// ================================================================================
// Main Poly-A Code
//   Expose main functionality of this module via public-facing function.
// * estimate_polya_for_single_read : performs all data-fetching and joins all
//     of the above functions; writes results to file in TSV format
// * polya_main : wrap `estimate_polya_for_single_read` with OpenMP directives
//     for easy multi-threading across reads
// ================================================================================
// Write Poly(A) region segmentation and tail length estimation data to TSV
void estimate_polya_for_single_read(const ReadDB& read_db,
                                    const faidx_t* fai,
                                    FILE* out_fp,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{
    //----- load a squiggle read:
    std::string read_name = bam_get_qname(record);
    std::string ref_name(hdr->target_name[record->core.tid]);
    size_t strand_idx = 0;

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

    //----- Perform pre-segmentation QC:
    if (pre_segmentation_qc(suffix_clip, prefix_clip, sequenced_transcript.length(), sr)) {
        return;
    }

    //----- perform HMM-based regional segmentation & post-segmentation QC:
    SegmentationHMM hmm(static_cast<float>(sr.scalings[0].scale),
                        static_cast<float>(sr.scalings[0].shift),
                        static_cast<float>(sr.scalings[0].var));
    Segmentation region_indices = hmm.segment_squiggle(sr);
    std::string post_segmentation_qc_flag = post_segmentation_qc(region_indices, sr);

    //----- compute duration profile for the read:
    double read_rate = estimate_unaligned_duration_profile(sr, fai, hdr, record, read_idx, strand_idx);

    //----- estimate number of nucleotides in poly-A tail & post-estimation QC:
    double polya_length = estimate_polya_length(sr, region_indices, read_rate);
    std::string post_estimation_qc_flag = post_estimation_qc(region_indices, sr, read_rate, polya_length);

    //----- Resolve QC flag based on priority:
    std::string qc_tag;
    if (post_segmentation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_segmentation_qc_flag;
    } else if (post_estimation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_estimation_qc_flag;
    } else {
        qc_tag = "PASS";
    }

    //----- print annotations to TSV:
    double leader_sample_start = region_indices.start+1;
    double adapter_sample_start = region_indices.leader+1;
    double polya_sample_start = region_indices.adapter+1;    
    double polya_sample_end = region_indices.polya;
    double transcr_sample_start = region_indices.polya+1;
    #pragma omp critical
    {
        fprintf(out_fp, "%s\t%s\t%zu\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%s\n",
                read_name.c_str(), ref_name.c_str(), record->core.pos,
                leader_sample_start, adapter_sample_start, polya_sample_start,
                transcr_sample_start, read_rate, polya_length, qc_tag.c_str());
        // if `verbose == 1`, print the samples (picoAmps) of the read,
        // up to the first 1000 samples of transcript region:
        if (opt::verbose == 1) {
            // copy 5'->3'-oriented samples from SquiggleRead and reverse back to 3'->5':
            std::vector<float> samples(sr.samples);
            std::reverse(samples.begin(), samples.end());
            const SegmentationHMM vhmm(static_cast<float>(sr.scalings[0].scale),
                                       static_cast<float>(sr.scalings[0].shift),
                                       static_cast<float>(sr.scalings[0].var));
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
                float s = samples[i];
                float scaled_s = (s - sr.scalings[0].shift) / sr.scalings[0].scale;
                std::vector<float> s_probas = vhmm.log_probas(s);
                fprintf(out_fp, "polya-samples\t%s\t%s\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                        read_name.substr(0,6).c_str(), ref_name.c_str(), i, s, scaled_s,
                        s_probas.at(0), s_probas.at(1), s_probas.at(2), s_probas.at(3), s_probas.at(4), s_probas.at(5),
                        tag.c_str());
            }
        }
        // if `verbose == 2`, print the raw event durations of the read:
        if (opt::verbose >= 2) {
            std::vector<double> raw_durations = fetch_event_durations(sr, fai, hdr, record, read_idx, strand_idx);
            for (size_t i = 0; i < raw_durations.size(); ++i) {
                double dura = raw_durations[i];
                fprintf(out_fp, "polya-durations\t%s\t%zu\t%f\t%s\n", read_name.c_str(), i, dura, qc_tag.c_str());
            }
        }
    }
}

// Wrap poly-A estimation code for parallelism
int polya_main(int argc, char** argv)
{
    parse_polya_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // print header line:
    fprintf(stdout, "readname\tcontig\tpos\tleader_start\tadapter_start\tpolya_start\ttranscr_start\tread_rate\tpolya_length\tqc_tag\n");

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(estimate_polya_for_single_read, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    // free allocated values:
    fai_destroy(fai);

    return EXIT_SUCCESS;
}
