//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_detect_polyi.cpp -- detect the presence of a
// poly(I) tail as in the nano-COP protocol.
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
#include "nanopolish_detect_polyi.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_emissions.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"


using namespace std::placeholders;

//
// Getopt
//
#define SUBPROGRAM "detect-polyi"

static const char *DETECT_POLYI_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *DETECT_POLYI_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Detect presence of poly(I) tails and estimate length of tails in direct RNA reads\n"
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

void parse_detect_polyi_options(int argc, char** argv)
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
                std::cout << DETECT_POLYI_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << DETECT_POLYI_VERSION_MESSAGE;
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
        std::cout << "\n" << DETECT_POLYI_USAGE_MESSAGE;
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
// Segmentation struct holds endpoints of distinct regions from a segmented squiggle:
struct DPISegmentation {
    size_t start;   // final index of S; might not exist if skipped over
    size_t leader;  // final index of L, as indicated by 3'->5' viterbi
    size_t adapter; // final index of A, as indicated by 3'->5' viterbi
    size_t polya;   // final index of P/C, as indicated by 3'->5' viterbi
    size_t cliffs;  // number of observed 'CLIFF' samples
};

// Basic HMM struct with fixed parameters and viterbi/segmentation methods.
// (N.B.: all of the below is relative to a **scaled & shifted** set of events.)
enum DPIHMMState
{
    DPI_HMM_START = 0,
    DPI_HMM_LEADER = 1,
    DPI_HMM_ADAPTER = 2,
    DPI_HMM_POLYA = 3,
    DPI_HMM_CLIFF = 4,
    DPI_HMM_TRANSCRIPT = 5,
    DPI_HMM_NUM_STATES = 6 // number of non-NULL states in HMM
};

// struct ViterbiOutputs composed of viterbi probs
// and a vector of integers from {0,1,2,3,4,5} == {S,L,A,P,C,T}.
struct DPIViterbiOutputs {
    std::vector<float> scores;
    std::vector<DPIHMMState> labels;
};

class DPISegmentationHMM {
private:
    // ----- state space parameters:
    // N.B.: `state transitions` is used to compute log probabilities, as viterbi decoding is done in log-space.
    // state transition probabilities (S->L->A->[P<->C]->T):
    float state_transitions[DPI_HMM_NUM_STATES][DPI_HMM_NUM_STATES] = {
        // S -> S (10%), S -> L (90%)
        {0.10f, 0.90f, 0.00f, 0.00f, 0.00f, 0.00f},
        // L -> A (10%), L -> L (90%)
        {0.00f, 0.90f, 0.10f, 0.00f, 0.00f, 0.00f},
        // A -> P (05%), A -> A (95%)
        {0.00f, 0.00f, 0.95f, 0.05f, 0.00f, 0.00f},
        // P -> P (89%), P -> C (01%), P -> T (10%)
        {0.00f, 0.00f, 0.00f, 0.89f, 0.01f, 0.10f},
        // C -> P (99%), C -> C (01%)
        {0.00f, 0.00f, 0.00f, 0.99f, 0.01f, 0.00f},
        // T -> T (100%)
        {0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f}
    };
    // All state sequences must start on S:
    float start_probs[DPI_HMM_NUM_STATES] = { 1.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f };

    // ----- emission parameters:
    // emission parameters, from empirical MLE on manually-flagged reads:
    // START has a mixture of Gaussian and Uniform emissions;
    // LEADER has a Gaussian emission;
    // ADAPTER, POLYA, TRANSCRIPT have Gaussian mixture emissions;
    // CLIFF has Uniform emissions.
    GaussianParameters s_emission = {70.2737f, 3.7743f};
    float s_begin = 40.0f;
    float s_end = 250.0f;
    float s_prob = 0.00476f; // == {1. / (250.0f - 40.0f)}
    float s_norm_coeff = 0.50f;
    float s_unif_coeff = 0.50f;
    GaussianParameters l_emission = {110.973f, 5.237f};
    GaussianParameters a0_emission = {79.347f, 8.3702f};
    GaussianParameters a1_emission = {63.3126f, 2.7464f};
    float a0_coeff = 0.874f;
    float a1_coeff = 0.126f;
    GaussianParameters p0_emission = {108.883f, 3.257f};
    GaussianParameters p1_emission = {108.498f, 5.257f};
    float p0_coeff = 0.500f;
    float p1_coeff = 0.500f;
    float c_begin = 70.0f;
    float c_end = 140.0f;
    float c_log_prob = -4.2485f; // natural log of [1/(140-70)]
    GaussianParameters t0_emission = {79.679f, 6.966f};
    GaussianParameters t1_emission = {105.784f, 16.022f};
    float t0_coeff = 0.346f;
    float t1_coeff = 0.654f;

    // log-probabilities are computed in the constructor:
    float log_state_transitions[DPI_HMM_NUM_STATES][DPI_HMM_NUM_STATES];
    float log_start_probs[DPI_HMM_NUM_STATES];

    // ----- inlined computation of emission log-probabilities:
    // Get the log-probability of seeing `x` given we're in state `state` of the HMM
    // N.B.: we scale the emission parameters (events are **not** scaled).
    inline float emit_log_proba(const float x, const DPIHMMState state) const
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
        if (state == DPI_HMM_START) {
            // START state:
            float norm_term = s_norm_coeff * normal_pdf(xx, this->s_emission);
            log_probs = std::log(norm_term + s_unif_coeff * s_prob);
        }
        if (state == DPI_HMM_LEADER) {
            // LEADER state:
            log_probs = log_normal_pdf(xx, this->l_emission);
        }
        if (state == DPI_HMM_ADAPTER) {
            // ADAPTER state: compute log of gaussian mixture probability
            float mixture_proba = (this->a0_coeff*normal_pdf(xx,this->a0_emission)) + \
                (this->a1_coeff*normal_pdf(xx, this->a1_emission));
            log_probs = std::log(mixture_proba);
        }
        if (state == DPI_HMM_POLYA) {
            // POLYA state:
	    float mixture_proba = (this->p0_coeff*normal_pdf(xx, this->p0_emission)) + \
		(this->p1_coeff*normal_pdf(xx, this->p1_emission));
            log_probs = std::log(mixture_proba);
        }
        if (state == DPI_HMM_CLIFF) {
            // CLIFF state: middle-out uniform distribution
            if ((xx > this->c_begin) && (xx <  this->c_end)) {
                log_probs = this->c_log_prob;
            } else {
                log_probs = -INFINITY;
            }
        }
        if (state == DPI_HMM_TRANSCRIPT) {
            // TRANSCRIPT state: compute log of gaussian mixture probability
            float mixture_proba = (this->t0_coeff*normal_pdf(xx, this->t0_emission)) + \
                (this->t1_coeff*normal_pdf(xx, this->t1_emission));
            log_probs = std::log(mixture_proba);
        }
        return log_probs;
    }

public:
    // ----- constructor: compute logs of params & scale/shift
    DPISegmentationHMM(float scale, float shift, float var)
    {
        // - - - initialize log-probabilities:
        for (int i = 0; i < DPI_HMM_NUM_STATES; ++i) {
            for (int j = 0; j < DPI_HMM_NUM_STATES; ++j) {
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
        this->p0_emission.mean = shift + scale*(this->p0_emission.mean);
        this->p0_emission.stdv = var * this->p0_emission.stdv;
        this->p0_emission.log_stdv = std::log(this->p0_emission.stdv);
        this->p1_emission.mean = shift + scale*(this->p1_emission.mean);
        this->p1_emission.stdv = var * this->p1_emission.stdv;
        this->p1_emission.log_stdv = std::log(this->p1_emission.stdv);
        // TRANSCRIPT emissions:
        this->t0_emission.mean = shift + scale*(this->t0_emission.mean);
        this->t0_emission.stdv = var * this->t0_emission.stdv;
        this->t0_emission.log_stdv = std::log(this->t0_emission.stdv);
        this->t1_emission.mean = shift + scale*(this->t1_emission.mean);
        this->t1_emission.stdv = var * this->t1_emission.stdv;
        this->t1_emission.log_stdv = std::log(this->t1_emission.stdv);
    }
    // ----- destructor: nothing to clean up
    ~DPISegmentationHMM() { }

    // ----- for a given sample value and shift/scale parameters, return log-probs for each state:
    std::vector<float> log_probas(const float x) const
    {
        std::vector<float> log_proba(DPI_HMM_NUM_STATES);
        for (uint8_t k = 0; k < DPI_HMM_NUM_STATES; ++k) {
            log_proba[k] = this->emit_log_proba(x, static_cast<DPIHMMState>(k));
        }
        return log_proba;
    }

    // ----- viterbi-decoding of a squiggle into region labels:
    // N.B.1: viterbi decoding happens in the 3'->5' direction.
    // N.B.2: this algorithm takes place in log-space for numerical stability;
    // the `scores` variable refers to log-prob scores.
    DPIViterbiOutputs viterbi(const SquiggleRead& sr) const
    {
        // count of raw samples:
        size_t num_samples = sr.samples.size();

        // create/initialize viterbi scores and backpointers:
        std::vector<float> init_scores(DPI_HMM_NUM_STATES, -std::numeric_limits<float>::infinity()); // log(0.0) == -INFTY
        std::vector<DPIHMMState> init_bptrs(DPI_HMM_NUM_STATES, DPI_HMM_NUM_STATES); // HMM_NUM_STATES used as a dummy value here
        std::vector< std::vector<float> > viterbi_scores(num_samples, init_scores);
        std::vector< std::vector<DPIHMMState> > viterbi_bptrs(num_samples, init_bptrs);

        // forward viterbi pass; fill up backpointers:
        // weight initially distributed between START and LEADER:
        viterbi_scores[0][DPI_HMM_START] = this->log_start_probs[DPI_HMM_START] + this->emit_log_proba(sr.samples[num_samples-1], DPI_HMM_START);
        viterbi_scores[0][DPI_HMM_LEADER] = this->log_start_probs[DPI_HMM_LEADER] + this->emit_log_proba(sr.samples[num_samples-1], DPI_HMM_LEADER);
        for (size_t i = 1; i < num_samples; ++i) {
            // get individual incoming state scores:
            float s_to_s = viterbi_scores.at(i-1)[DPI_HMM_START] + this->log_state_transitions[DPI_HMM_START][DPI_HMM_START];
            float s_to_l = viterbi_scores.at(i-1)[DPI_HMM_START] + this->log_state_transitions[DPI_HMM_START][DPI_HMM_LEADER];
            float l_to_l = viterbi_scores.at(i-1)[DPI_HMM_LEADER] + this->log_state_transitions[DPI_HMM_LEADER][DPI_HMM_LEADER];
            float l_to_a = viterbi_scores.at(i-1)[DPI_HMM_LEADER] + this->log_state_transitions[DPI_HMM_LEADER][DPI_HMM_ADAPTER];
            float a_to_a = viterbi_scores.at(i-1)[DPI_HMM_ADAPTER] + this->log_state_transitions[DPI_HMM_ADAPTER][DPI_HMM_ADAPTER];
            float a_to_p = viterbi_scores.at(i-1)[DPI_HMM_ADAPTER] + this->log_state_transitions[DPI_HMM_ADAPTER][DPI_HMM_POLYA];
            float p_to_p = viterbi_scores.at(i-1)[DPI_HMM_POLYA] + this->log_state_transitions[DPI_HMM_POLYA][DPI_HMM_POLYA];
            float p_to_c = viterbi_scores.at(i-1)[DPI_HMM_POLYA] + this->log_state_transitions[DPI_HMM_POLYA][DPI_HMM_CLIFF];
            float p_to_t = viterbi_scores.at(i-1)[DPI_HMM_POLYA] + this->log_state_transitions[DPI_HMM_POLYA][DPI_HMM_TRANSCRIPT];
            float c_to_c = viterbi_scores.at(i-1)[DPI_HMM_CLIFF] + this->log_state_transitions[DPI_HMM_CLIFF][DPI_HMM_CLIFF];
            float c_to_p = viterbi_scores.at(i-1)[DPI_HMM_CLIFF] + this->log_state_transitions[DPI_HMM_CLIFF][DPI_HMM_POLYA];
            float t_to_t = viterbi_scores.at(i-1)[DPI_HMM_TRANSCRIPT] + this->log_state_transitions[DPI_HMM_TRANSCRIPT][DPI_HMM_TRANSCRIPT];

            // update the viterbi scores for each state at this timestep:
            viterbi_scores.at(i)[DPI_HMM_START] = s_to_s + this->emit_log_proba(sr.samples[i], DPI_HMM_START);
            viterbi_scores.at(i)[DPI_HMM_LEADER] = std::max(l_to_l, s_to_l) + this->emit_log_proba(sr.samples[i], DPI_HMM_LEADER);
            viterbi_scores.at(i)[DPI_HMM_ADAPTER] = std::max(a_to_a, l_to_a) + this->emit_log_proba(sr.samples[i], DPI_HMM_ADAPTER);
            viterbi_scores.at(i)[DPI_HMM_POLYA] = std::max(p_to_p, std::max(a_to_p, c_to_p)) + this->emit_log_proba(sr.samples[i], DPI_HMM_POLYA);
            viterbi_scores.at(i)[DPI_HMM_CLIFF] = std::max(c_to_c, p_to_c) + this->emit_log_proba(sr.samples[i], DPI_HMM_CLIFF);
            viterbi_scores.at(i)[DPI_HMM_TRANSCRIPT] = std::max(p_to_t, t_to_t) + this->emit_log_proba(sr.samples[i], DPI_HMM_TRANSCRIPT);

            // backpointers:
            // START: S can only come from S
            viterbi_bptrs.at(i)[DPI_HMM_START] = DPI_HMM_START;
            // LEADER: L->L or S->L
            if (s_to_l < l_to_l) {
                viterbi_bptrs.at(i)[DPI_HMM_LEADER] = DPI_HMM_LEADER;
            } else {
                viterbi_bptrs.at(i)[DPI_HMM_LEADER] = DPI_HMM_START;
            }
            // ADAPTER:
            if (l_to_a < a_to_a) {
                viterbi_bptrs.at(i)[DPI_HMM_ADAPTER] = DPI_HMM_ADAPTER;
            } else {
                viterbi_bptrs.at(i)[DPI_HMM_ADAPTER] = DPI_HMM_LEADER;
            }
            // POLYA:
            if ((a_to_p < p_to_p) && (c_to_p < p_to_p)) {
                viterbi_bptrs.at(i)[DPI_HMM_POLYA] = DPI_HMM_POLYA;
            } else if ((p_to_p < a_to_p) && (c_to_p < a_to_p)) {
                viterbi_bptrs.at(i)[DPI_HMM_POLYA] = DPI_HMM_ADAPTER;
            } else {
                viterbi_bptrs.at(i)[DPI_HMM_POLYA] = DPI_HMM_CLIFF;
            }
            // CLIFF:
            if (p_to_c < c_to_c) {
                viterbi_bptrs.at(i)[DPI_HMM_CLIFF] = DPI_HMM_CLIFF;
            } else {
                viterbi_bptrs.at(i)[DPI_HMM_CLIFF] = DPI_HMM_POLYA;
            }
            // TRANSCRIPT:
            if (p_to_t < t_to_t) {
                viterbi_bptrs.at(i)[DPI_HMM_TRANSCRIPT] = DPI_HMM_TRANSCRIPT;
            } else {
                viterbi_bptrs.at(i)[DPI_HMM_TRANSCRIPT] = DPI_HMM_POLYA;
            }
        }

        // backwards viterbi pass:
        // allocate `regions` vector of same dimensions as sample sequence;
        // clamp final state to 'T' ~ transcript:
        std::vector<DPIHMMState> regions(num_samples, DPI_HMM_START);
        std::vector<float> scores(num_samples, 0);
        regions[num_samples-1] = DPI_HMM_TRANSCRIPT;
        scores[num_samples-1] = viterbi_scores.at(num_samples-1)[DPI_HMM_TRANSCRIPT];
        // loop backwards and keep appending best states:
        for (size_t j=(num_samples-2); j > 0; --j) {
            regions[j] = viterbi_bptrs.at(j)[regions.at(j+1)];
            scores[j] = viterbi_scores.at(j)[regions.at(j+1)];
        }

        // format as DPIViterbiOutputs struct and return:
        DPIViterbiOutputs output_vectors = { scores, regions };
        return output_vectors;
    }

    // ----- parse a squiggle's viterbi labels into a regional segmentation:
    DPISegmentation segment_squiggle(const SquiggleRead& sr) const
    {
        DPIViterbiOutputs viterbi_outs = this->viterbi(sr);

        // compute final sample indices of each region:
        std::vector<DPIHMMState>& region_labels = viterbi_outs.labels;

        // initial values for indices should preserve expected order:
        DPISegmentation ixs = { 0, 1, 2, 3, 0 };

        // loop through sequence and collect values:
        for (std::vector<uint8_t>::size_type i = 0; i < region_labels.size(); ++i) {
            // call end of START:
            if (region_labels[i] == DPI_HMM_START && region_labels[i+1] == DPI_HMM_LEADER) {
                ixs.start = static_cast<size_t>(i);
            }
            // call end of leader:
            if (region_labels[i] == DPI_HMM_LEADER && region_labels[i+1] == DPI_HMM_ADAPTER) {
                ixs.leader = static_cast<size_t>(i);
            }
            // call end of adapter:
            if (region_labels[i] == DPI_HMM_ADAPTER && region_labels[i+1] == DPI_HMM_POLYA) {
                ixs.adapter = static_cast<size_t>(i);
            }
            // call end of polya:
            if (region_labels[i] == DPI_HMM_POLYA && region_labels[i+1] == DPI_HMM_TRANSCRIPT) {
                ixs.polya = static_cast<size_t>(i);
            }
            // increment cliff counter:
            if (region_labels[i] == DPI_HMM_CLIFF) {
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
// Bernoulli Hidden Markov Model
//   The `BernoulliHMM` class is a two-state hidden markov model designed to find
//   a (potentially nonexistent) switchpoint in a region with two very similar
//   Gaussian emissions by discretizing the log-likelihood ratio.
// * struct BernoulliOutputs: contains log-prob scores and inferred state sequence
//   from a run of the viterbi algorithm on the bernoulli-distributed sequence.
// * struct BernoulliSegmentation: contains ending sample indices for each region
//   of a squiggle's segmentation.
// * BernoulliHMM: class defining a hidden markov model for segmentation of
//   a squiggle.
// ================================================================================
// Contains the (possibly nonexistent) locations of the switchpoints.
struct BernoulliSegmentation {
    int polyi;    // **LAST** index of p(I) region (or -1 if not found).
    int polya;    // **FIRST** index of p(A) region (or -1 if not found).
};

// Descriptive shorthands for the states of the bernoulli HMM.
enum BernoulliState {
    BERN_POLYI = 0,
    BERN_POLYA = 1,
    BERN_NUM_STATES = 2
};

// Container struct for the output sequences of the viterbi algorithm on the model.
struct BernoulliOutputs {
    std::vector<float> scores;
    std::vector<BernoulliState> labels;
};

class BernoulliHMM {
private:
    float state_transitions[BERN_NUM_STATES][BERN_NUM_STATES] = {
        {0.90f, 0.10f},
        {0.00f, 1.00f}
    };
    // initial state probabilities:
    float start_probs[BERN_NUM_STATES] = { 1.00f, 0.00f };


    // Normal distributions for poly(I) and poly(A) states:
    // (N.B.: these should *not* be scaled by the SquiggleRead's linear parameters.)
    GaussianParameters pI_gaussian = {108.498f, 5.257f};
    GaussianParameters pA_gaussian = {108.883f, 3.257f};

    // Global mean sample value for recentering:
    float global_mean = 108.000f;

    // Bernoulli parameters for poly(I) and poly(A) binary log-likelihood ratios:
    // (These represent the probability of observing a 1 in a {0,1}-supported bernoulli)
    float pI_bernoulli = 0.72304f;
    float pA_bernoulli = 0.92154f;

    // log-probabilities are computed in the constructor:
    float log_state_transitions[BERN_NUM_STATES][BERN_NUM_STATES];
    float log_start_probs[BERN_NUM_STATES];

    // ----- inlined computation of emission log-probabilities:
    // Get the log-probability of seeing `x` given we're in state `state` of the HMM
    // N.B.: we scale the emission parameters (events are **not** scaled).
    inline float emit_log_proba_gaussian(const float x, const BernoulliState state) const
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
        if (state == BERN_POLYI) {
            // POLY(I) state:
            log_probs = std::log(normal_pdf(xx, this->pI_gaussian));
        }
        if (state == BERN_POLYA) {
            // POLY(A) state:
            log_probs = std::log(normal_pdf(xx, this->pA_gaussian));
        }
        return log_probs;
    }

    // ----- inlined computation of emission log-probabilities, for bernoulli distributions:
    inline float emit_log_proba_bernoulli(const int n, const BernoulliState state) const
    {
        float log_probs;
        if (state == BERN_POLYI) {
            if (n == 1) {
                log_probs = std::log(this->pI_bernoulli);
            } else {
                log_probs = std::log(1.0f - this->pI_bernoulli);
            }
        }
        if (state == BERN_POLYA) {
            if (n == 1) {
                log_probs = std::log(this->pA_bernoulli);
            } else {
                log_probs = std::log(1.0f - this->pA_bernoulli);
            }
        }
        return log_probs;
    }

public:
    // Constructor method for BernoulliHMM
    BernoulliHMM(float scale, float shift, float var)
    {
        // initialize log probabilities:
        for (int i = 0; i < BERN_NUM_STATES; ++i) {
            for (int j = 0; j < BERN_NUM_STATES; ++j) {
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
    }

    // Destructor method: nothing to do
    ~BernoulliHMM() { }

    // Take a vector of floats and return a vector of {0,1}-valued ints
    // by computing log-likelihood ratios { llkd(polyI) / llkd(polyA) } for each float.
    std::vector<int> log_lkhd_ratio_sequence(const std::vector<float>& signal) const
    {
        // --- compute mean pico-ampere value for this particular sequence:
        float instance_mean = 0.00f;
        for (size_t i = 0; i < signal.size(); ++i) {
            instance_mean = instance_mean + signal.at(i);
        }
        instance_mean = instance_mean / static_cast<float>(signal.size());

        // re-center this sequence to the global mean picoamp value and then binarize the log-likelihood ratio values:
        std::vector<int> bernoullis(signal.size());
        for (size_t i = 0; i < signal.size(); ++i) {
            float s = signal.at(i) - instance_mean + this->global_mean;
            float loglkhd_pI = this->emit_log_proba_gaussian(s, BERN_POLYI);
            float loglkhd_pA = this->emit_log_proba_gaussian(s, BERN_POLYA);
            if ((loglkhd_pI / loglkhd_pA) > 1.0f) {
                bernoullis.at(i) = 1;
            } else {
                bernoullis.at(i) = 0;
            }
        }

        return bernoullis;
    }

    // Viterbi implementation for BernoulliHMM with a 0-1 signal as input.
    BernoulliOutputs viterbi(const std::vector<int>& bernoullis) const {
        // --- Initialize viterbi score and backpointer vectors:
        std::vector<float> _init_scores(BERN_NUM_STATES, -std::numeric_limits<float>::infinity());
        std::vector<BernoulliState> _init_bptrs(BERN_NUM_STATES, BERN_NUM_STATES);
        std::vector< std::vector<float> > viterbi_scores(bernoullis.size(), _init_scores);
        std::vector< std::vector<BernoulliState> > viterbi_bptrs(bernoullis.size(), _init_bptrs);

        // --- forward viterbi pass: compute viterbi scores & backpointers.
        viterbi_scores[0][BERN_POLYI] = this->log_start_probs[BERN_POLYI] + this->emit_log_proba_bernoulli(bernoullis[0], BERN_POLYI);
        viterbi_scores[0][BERN_POLYA] = this->log_start_probs[BERN_POLYA] + this->emit_log_proba_bernoulli(bernoullis[0], BERN_POLYA);
        for (size_t i = 1; i < bernoullis.size(); ++i) {
            // compute transition probabilities:
            float i2i = viterbi_scores.at(i-1)[BERN_POLYI] + this->log_state_transitions[BERN_POLYI][BERN_POLYI];
            float i2a = viterbi_scores.at(i-1)[BERN_POLYI] + this->log_state_transitions[BERN_POLYI][BERN_POLYA];
            float a2a = viterbi_scores.at(i-1)[BERN_POLYA] + this->log_state_transitions[BERN_POLYA][BERN_POLYA];

            // update viterbi scores:
            viterbi_scores.at(i)[BERN_POLYI] = i2i + this->emit_log_proba_bernoulli(bernoullis[i], BERN_POLYI);
            viterbi_scores.at(i)[BERN_POLYA] = std::max(i2a, a2a) + this->emit_log_proba_bernoulli(bernoullis[i], BERN_POLYA);

            // update backpointers:
            viterbi_bptrs.at(i)[BERN_POLYI] = BERN_POLYI;
            if (a2a < i2a) {
                viterbi_bptrs.at(i)[BERN_POLYA] = BERN_POLYI;
            } else {
                viterbi_bptrs.at(i)[BERN_POLYA] = BERN_POLYA;
            }
        }

        // --- backwards viterbi pass:
        std::vector<BernoulliState> regions(bernoullis.size(), BERN_POLYI);
        std::vector<float> scores(bernoullis.size(), 0.0f);
        if (viterbi_scores.at(bernoullis.size()-1)[BERN_POLYI] < viterbi_scores.at(bernoullis.size()-1)[BERN_POLYA]) {
            regions[bernoullis.size()-1] = BERN_POLYA;
            scores[bernoullis.size()-1] = viterbi_scores.at(bernoullis.size()-1)[BERN_POLYA];
        } else {
            regions[bernoullis.size()-1] = BERN_POLYI;
            scores[bernoullis.size()-1] = viterbi_scores.at(bernoullis.size()-1)[BERN_POLYI];
        }
        for (size_t j=(bernoullis.size()-2); j > 0; --j) {
            regions[j] = viterbi_bptrs.at(j)[regions.at(j+1)];
            scores[j] = viterbi_scores.at(j)[regions.at(j+1)];
        }

        // --- format BernoulliOutputs structure and return:
        BernoulliOutputs output_vectors = { scores, regions };
        return output_vectors;
    }

    // Compute the Bernoulli HMM segmentation; this is the final public interface that gets called.
    BernoulliSegmentation segmentation(const SquiggleRead& sr, int start, int stop) const {
        // --- initialize BernoulliSegmentation (indices of -1 mean that the respective regions were not found):
        BernoulliSegmentation segmentation = { -1, -1 };
        // --- guard: if fewer than 100 samples in the region, return 'not found':
        if (stop - start < 100) {
            return segmentation;
        }

        // --- subset the squiggleread sequence, perform linear adjustment, and binarize via log-likelihood test:
        std::vector<float> squiggle(&sr.samples[start], &sr.samples[stop]);
        for (size_t i = 0; i < squiggle.size(); ++i) {
            squiggle.at(i) = (squiggle.at(i) - sr.scalings[0].shift) / sr.scalings[0].scale;
        }
        std::vector<int> bernoullis = this->log_lkhd_ratio_sequence(squiggle);

        // --- run viterbi algorithm:
        BernoulliOutputs viterbi_results = this->viterbi(bernoullis);

        // --- parse viterbi labels to find segmentation (keep indices at -1 if no region found):
        for (size_t i = 0; i < viterbi_results.labels.size(); ++i) {
            if (viterbi_results.labels.at(i) == BERN_POLYI) {
                segmentation.polyi = i;
            }
            if ((viterbi_results.labels.at(i) == BERN_POLYA) && (segmentation.polya < 0)) {
                segmentation.polya = i;
            }
        }
        return segmentation;
    }
};

// ================================================================================
// Estimate the duration profile for a single read.
//   Estimate the underlying read rate.
// * dpi_estimate_eventalign_duration_profile : compute median read rate via collapsed-
//     duration event-alignment.
// * dpi_estimate_unaligned_duration_profile : compute median read rate via collapsed-
//     durations, without event-alignment.
// ================================================================================
// Compute a read-rate based on event-alignment, collapsed by consecutive 5mer identity
// (N.B.: deprecated; using non-eventaligned durations seems to work just as well
// while being faster to run.)
double dpi_estimate_eventalign_duration_profile(SquiggleRead& sr,
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
double dpi_estimate_unaligned_duration_profile(const SquiggleRead& sr,
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
std::vector<double> dpi_fetch_event_durations(const SquiggleRead& sr,
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
// * dpi_estimate_polya_length : return an estimate of the read rate for this read.
// ================================================================================
// Compute an estimate of the number of nucleotides in the poly-A tail
double dpi_estimate_polya_length(const SquiggleRead& sr, const DPISegmentation& region_indices, const double read_rate)
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
// * dpi_pre_segmentation_qc: return a bool (true ~> FAIL) performing basic pre-seg QC.
// * dpi_post_segmentation_qc: check the segmentation results for failures.
// * dpi_post_estimation_qc: sanity check for estimates.
// ================================================================================
// QC before segmentation; check if event-alignment passes.
std::string dpi_pre_segmentation_qc(uint32_t suffix_clip, uint32_t prefix_clip, double transcript_length, const SquiggleRead& sr)
{
    std::string qc_tag;
    if (suffix_clip > 200) {
        // fail if this read has a long skip at end:
        qc_tag = "SUFFCLIP";
    } else {
        // pass if none of the above fail:
        qc_tag = "PASS";
    }
    return qc_tag;
}

// QC pass after constructing a segmentation; returns a QC flag represented as a string,
// either "PASS" or "NOREGION".
// These tests indicate that something went wrong in the segmentation algorithm.
std::string dpi_post_segmentation_qc(const DPISegmentation& region_indices, const SquiggleRead& sr)
{
    // fetch sizes of ADAPTER and POLYA regions:
    double num_adapter_samples = (region_indices.adapter+1) - region_indices.leader;
    double num_polya_samples = region_indices.polya - (region_indices.adapter+1);

    // check for NOREGION:
    std::string qc_tag;
    if (num_adapter_samples < 200.0 || num_polya_samples < 200.0) {
        qc_tag = "NOREGION";
    } else {
        qc_tag = "PASS";
    }
    return qc_tag;
}

// QC pass after performing estimation; returns a QC flag represented as a string.
// Currently returns either "PASS" or "ADAPTER".
// These tests indicate that something went wrong in the estimation algorithm.
std::string dpi_post_estimation_qc(const DPISegmentation& region_indices, const SquiggleRead& sr, double read_rate, double polya_length)
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

// QC pass to remove erroneous POLY{A,I} detection calls.
std::string post_boolhmm_detetection_qc(const BernoulliSegmentation& segmentation, int region_length) {
    // empirically-discovered threshold for number of samples needed for POLYA or POLYI region to be detected:
    double cutoff = 200;
    // detect regions:
    bool polyi_found = false;
    if (segmentation.polyi > cutoff) {
        polyi_found = true;
    }
    bool polya_found = false;
    if ((segmentation.polya > 0) && (region_length - segmentation.polya > cutoff)) {
        polya_found = true;
    }
    // compute output tag:
    std::string qc_tag;
    if (polyi_found && polya_found) {
        qc_tag = "A+I";
    } else if (!polyi_found && polya_found) {
        qc_tag = "POLYA-ONLY";
    } else if (polyi_found && !polya_found) {
        qc_tag = "POLYI-ONLY";
    } else {
        qc_tag = "NONE";
    }
    return qc_tag;
}

// ================================================================================
// Main Poly-A Code
//   Expose main functionality of this module via public-facing function.
// * dpi_estimate_polya_for_single_read : performs all data-fetching and joins all
//     of the above functions; writes results to file in TSV format
// * detect_polyi_main : wrap `dpi_estimate_polya_for_single_read` with OpenMP directives
//     for easy multi-threading across reads
// ================================================================================
// Write Poly(A) region segmentation and tail length estimation data to TSV
void dpi_estimate_polya_for_single_read(const ReadDB& read_db,
                                        const faidx_t* fai,
                                        FILE* out_fp,
                                        const bam_hdr_t* hdr,
                                        const bam1_t* record,
                                        size_t read_idx,
                                        int region_start,
                                        int region_end)
{
    //----- check if primary or secondary alignment by inspecting FLAG; skip read if secondary:
    if (record->core.flag & BAM_FSECONDARY) {
        return;
    }

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

    //----- construct SquiggleRead; if there are load issues, print -1's and skip compute:
    SquiggleRead sr(read_name, read_db, SRF_LOAD_RAW_SAMPLES);
    if (sr.fast5_path == "" || sr.events[0].empty()) {
        #pragma omp critical
        {
            fprintf(out_fp, "%s\t%s\t%zu\t-1.0\t-1.0\t-1.0\t-1.0\t-1.00\t-1.00\tREAD_FAILED_LOAD\n",
                read_name.c_str(), ref_name.c_str(), record->core.pos);
            if (opt::verbose == 1) {
                fprintf(out_fp,
                    "polya-samples\t%s\t%s\t-1\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\tREAD_FAILED_LOAD\n",
                    read_name.substr(0,6).c_str(), ref_name.c_str());
            }
            if (opt::verbose == 2) {
                fprintf(out_fp, "polya-durations\t%s\t-1\t-1.0\tREAD_FAILED_LOAD\n", read_name.substr(0,6).c_str());
            }
        }
        return;
    }

    //----- print clipping data if `verbose > 2` set:
    if (opt::verbose > 2) {
        fprintf(stderr, "[polya] read: %s length: %zu prefix clip: %zu suffix clip %zu\n",
                read_name.c_str(), sr.read_sequence.length(), prefix_clip, suffix_clip);
    }
    std::string sequenced_transcript = sr.read_sequence;

    //----- Perform pre-segmentation QC:
    std::string pre_segmentation_qc_flag = dpi_pre_segmentation_qc(suffix_clip, prefix_clip, sequenced_transcript.length(), sr);

    //----- perform HMM-based regional segmentation & post-segmentation QC:
    DPISegmentationHMM hmm(static_cast<float>(sr.scalings[0].scale),
                           static_cast<float>(sr.scalings[0].shift),
                           static_cast<float>(sr.scalings[0].var));
    DPISegmentation region_indices = hmm.segment_squiggle(sr);
    std::string post_segmentation_qc_flag = dpi_post_segmentation_qc(region_indices, sr);

    //----- compute duration profile for the read:
    double read_rate = dpi_estimate_unaligned_duration_profile(sr, fai, hdr, record, read_idx, strand_idx);

    //----- estimate number of nucleotides in poly-A tail & post-estimation QC:
    double polya_length = dpi_estimate_polya_length(sr, region_indices, read_rate);
    std::string post_estimation_qc_flag = dpi_post_estimation_qc(region_indices, sr, read_rate, polya_length);

    //----- Resolve QC flag based on priority:
    std::string qc_tag;
    if (post_segmentation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_segmentation_qc_flag;
    } else if (post_estimation_qc_flag.compare("PASS") != 0) {
        qc_tag = post_estimation_qc_flag;
    } else if (pre_segmentation_qc_flag.compare("PASS") != 0) {
        qc_tag = pre_segmentation_qc_flag;
    } else {
        qc_tag = "PASS";
    }

    //----- Detect POLY{A,I} region with BernoulliHMM:
    BernoulliHMM boolhmm(static_cast<float>(sr.scalings[0].scale),
                         static_cast<float>(sr.scalings[0].shift),
                         static_cast<float>(sr.scalings[0].var));
    BernoulliSegmentation poly_segmentation = boolhmm.segmentation(sr, region_indices.adapter+1, region_indices.polya);
    std::string poly_detect_tag = post_boolhmm_detetection_qc(poly_segmentation, (region_indices.polya-(region_indices.adapter+1)));

    //----- print annotations to TSV:
    double leader_sample_start = region_indices.start+1;
    double adapter_sample_start = region_indices.leader+1;
    double polya_sample_start = region_indices.adapter+1;
    double polya_sample_end = region_indices.polya;
    double transcr_sample_start = region_indices.polya+1;
    #pragma omp critical
    {
        fprintf(out_fp, "%s\t%s\t%zu\t%.1lf\t%.1lf\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%s\t%s\n",
                read_name.c_str(), ref_name.c_str(), record->core.pos,
                leader_sample_start, adapter_sample_start, polya_sample_start, transcr_sample_start,
                read_rate, polya_length, poly_detect_tag.c_str(), qc_tag.c_str());
        // if `verbose == 1`, print the samples (picoAmps) of the read,
        // up to the first 1000 samples of transcript region:
        if (opt::verbose == 1) {
            for (size_t i = 0; i < std::min(static_cast<size_t>(polya_sample_end)+1000, sr.samples.size()); ++i) {
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
                float s = sr.samples[i];
                float scaled_s = (s - sr.scalings[0].shift) / sr.scalings[0].scale;
                std::vector<float> s_probas = hmm.log_probas(s);
                fprintf(out_fp, "polya-samples\t%s\t%s\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",
                        read_name.substr(0,6).c_str(), ref_name.c_str(), i, s, scaled_s,
                        s_probas.at(0), s_probas.at(1), s_probas.at(2), s_probas.at(3), s_probas.at(4), s_probas.at(5),
                        tag.c_str());
            }
        }
        // if `verbose == 2`, print the raw event durations of the read:
        if (opt::verbose == 2) {
            std::vector<double> raw_durations = dpi_fetch_event_durations(sr, fai, hdr, record, read_idx, strand_idx);
            for (size_t i = 0; i < raw_durations.size(); ++i) {
                double dura = raw_durations[i];
                fprintf(out_fp, "polya-durations\t%s\t%zu\t%f\t%s\n",
                    read_name.substr(0,6).c_str(), i, dura, qc_tag.c_str());
            }
        }
    }
}

// Wrap poly-A estimation code for parallelism
int detect_polyi_main(int argc, char** argv)
{
    parse_detect_polyi_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // print header line:
    fprintf(stdout, "readname\tcontig\tposition\tleader_start\tadapter_start\tpolya_start\ttranscript_start\tread_rate\tpolya_length\tdetected\tqc_tag\n");

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(dpi_estimate_polya_for_single_read, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    // free allocated values:
    fai_destroy(fai);

    return EXIT_SUCCESS;
}
