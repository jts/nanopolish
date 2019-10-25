//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_methyltrain -- train a methylation model
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <map>
#include <omp.h>
#include <getopt.h>
#include <cstddef>
#include <random>
#include <functional>
#include <unordered_map>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "htslib/faidx.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_model_names.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_read_db.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_bedmethyl.h"
#include "training_core.hpp"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"
#include "logger.hpp"

#include "nanopolish_scorereads.h"
#include "Eigen/Dense"

using namespace std::placeholders;

extern float g_p_skip, g_p_skip_self, g_p_bad, g_p_bad_self;

//
// Enums
//
enum TrainingTarget
{
    TT_UNMETHYLATED_KMERS,
    TT_METHYLATED_KMERS,
    TT_ALL_KMERS
};

//
// Structs
//

struct StateSummary
{
    StateSummary() { num_matches = 0; num_skips = 0; num_stays = 0; }
    std::vector<StateTrainingData> events;

    int num_matches;
    int num_skips;
    int num_stays;
};

//
const Alphabet* train_alphabet_ptr = NULL;

//
// Typedefs
//

//
// Getopt
//
#define SUBPROGRAM "train"

static const char *TRAIN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *TRAIN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Train a methylation model\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"      --train-kmers=STR                train methylated, unmethylated or all kmers\n"
//"  -c  --calibrate                      recalibrate aligned reads to model before training\n"
"      --no-update-models               do not write out trained models\n"
//"      --output-scores                  optionally output read scores during training\n"
"  -r, --reads=FILE                     the ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -i, --input-model-filename           start training from this model\n"
"  -d, --output-directory               write results to this output directory\n"
"      --bedmethyl=FILE                 infer the presence/absence of methylation sites using FILE as a map\n"
"      --filter-policy=STR              filter data using parameters for STR data (R7, R9, R94)\n"
"      --rounds=NUM                     number of training rounds to perform\n"
"      --max-reads=NUM                  stop after processing NUM reads in each round\n"
"      --max-read-length=NUM            do not use reads exceeding NUM bp in length\n"
"      --progress                       print out a progress message\n"
//"      --stdv                           enable stdv modelling\n"
"      --max-events=NUM                 use NUM events for training (default: 1000)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static unsigned int calibrate=0;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;
    static std::string input_model_filename;
    static std::string output_directory;
    static std::string initial_model_type = "ONT";
    static std::string trained_model_type = "reftrained";
    static std::string bedmethyl_filename; 
    static TrainingTarget training_target = TT_METHYLATED_KMERS;
    static bool write_models = true;
    static bool output_scores = false;
    static unsigned progress = 0;
    static unsigned num_threads = 1;
    static unsigned batch_size = 128;
    static unsigned max_reads = -1;
    static unsigned max_read_length = -1;

    // Constants that determine which events to use for training
    static float min_event_duration = 0.002;
    static unsigned min_distance_from_alignment_end = 5;
    static unsigned min_number_of_events_to_train = 100;
    static unsigned num_training_rounds = 5;
    static unsigned max_events = 10000;
}

static const char* shortopts = "r:b:g:t:i:d:vnc";

enum { OPT_HELP = 1,
       OPT_VERSION,
       OPT_PROGRESS,
       OPT_NO_UPDATE_MODELS,
       OPT_TRAIN_KMERS,
       OPT_OUTPUT_SCORES,
       OPT_OUT_FOFN,
       OPT_STDV,
       OPT_LOG_LEVEL,
       OPT_FILTER_POLICY,
       OPT_NUM_ROUNDS,
       OPT_P_SKIP,
       OPT_P_SKIP_SELF,
       OPT_P_BAD,
       OPT_P_BAD_SELF,
       OPT_MAX_READS,
       OPT_MAX_READ_LENGTH,
       OPT_MAX_EVENTS,
       OPT_BEDMETHYL
     };

static const struct option longopts[] = {
    { "verbose",              no_argument,       NULL, 'v' },
    { "calibrate",            no_argument,       NULL, 'c' },
    { "reads",                required_argument, NULL, 'r' },
    { "bam",                  required_argument, NULL, 'b' },
    { "genome",               required_argument, NULL, 'g' },
    { "threads",              required_argument, NULL, 't' },
    { "input-model-filename", required_argument, NULL, 'i' },
    { "output-directory",     required_argument, NULL, 'd' },
    { "out-suffix",           required_argument, NULL, 's' },
    { "stdv",                 no_argument,       NULL, OPT_STDV },
    { "train-kmers",          required_argument, NULL, OPT_TRAIN_KMERS },
    { "p-skip",               required_argument, NULL, OPT_P_SKIP },
    { "p-skip-self",          required_argument, NULL, OPT_P_SKIP_SELF },
    { "p-bad",                required_argument, NULL, OPT_P_BAD },
    { "p-bad-self",           required_argument, NULL, OPT_P_BAD_SELF },
    { "output-scores",        no_argument,       NULL, OPT_OUTPUT_SCORES },
    { "no-update-models",     no_argument,       NULL, OPT_NO_UPDATE_MODELS },
    { "progress",             no_argument,       NULL, OPT_PROGRESS },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { "log-level",            required_argument, NULL, OPT_LOG_LEVEL },
    { "filter-policy",        required_argument, NULL, OPT_FILTER_POLICY },
    { "rounds",               required_argument, NULL, OPT_NUM_ROUNDS },
    { "max-reads",            required_argument, NULL, OPT_MAX_READS },
    { "max-read-length",      required_argument, NULL, OPT_MAX_READ_LENGTH },
    { "max-events",           required_argument, NULL, OPT_MAX_EVENTS },
    { "bedmethyl",            required_argument, NULL, OPT_BEDMETHYL },
    { NULL, 0, NULL, 0 }
};

void parse_train_options(int argc, char** argv)
{
    std::string training_target_str = "";
    std::string filter_policy_str = "";
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'i': arg >> opt::input_model_filename; break;
            case 'd': arg >> opt::output_directory; break;
            case 'v': opt::verbose++; break;
            case 'c': opt::calibrate = 1; break;
            case OPT_STDV: /*model_stdv() = true;*/ break;
            case OPT_NUM_ROUNDS: arg >> opt::num_training_rounds; break;
            case OPT_OUTPUT_SCORES: opt::output_scores = true; break;
            case OPT_TRAIN_KMERS: arg >> training_target_str; break;
            case OPT_FILTER_POLICY: arg >> filter_policy_str; break;
            case OPT_NO_UPDATE_MODELS: opt::write_models = false; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_P_SKIP: arg >> g_p_skip; break;
            case OPT_P_SKIP_SELF: arg >> g_p_skip_self; break;
            case OPT_P_BAD: arg >> g_p_bad; break;
            case OPT_P_BAD_SELF: arg >> g_p_bad_self; break;
            case OPT_MAX_READS: arg >> opt::max_reads; break;
            case OPT_MAX_READ_LENGTH: arg >> opt::max_read_length; break;
            case OPT_MAX_EVENTS: arg >> opt::max_events; break;
            case OPT_BEDMETHYL: arg >> opt::bedmethyl_filename; break;
            case OPT_HELP:
                std::cout << TRAIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TRAIN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                logger::Logger::set_level_from_option(arg.str());
                break;
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

    if(opt::output_directory.empty()) {
        std::cerr << SUBPROGRAM ": an --output-directory name must be provided\n";
        die = true;
    }

    if(opt::input_model_filename.empty()) {
        std::cerr << SUBPROGRAM ": an --input-model-filename must be provided\n";
        die = true;
    }

    // Parse the training target string
    if(training_target_str != "") {
        if(training_target_str == "unmethylated") {
            opt::training_target = TT_UNMETHYLATED_KMERS;
        } else if(training_target_str == "methylated") {
            opt::training_target = TT_METHYLATED_KMERS;
        } else if(training_target_str == "all") {
            opt::training_target = TT_ALL_KMERS;
        } else {
            std::cerr << SUBPROGRAM ": unknown --train-kmers string\n";
            die = true;
        }
    }

    // Parse the training target string
    if(filter_policy_str != "") {
        if(filter_policy_str == "R9") {
            opt::min_event_duration = 0.002;
        } else if(filter_policy_str == "R9.4") {
            opt::min_event_duration = 0.000;
        } else if(filter_policy_str == "R7") {
            opt::min_event_duration = 0.005;
        } else {
            std::cerr << SUBPROGRAM ": unknown --filter-policy\n";
            die = true;
        }
    }

    if (die) {
        std::cout << "\n" << TRAIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}


// recalculate shift, scale, drift, scale_sd from an alignment and the read
// returns true if the recalibration was performed
// in either case, sets residual to the L1 norm of the residual
extern bool recalibrate_model(SquiggleRead &sr,
                       const PoreModel& pore_model,
                       const int strand_idx,
                       const std::vector<EventAlignment> &alignment_output,
                       const bool scale_var,
                       const bool scale_drift);

// Use reservoir sampling to limit the amount of events for each kmer
void add_event_tmp(StateSummary& kmer_summary, StateTrainingData std, int event_count)
{
    // Add all events for each kmer until we have the maximum number of events to train
    if(event_count <= opt::max_events) {
        kmer_summary.events.push_back(std);
    } else {
        // Select a random number between 0 and the total number of events for this kmer
        std::mt19937 rng(event_count);
        std::uniform_int_distribution<int> gen(0, event_count);
        int rs_location = gen(rng);

        // Only update the training data if the random number is below the maximum number of events
        if(rs_location < opt::max_events) {
            kmer_summary.events[rs_location] = std;
        }
    }
}

std::string infer_read_sequence_from_variants_and_methylation(SquiggleRead& sr,
                                                              const PoreModel* pore_model,
                                                              const faidx_t* fai,
                                                              const bam_hdr_t* hdr,
                                                              const bam1_t* record,
                                                              const std::vector<BedMethylRecord>& methylation_records)
{
    int strand_idx = 0;
    uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

    // 
    std::string ref_name = hdr->target_name[record->core.tid];
    int ref_start_pos = record->core.pos;
    int ref_end_pos = bam_endpos(record);

    int reference_length;
    std::string reference_seq =
        get_reference_region_ts(fai, ref_name.c_str(), ref_start_pos, ref_end_pos, &reference_length);
    std::transform(reference_seq.begin(), reference_seq.end(), reference_seq.begin(), ::toupper);

    std::string read_outseq = reference_seq;
    
    // Search the variant collection for the index of the first/last variants to phase
    BedMethylRecord lower_search;
    lower_search.ref_name = ref_name;
    lower_search.ref_start_position = ref_start_pos;
    auto lower_iter = std::lower_bound(methylation_records.begin(), methylation_records.end(), lower_search, sortBedMethylByPosition);

    BedMethylRecord upper_search;
    upper_search.ref_name = ref_name;
    upper_search.ref_start_position = ref_end_pos;
    auto upper_iter = std::upper_bound(methylation_records.begin(), methylation_records.end(), upper_search, sortBedMethylByPosition);
    
    if(opt::verbose >= 1) {
        fprintf(stderr, "Inferring read %s %s:%u-%u %zu\n", 
            sr.read_name.c_str(), ref_name.c_str(), ref_start_pos, ref_end_pos, upper_iter - lower_iter);
    }

    SequenceAlignmentRecord seq_align_record(record);
    EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);
    Haplotype reference_haplotype(ref_name, ref_start_pos, reference_seq);

    //
    int min_flanking_sequence = 30;
    for(; lower_iter < upper_iter; ++lower_iter) {

        const BedMethylRecord& bmr = *lower_iter;
        if(bmr.strand == "-") {
            continue; // TODO: handle properly
        }

        // Convert BedMethyl to VCF as this is what Haplotype wants
        Variant v;
        v.ref_name = bmr.ref_name;
        v.ref_position = bmr.ref_start_position;
        v.ref_seq = reference_seq.substr(bmr.ref_start_position - ref_start_pos, 1);
        v.alt_seq = "M";

        // Calculate methylated/unmethylated priors
        double prior_methylated = bmr.coverage > 10 ? bmr.percent_methylated / 100.0 : 0.5;
        double log_prior_methylated = log(prior_methylated);
        double log_prior_unmethylated = log(1.0 - prior_methylated);

        // set up HMM data
        int calling_start = v.ref_position - min_flanking_sequence;
        int calling_end = v.ref_position + min_flanking_sequence;

        HMMInputData data;
        data.read = event_align_record.sr;
        data.strand = event_align_record.strand;
        data.rc = event_align_record.rc;
        data.event_stride = event_align_record.stride;
        data.pore_model = pore_model;

        int e1,e2;
        bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events,
                calling_start,
                calling_end,
                e1,
                e2);

        // The events of this read do not span the calling window, skip
        if(!bounded || fabs(e2 - e1) / (calling_start - calling_end) > MAX_EVENT_TO_BP_RATIO) {
            continue;
        }

        data.event_start_idx = e1;
        data.event_stop_idx = e2;

        Haplotype calling_haplotype =
            reference_haplotype.substr_by_reference(calling_start, calling_end);

        HMMInputSequence ref_hmm_input(calling_haplotype.get_sequence(), pore_model->pmalphabet);
        double ref_score = profile_hmm_score(ref_hmm_input, data, hmm_flags) + log_prior_unmethylated;

        bool good_haplotype = calling_haplotype.apply_variant(v);
        if(good_haplotype) {
            HMMInputSequence alt_hmm_input(calling_haplotype.get_sequence(), pore_model->pmalphabet);
            double alt_score = profile_hmm_score(alt_hmm_input, data, hmm_flags) + log_prior_methylated;
            double log_sum = add_logs(alt_score, ref_score);
            double log_p_ref = ref_score - log_sum;
            double log_p_alt = alt_score - log_sum;
            char call;
            double log_p_wrong;

            if(alt_score > ref_score) {
                call = v.alt_seq[0];
                log_p_wrong = log_p_ref;
            } else {
                call = v.ref_seq[0];
                log_p_wrong = log_p_alt;
            }

            double q_score = -10 * log_p_wrong / log(10);
            char q_char = (int)q_score + 33;

            /*
            fprintf(stderr, "\t%s prior meth: %.2lf score: %.2lf %.2lf %c p_wrong: %.3lf Q: %d QC: %c\n", 
                v.key().c_str(), prior_methylated, ref_score, alt_score, call, log_p_wrong, (int)q_score, q_char);
            */

            int out_position = v.ref_position - ref_start_pos;
            if(read_outseq[out_position] != v.ref_seq[0]) {
                fprintf(stderr, "warning: reference base at position %d does not match variant record (%c != %c)\n",
                        v.ref_position, v.ref_seq[0], read_outseq[out_position]);
            }

            read_outseq[out_position] = call;
        }
    }


    return read_outseq;
}

// Update the training data with aligned events from a read
void add_aligned_events_for_read(const ReadDB& read_db,
                                 const faidx_t* fai,
                                 std::vector<StateSummary>& training_data,
                                 const std::string& training_kit,
                                 const std::string& training_alphabet,
                                 size_t training_k,
                                 std::unordered_map<uint32_t, int>& event_count,
                                 const std::vector<BedMethylRecord>& methylation_records,
                                 const bam_hdr_t* hdr,
                                 const bam1_t* record,
                                 size_t read_idx,
                                 int region_start,
                                 int region_end)
{

    // only support training template strand
    size_t strand_idx = 0;

    // optionally skip long reads
    if(opt::max_read_length > 0 && record->core.l_qseq > opt::max_read_length) {
        return;
    }

    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);

    // load read
    SquiggleRead sr(read_name, read_db);

    // skip if 1D reads and this is the wrong strand
    if(!sr.has_events_for_strand(strand_idx)) {
        return;
    }

    assert(training_kit == sr.get_model_kit_name(strand_idx));
    assert(training_k == sr.get_model_k(strand_idx));

    // set k
    uint32_t k = sr.get_model_k(strand_idx);

    // Get the reference sequence for this read
    std::string ref_name = hdr->target_name[record->core.tid];
    int alignment_start_pos = record->core.pos;
    int alignment_end_pos = bam_endpos(record);

    // skip if region too short
    // TODO: make parameter
    if(alignment_end_pos - alignment_start_pos < 1000) {
        return;
    }
    
    // get the pore model we want for the type of training data
    const PoreModel* pore_model = sr.get_model(strand_idx, train_alphabet_ptr->get_name());

    // infer the truth sequence for this read using the input  genomic variants and methylation sites
    std::string reference_seq = infer_read_sequence_from_variants_and_methylation(sr, pore_model, fai, hdr, record, methylation_records);

    // Reverse complement if necessary
    if(bam_is_rev(record)) {
        reference_seq = train_alphabet_ptr->reverse_complement(reference_seq);
    }

    // align events to the reference sequence


    // prepare data structures for the event-guided forward/backward to calculate state posteriors
    AdaBandedParameters alignment_parameters;
    alignment_parameters.bandwidth = 25;
    alignment_parameters.p_skip = 0.005;
    //alignment_parameters.verbose = 1;
    alignment_parameters.min_average_log_emission = -3.5;
    alignment_parameters.max_stay_threshold = 200;
    alignment_parameters.min_posterior = 1e-1;
    Haplotype reference_haplotype(ref_name, alignment_start_pos, reference_seq);
    SequenceAlignmentRecord seq_align_record(record);
    EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);
    
    // calculate event->kmer posteriors
    std::vector<EventKmerPosterior> state_assignments = 
        guide_banded_simple_posterior(sr, *pore_model, reference_haplotype, event_align_record, alignment_parameters);

    size_t edge_ignore = 50;
    if(state_assignments.size() < 2*edge_ignore) {
        return;
    }

    for(size_t i = edge_ignore; i < state_assignments.size() - edge_ignore; ++i) {
        const auto& ekp = state_assignments[i];
        
        size_t kmer_idx = ekp.kmer_idx;
        size_t event_idx = ekp.event_idx;
        double log_posterior = ekp.log_posterior;

        /*
        size_t kmer_idx = ekp.ref_pos;
        size_t event_idx = ekp.read_pos;
        double log_posterior = 0.0;
        */
        // Get the rank of the kmer that we aligned to (on the sequencing strand, = model_kmer)
        uint32_t rank = train_alphabet_ptr->kmer_rank(reference_seq.c_str() + kmer_idx, k);
        assert(rank < training_data.size());
        auto& kmer_summary = training_data[rank];

        // We only use this event for training if its not at the end of the alignment
        // (to avoid bad alignments around the read edges) and if its not too short (to
        // avoid bad measurements from effecting the levels too much)
        bool use_for_training = sr.get_duration( event_idx, strand_idx) >= opt::min_event_duration &&
                                sr.get_fully_scaled_level(event_idx, strand_idx) >= 1.0;

        if(use_for_training) {
            StateTrainingData std(sr.get_fully_scaled_level(event_idx, strand_idx),
                                  sr.get_scaled_stdv(event_idx, strand_idx),
                                  sr.scalings[strand_idx].var,
                                  sr.scalings[strand_idx].scale,
                                  log_posterior);

            #pragma omp critical(kmer)
            {
                event_count[rank]++;
                add_event_tmp(kmer_summary, std, event_count[rank]);
            }
        }
    }
}

struct TrainingResult
{
    PoreModel trained_model;
    size_t num_kmers_to_train;
    size_t num_kmers_trained;
};


TrainingResult train_model_from_events(const PoreModel& current_model,
                                       const std::vector<StateSummary>& summaries,
                                       std::ofstream& training_ofs,
                                       FILE* summary_fp)
{
    TrainingResult result;

    // Initialize the new model from the current model
    result.trained_model = current_model;
    result.num_kmers_trained = 0;
    result.num_kmers_to_train = 0;

    uint32_t k = result.trained_model.k;
    std::string model_key = PoreModelSet::get_model_key(result.trained_model);
    std::string model_short_name = result.trained_model.metadata.get_short_name();

    // Generate the complete set of kmers
    std::string gen_kmer(k, 'A');
    std::vector<std::string> all_kmers;
    for(size_t ki = 0; ki < summaries.size(); ++ki) {
        all_kmers.push_back(gen_kmer);
        train_alphabet_ptr->lexicographic_next(gen_kmer);
    }
    assert(gen_kmer == std::string(k, 'A'));
    assert(all_kmers.front() == std::string(k, 'A'));
    assert(all_kmers.back() == std::string(k, 'T'));

    Progress training_timer("");

    // Update means for each kmer
    omp_set_num_threads(opt::num_threads);

    #pragma omp parallel for
    for(size_t ki = 0; ki < summaries.size(); ++ki) {
        assert(ki < all_kmers.size());
        std::string kmer = all_kmers[ki];

        /*
        // write the observed values to a tsv file
        #pragma omp critical
        {
            for(size_t ei = 0; ei < summaries[ki].events.size(); ++ei) {
                summaries[ki].events[ei].write_tsv(training_ofs, model_short_name, kmer);
            }
        }
        */

        bool is_m_kmer = kmer.find('M') != std::string::npos;
        bool update_kmer = opt::training_target == TT_ALL_KMERS ||
                           (is_m_kmer && opt::training_target == TT_METHYLATED_KMERS) ||
                           (!is_m_kmer && opt::training_target == TT_UNMETHYLATED_KMERS);

        #pragma omp atomic
        result.num_kmers_to_train += update_kmer;

        bool trained = false;

        // only train if there are a sufficient number of events for this kmer
        if(update_kmer && summaries[ki].events.size() >= opt::min_number_of_events_to_train) {

            // train a mixture model where a minority of k-mers aren't methylated
            ParamMixture mixture;

            float incomplete_methylation_rate = 0.05f;
            std::string um_kmer = train_alphabet_ptr->unmethylate(kmer);
            size_t um_ki = train_alphabet_ptr->kmer_rank(um_kmer.c_str(), k);

            // Initialize the training parameters. If this is a kmer containing
            // a methylation site we train a two component mixture, otherwise
            // just fit a gaussian
            float major_weight = is_m_kmer ? 1 - incomplete_methylation_rate : 1.0f;
            mixture.log_weights.push_back(log(major_weight));
            mixture.params.push_back(current_model.get_parameters(ki));
            
            /*
            if(is_m_kmer) {
                // add second unmethylated component
                mixture.log_weights.push_back(std::log(incomplete_methylation_rate));
                mixture.params.push_back(current_model.get_parameters(um_ki));
            }
            */

            if(opt::verbose > 1) {
                fprintf(stderr, "INIT__MIX %s\t%s\t[%.2lf %.2lf %.2lf]\t[%.2lf %.2lf %.2lf]\n", model_key.c_str(), kmer.c_str(),
                    std::exp(mixture.log_weights[0]), mixture.params[0].level_mean, mixture.params[0].level_stdv,
                    std::exp(mixture.log_weights[1]), mixture.params[1].level_mean, mixture.params[1].level_stdv);
            }

            ParamMixture trained_mixture = train_gaussian_mixture(summaries[ki].events, mixture);

            if(opt::verbose > 1) {
                fprintf(stderr, "TRAIN_MIX %s\t%s\t[%.2lf %.2lf %.2lf]\t[%.2lf %.2lf %.2lf]\n", model_key.c_str(), kmer.c_str(),
                    std::exp(trained_mixture.log_weights[0]), trained_mixture.params[0].level_mean, trained_mixture.params[0].level_stdv,
                    std::exp(trained_mixture.log_weights[1]), trained_mixture.params[1].level_mean, trained_mixture.params[1].level_stdv);
            }

            #pragma omp critical
            {
                result.trained_model.states[ki] = trained_mixture.params[0];
                result.num_kmers_trained += 1;
                trained = true;
            }

#if 0
            if (false && model_stdv()) {
                ParamMixture ig_mixture;
                // weights
                ig_mixture.log_weights = trained_mixture.log_weights;
                // states
                ig_mixture.params.emplace_back(trained_mixture.params[0]);

                if(is_m_kmer) {
                    ig_mixture.params.emplace_back(current_model.get_parameters(um_ki));
                }
                // run training
                auto trained_ig_mixture = train_invgaussian_mixture(summaries[ki].events, ig_mixture);

                LOG("methyltrain", debug)
                    << "IG_INIT__MIX " << model_key << " " << kmer.c_str() << " ["
                    << std::fixed << std::setprecision(5) << ig_mixture.params[0].sd_mean << " "
                    << ig_mixture.params[1].sd_mean << "]" << std::endl
                    << "IG_TRAIN_MIX " << model_key << " " << kmer.c_str() << " ["
                    << trained_ig_mixture.params[0].sd_mean << " "
                    << trained_ig_mixture.params[1].sd_mean << "]" << std::endl;

                // update state
                #pragma omp critical
                {
                    result.trained_model.states[ki] = trained_ig_mixture.params[0];
                }
            }
#endif
        }

        #pragma omp critical
        {
            fprintf(summary_fp, "%s\t%s\t%d\t%d\t%d\t%zu\t%d\t%.2lf\t%.2lf\n",
                                    model_short_name.c_str(), kmer.c_str(),
                                    summaries[ki].num_matches, summaries[ki].num_skips, summaries[ki].num_stays,
                                    summaries[ki].events.size(), trained, result.trained_model.states[ki].level_mean, result.trained_model.states[ki].level_stdv);
        }
    }

    fprintf(stderr, "[train] trained %zu/%zu k-mers in %.2lfs\n", result.num_kmers_trained, result.num_kmers_to_train, training_timer.get_elapsed_seconds());

    return result;
}

PoreModel train_round(const ReadDB& read_db,
                      const std::string& kit_name,
                      const std::string& alphabet,
                      size_t k,
                      size_t round,
                      const PoreModel& current_model,
                      const std::vector<BedMethylRecord>& methylation_records)
{
    // initialize summary data
    std::vector<StateSummary> model_training_data(current_model.get_num_states());
    std::unordered_map<uint32_t, int> event_count;

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // Open the BAM and iterate over reads
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.set_min_mapping_quality(20);

    auto f = std::bind(add_aligned_events_for_read,
                       std::ref(read_db),
                       std::ref(fai),
                       std::ref(model_training_data),
                       std::ref(kit_name),
                       std::ref(alphabet),
                       k,
                       std::ref(event_count),
                       std::ref(methylation_records),
                       _1, _2, _3, _4, _5); // parameters BamProcessor passes in

    processor.parallel_run(f);

    // Set up output files
    std::string model_base_name = PoreModelSet::get_model_key(current_model);

    // open the summary file
    std::stringstream summary_fn;
    summary_fn << opt::output_directory << "/" << model_base_name << ".round" << round << ".summmary.tsv";
    FILE* summary_fp = fopen(summary_fn.str().c_str(), "w");
    fprintf(summary_fp, "model_short_name\tkmer\tnum_matches\tnum_skips\t"
                         "num_stays\tnum_events_for_training\twas_trained\t"
                         "trained_level_mean\ttrained_level_stdv\n");

    // open the tsv file with the raw training data
    std::stringstream events_fn;
    events_fn << opt::output_directory << "/" << model_base_name << ".round" << round << ".events.tsv";
    std::ofstream training_ofs(events_fn.str());

    // write out a header for the training data
    StateTrainingData::write_header(training_ofs);

    //
    TrainingResult result = train_model_from_events(current_model, model_training_data, training_ofs, summary_fp);

    // add the updated model into the collection (or replace what is already there)
    PoreModelSet::add_model(result.trained_model);

    // write the updated model to disk
    if(opt::write_models && result.num_kmers_trained > 0) {
        std::stringstream model_fn;
        model_fn << opt::output_directory << "/" << model_base_name << ".round" << round << ".model";
        result.trained_model.write(model_fn.str(), model_base_name + ".model");
    }

    // cleanup
    fai_destroy(fai);
    fclose(summary_fp);

    return result.trained_model;
}

extern void write_models(const std::map<std::string, PoreModel>& models, int round);

int train_main(int argc, char** argv)
{
    parse_train_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // Read the input model and initialize the pore model set to contain it
    PoreModel current_model(opt::input_model_filename);
    PoreModelSet::add_model(current_model);

    // Create the output directory, if it doesn't exist
    struct stat st = {0};
    if(stat(opt::output_directory.c_str(), &st) == -1) {
            mkdir(opt::output_directory.c_str(), 0700);
    }

    // optionally load in methylation calls for this sample
    // these calls will be used as a map to determine read-specific
    // methylation patterns for training
    std::vector<BedMethylRecord> bedmethyl_records;
    if(!opt::bedmethyl_filename.empty()) {
        if(opt::region.empty()) {
            fprintf(stderr, "error: a training region must be provided when using the bedmethyl option\n");
            exit(EXIT_FAILURE);
        }

        std::string contig;
        int start_base;
        int end_base;
        parse_region_string(opt::region, contig, start_base, end_base);

        bedmethyl_records = read_bedmethyl_for_region(opt::bedmethyl_filename, contig, start_base, end_base);
        fprintf(stderr, "[train] read %zu methylation records in region\n", bedmethyl_records.size());
        
        // Sort variants by reference coordinate
        std::sort(bedmethyl_records.begin(), bedmethyl_records.end(), sortBedMethylByPosition);
    }

    //
    std::string training_kit = current_model.metadata.get_kit_name();
    train_alphabet_ptr = current_model.pmalphabet;
    size_t training_k = current_model.k;
    fprintf(stderr, "[train] initialized %s for alphabet %s for %zu-mers\n", training_kit.c_str(), train_alphabet_ptr->get_name(), training_k);

    for(size_t round = 0; round < opt::num_training_rounds; round++) {

        fprintf(stderr, "[train] round %zu\n", round);
        current_model = train_round(read_db, training_kit, train_alphabet_ptr->get_name(), training_k, round, current_model, bedmethyl_records);
    }

    return EXIT_SUCCESS;
}

