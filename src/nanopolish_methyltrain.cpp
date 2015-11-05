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
#include "htslib/faidx.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "training_core.hpp"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"
#include "logger.hpp"
#include "pfor.hpp"

#include "nanopolish_scorereads.h"
#include "../eigen/Eigen/Dense"

//
// Enus
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
const Alphabet* mtrain_alphabet = NULL;

//
// Typedefs
//
typedef std::map<std::string, std::vector<StateSummary>> ModelTrainingMap;

//
// Getopt
//
#define SUBPROGRAM "methyltrain"

static const char *METHYLTRAIN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *METHYLTRAIN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Train a methylation model\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -m, --models-fofn=FILE               read the models to be trained from the FOFN\n"
"      --train-kmers=STR                train methylated, unmethylated or all kmers\n"
"  -c  --calibrate                      recalibrate aligned reads to model before training\n"
"      --no-update-models               do not write out trained models\n"
"      --output-scores                  optionally output read scores during training\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -s, --out-suffix=STR                 name output files like <strand>.out_suffix\n"
"      --out-fofn=FILE                  write the names of the output models into FILE\n"
"      --progress                       print out a progress message\n"
"      --stdv                           enable stdv modelling\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static unsigned int calibrate=0;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string models_fofn;
    static std::string region;
    static std::string out_suffix = ".trained";
    static std::string out_fofn = "trained.fofn";
    static TrainingTarget training_target = TT_METHYLATED_KMERS;
    static bool write_models = true;
    static bool output_scores = false;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:m:vnc";

enum { OPT_HELP = 1, 
       OPT_VERSION, 
       OPT_PROGRESS, 
       OPT_NO_UPDATE_MODELS, 
       OPT_TRAIN_KMERS, 
       OPT_OUTPUT_SCORES,
       OPT_OUT_FOFN,
       OPT_STDV,
       OPT_LOG_LEVEL
     };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "calibrate",          no_argument,       NULL, 'c' },
    { "reads",              required_argument, NULL, 'r' },
    { "bam",                required_argument, NULL, 'b' },
    { "genome",             required_argument, NULL, 'g' },
    { "window",             required_argument, NULL, 'w' },
    { "threads",            required_argument, NULL, 't' },
    { "models-fofn",        required_argument, NULL, 'm' },
    { "out-suffix",         required_argument, NULL, 's' },
    { "stdv",               no_argument,       NULL, OPT_STDV },
    { "out-fofn",           required_argument, NULL, OPT_OUT_FOFN },
    { "train-kmers",        required_argument, NULL, OPT_TRAIN_KMERS },
    { "output-scores",      no_argument,       NULL, OPT_OUTPUT_SCORES },
    { "no-update-models",   no_argument,       NULL, OPT_NO_UPDATE_MODELS },
    { "progress",           no_argument,       NULL, OPT_PROGRESS },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "log-level",          required_argument, NULL, OPT_LOG_LEVEL },
    { NULL, 0, NULL, 0 }
};

std::string get_model_short_name(const std::string& model_name)
{
    static std::map< std::string, std::string > model_name_map = {
        { "r7.3_template_median68pA.model", "t.005" },
        { "r7.3_complement_median68pA_pop1.model", "c.p1.005" },
        { "r7.3_complement_median68pA_pop2.model", "c.p2.005" },
        { "r7.3_e6_70bps_6mer_template_median68pA.model", "t.006" },
        { "r7.3_e6_70bps_6mer_complement_median68pA_pop1.model", "c.p1.006" },
        { "r7.3_e6_70bps_6mer_complement_median68pA_pop2.model", "c.p2.006" }
    };
    auto iter = model_name_map.find(model_name);
    if(iter == model_name_map.end()) {
        fprintf(stderr, "Error: unknown model %s\n", model_name.c_str());
        exit(EXIT_FAILURE);
    }
    return iter->second;
}

// recalculate shift, scale, drift, scale_sd from an alignment and the read
void recalibrate_model(SquiggleRead &sr,
                       const int strand_idx,
                       const std::vector<EventAlignment> &alignment_output, 
                       bool scale_var) 
{
    std::vector<double> raw_events, times, level_means, level_stdvs;
    uint32_t k = sr.pore_model[strand_idx].k;

    //std::cout << "Previous pore model parameters: " << sr.pore_model[strand_idx].shift << ", " 
    //                                                << sr.pore_model[strand_idx].scale << ", " 
    //                                                << sr.pore_model[strand_idx].drift << ", " 
    //                                                << sr.pore_model[strand_idx].var   << std::endl;

    // extract necessary vectors from the read and the pore model; note do not want scaled values
    for ( const auto &ea : alignment_output ) {
        if(ea.hmm_state == 'M') {
            std::string model_kmer = ea.rc ? mtrain_alphabet->reverse_complement(ea.ref_kmer) : ea.ref_kmer;
            uint32_t rank = mtrain_alphabet->kmer_rank(model_kmer.c_str(), k);

            raw_events.push_back ( sr.get_uncorrected_level(ea.event_idx, strand_idx) );
            times.push_back      ( sr.get_time(ea.event_idx, strand_idx) );
            level_means.push_back( sr.pore_model[strand_idx].states[rank].level_mean );
            level_stdvs.push_back( sr.pore_model[strand_idx].states[rank].level_stdv );
        }
    }

    const int minNumEventsToRescale = 500;
    if (raw_events.size() < minNumEventsToRescale) 
        return;

    // Assemble linear system corresponding to weighted least squares problem
    // Can just directly call a weighted least squares solver, but there's enough
    // structure in our problem it's a little faster just to build the normal eqn
    // matrices ourselves
    Eigen::Matrix3d A;
    Eigen::Vector3d b;

    for (int i=0; i<3; i++) {
        b(i) = 0.;
        for (int j=0; j<3; j++)
            A(i,j) = 0.;
    }

    for (size_t i=0; i<raw_events.size(); i++) {
        double inv_var = 1./(level_stdvs[i]*level_stdvs[i]);
        double mu = level_means[i];
        double t  = times[i];
        double e  = raw_events[i];

        A(0,0) += inv_var;  A(0,1) += mu*inv_var;    A(0,2) += t*inv_var;
                            A(1,1) += mu*mu*inv_var; A(1,2) += mu*t*inv_var;
                                                     A(2,2) += t*t*inv_var;

        b(0) += e*inv_var;
        b(1) += mu*e*inv_var;
        b(2) += t*e*inv_var;
    }
    A(1,0) = A(0,1);
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);

    // perform the linear solve
    Eigen::Vector3d x = A.fullPivLu().solve(b);

    double shift = x(0);
    double scale = x(1);
    double drift = x(2);

    sr.pore_model[strand_idx].shift = shift;
    sr.pore_model[strand_idx].scale = scale;
    sr.pore_model[strand_idx].drift = drift;

    if (scale_var) {
        double var = 0.;
        for (size_t i=0; i<raw_events.size(); i++) {
            double yi = (raw_events[i] - shift - scale*level_means[i] - drift*times[i]);
            var+= yi*yi/(level_stdvs[i]*level_stdvs[i]);
        }
        var /= raw_events.size();

        sr.pore_model[strand_idx].var   = var;
    }

    if (sr.pore_model[strand_idx].is_scaled)
        sr.pore_model[strand_idx].bake_gaussian_parameters();

    //std::cout << "Updated pore model parameters:  " << sr.pore_model[strand_idx].shift << ", " 
    //                                                << sr.pore_model[strand_idx].scale << ", " 
    //                                                << sr.pore_model[strand_idx].drift << ", " 
    //                                                << sr.pore_model[strand_idx].var   << std::endl;
}

// Update the training data with aligned events from a read
void add_aligned_events(const ModelMap& model_map,
                        const Fast5Map& name_map,
                        const faidx_t* fai,
                        const bam_hdr_t* hdr,
                        const bam1_t* record,
                        size_t read_idx,
                        int region_start,
                        int region_end,
                        size_t round,
                        ModelTrainingMap& training)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    std::string fast5_path = name_map.get_path(read_name);

    // load read
    SquiggleRead sr(read_name, fast5_path);

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        
        // replace model with the training model
        std::string curr_model = sr.pore_model[strand_idx].name;
        auto model_iter = model_map.find(curr_model);

        if(model_iter != model_map.end()) {
            sr.pore_model[strand_idx].update_states(model_iter->second);
        } else {
            printf("Error: model %s not found\n", curr_model.c_str());
            assert(false && "Model not found");
        }
        
        // set k
        uint32_t k = sr.pore_model[strand_idx].k;

        // Align to the new model
        EventAlignmentParameters params;
        params.sr = &sr;
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
        params.strand_idx = strand_idx;
 
        params.alphabet = mtrain_alphabet;
        params.read_idx = read_idx;
        params.region_start = region_start;
        params.region_end = region_end;
        std::vector<EventAlignment> alignment_output = align_read_to_ref(params);
        if (alignment_output.size() == 0)
            return;

        // Update pore model based on alignment
        double orig_score;
        if (opt::output_scores) {
            orig_score = model_score(sr, strand_idx, fai, alignment_output, 500);
            #pragma omp critical(print)
            std::cout << round << " " << curr_model << " " << read_idx << " " << strand_idx << " Original " << orig_score << std::endl;
        }

        if ( opt::calibrate ) {
            recalibrate_model(sr, strand_idx, alignment_output, false);

            if (opt::output_scores) {
                double rescaled_score = model_score(sr, strand_idx, fai, alignment_output, 500);
                #pragma omp critical(print)
                {
                    std::cout << round << " " << curr_model << " " << read_idx << " " << strand_idx << " Rescaled " << rescaled_score << std::endl;
                    std::cout << round << " " << curr_model << " " << read_idx << " " << strand_idx << " Delta " << rescaled_score-orig_score << std::endl;
                }
            }
        }
        // Update model observations
        //emit_event_alignment_tsv(stdout, sr, params, alignment_output);

        // Get the training data for this model
        auto& emission_map = training[curr_model];

        for(size_t i = 0; i < alignment_output.size(); ++i) {
            const EventAlignment& ea = alignment_output[i];
            std::string model_kmer = ea.model_kmer;

            // Grab the previous/next model kmer
            // If the read is from the same strand as the reference
            // the next kmer comes from the next alignment_output (and vice-versa)
            // other the indices are swapped
            int next_stride = ea.rc ? -1 : 1;

            std::string prev_kmer = "";
            std::string next_kmer = "";

            if(i > 0 && i < alignment_output.size() - 1) {
                assert(alignment_output[i + next_stride].event_idx - ea.event_idx == 1);
                assert(alignment_output[i - next_stride].event_idx - ea.event_idx == -1);

                // check for exactly one base of movement along the reference
                if( std::abs(alignment_output[i + next_stride].ref_position - ea.ref_position) == 1) {
                    next_kmer = alignment_output[i + next_stride].model_kmer;
                }

                if( std::abs(alignment_output[i - next_stride].ref_position - ea.ref_position) == 1) {
                    prev_kmer = alignment_output[i - next_stride].model_kmer;
                }
            }

            uint32_t rank = mtrain_alphabet->kmer_rank(model_kmer.c_str(), k);
            assert(rank < emission_map.size());
            auto& kmer_summary = emission_map[rank];

            // Should we use this event for training?
            bool use_for_training = i > 5 && 
                i + 5 < alignment_output.size() &&
                alignment_output[i].hmm_state == 'M';

            if(use_for_training) {
                StateTrainingData std(sr, ea, rank, prev_kmer, next_kmer);
                #pragma omp critical(kmer)
                kmer_summary.events.push_back(std);
            }

            if(ea.hmm_state == 'M')  {
                #pragma omp atomic
                kmer_summary.num_matches += 1;
            } else if(ea.hmm_state == 'E') {
                #pragma omp atomic
                kmer_summary.num_stays += 1;
            }
        }
    } // for strands
}

void parse_methyltrain_options(int argc, char** argv)
{
    std::string training_target_str = "";
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'm': arg >> opt::models_fofn; break;
            case 's': arg >> opt::out_suffix; break;
            case 'v': opt::verbose++; break;
            case 'c': opt::calibrate = 1; break;
            case OPT_STDV: model_stdv() = true; break;
            case OPT_OUT_FOFN: arg >> opt::out_fofn; break;
            case OPT_OUTPUT_SCORES: opt::output_scores = true; break;
            case OPT_TRAIN_KMERS: arg >> training_target_str; break;
            case OPT_NO_UPDATE_MODELS: opt::write_models = false; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << METHYLTRAIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << METHYLTRAIN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                Logger::set_level_from_option(arg.str());
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
    
    if(opt::models_fofn.empty()) {
        std::cerr << SUBPROGRAM ": a --models-fofn file must be provided\n";
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

    if (die) {
        std::cout << "\n" << METHYLTRAIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

struct Training_PFor_Structs
{
    struct Input
    {
        std::string kmer;
        size_t ki;
    };
    struct Chunk_Output
    {
        ~Chunk_Output() {
            if (training_sink_osp()) (*training_sink_osp()) << training_os.str();
            if (summary_sink_osp()) (*summary_sink_osp()) << summary_os.str();
        }
        std::ostringstream training_os;
        std::ostringstream summary_os;
        static std::ostream*& training_sink_osp() {
            static std::ostream* _training_sink_osp = nullptr;
            return _training_sink_osp;
        }
        static std::ostream*& summary_sink_osp() {
            static std::ostream* _summary_sink_osp = nullptr;
            return _summary_sink_osp;
        }
    };
};

ModelMap train_one_round(const ModelMap& models, const Fast5Map& name_map, size_t round)
{
    // Initialize the training summary stats for each kmer for each model
    ModelTrainingMap model_training_data;
    for(auto model_iter = models.begin(); model_iter != models.end(); model_iter++) {
        std::vector<StateSummary> summaries(model_iter->second.get_num_states()); // one per kmer in the model
        model_training_data[model_iter->first] = summaries;
    }

    // Open the BAM and iterate over reads

    // load bam file
    htsFile* bam_fh = sam_open(opt::bam_file.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = opt::bam_file + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    hts_itr_t* itr;

    // If processing a region of the genome, only emit events aligned to this window
    int clip_start = -1;
    int clip_end = -1;

    if(opt::region.empty()) {
        // TODO: is this valid?
        itr = sam_itr_queryi(bam_idx, HTS_IDX_START, 0, 0);
    } else {

        fprintf(stderr, "Region: %s\n", opt::region.c_str());
        itr = sam_itr_querys(bam_idx, hdr, opt::region.c_str());
        hts_parse_reg(opt::region.c_str(), &clip_start, &clip_end);
    }

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    // Initialize iteration
    std::vector<bam1_t*> records(opt::batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_realigned = 0;
    size_t num_records_buffered = 0;
    Progress progress("[methyltrain]");

    do {
        assert(num_records_buffered < records.size());
        
        // read a record into the next slot in the buffer
        result = sam_itr_next(bam_fh, itr, records[num_records_buffered]);
        num_records_buffered += result >= 0;

        // realign if we've hit the max buffer size or reached the end of file
        if(num_records_buffered == records.size() || result < 0) {
            #pragma omp parallel for            
            for(size_t i = 0; i < num_records_buffered; ++i) {
                bam1_t* record = records[i];
                size_t read_idx = num_reads_realigned + i;
                if( (record->core.flag & BAM_FUNMAP) == 0) {
                    add_aligned_events(models, name_map, fai, hdr, record, read_idx, clip_start, clip_end, round, model_training_data);
                }
            }

            num_reads_realigned += num_records_buffered;
            num_records_buffered = 0;
        }

        if(opt::progress) {
            fprintf(stderr, "Realigned %zu reads in %.1lfs\r", num_reads_realigned, progress.get_elapsed_seconds());
        }
    } while(result >= 0);
    
    assert(num_records_buffered == 0);
    progress.end();

    std::stringstream training_fn;
    training_fn << opt::bam_file << ".round" << round << ".methyltrain.tsv";
    
    std::stringstream summary_fn;
    summary_fn << opt::bam_file << ".methyltrain.summary";

    std::ofstream summary_ofs(summary_fn.str());
    std::ofstream training_ofs(training_fn.str());

    // training header
    StateTrainingData::write_header(training_ofs);

    // Process the training results
    ModelMap trained_models;
    
    for(auto model_training_iter = model_training_data.begin(); 
             model_training_iter != model_training_data.end(); model_training_iter++) {
        
        // Initialize trained model from input model
        auto model_iter = models.find(model_training_iter->first);
        assert(model_iter != models.end());

        std::string model_name = model_training_iter->first;
        std::string model_short_name = get_model_short_name(model_name);

        trained_models[model_training_iter->first] = model_iter->second;
        PoreModel& new_pm = trained_models[model_training_iter->first];

        const std::vector<StateSummary>& summaries = model_training_iter->second;

        // Update means for each kmer
        uint32_t k = new_pm.k;
        //std::string kmer(k, 'A');
        //for(size_t ki = 0; ki < summaries.size(); ++ki, mtrain_alphabet->lexicographic_next(kmer)) {

        //
        // parallel for
        //
        Training_PFor_Structs::Input crt_input;
        crt_input.ki = 0;
        crt_input.kmer = std::string(k, 'A');
        Training_PFor_Structs::Chunk_Output::training_sink_osp() = &training_ofs;
        Training_PFor_Structs::Chunk_Output::summary_sink_osp() = &summary_ofs;
        pfor< Training_PFor_Structs >(
             nullptr, // do_before
             [&] (Training_PFor_Structs::Input& in) {
                 if (crt_input.ki >= summaries.size())
                 {
                     return false;
                 }
                 else
                 {
                     in = crt_input;
                     ++crt_input.ki;
                     mtrain_alphabet->lexicographic_next(crt_input.kmer);
                     return true;
                 }
             },
             [&] (Training_PFor_Structs::Input& in, Training_PFor_Structs::Chunk_Output& out)
             {
                 size_t& ki = in.ki;
                 std::string kmer = in.kmer;

            // write a training file
            for(size_t ei = 0; ei < summaries[ki].events.size(); ++ei) {
                summaries[ki].events[ei].write_tsv(out.training_os, model_short_name, kmer);
            }

            // write to the summary file
            out.summary_os << model_short_name << '\t'
                           << kmer << '\t'
                           << summaries[ki].num_matches << '\t'
                           << summaries[ki].num_skips << '\t'
                           << summaries[ki].num_stays << std::endl;

            bool is_m_kmer = kmer.find('M') != std::string::npos;
            bool update_kmer = opt::training_target == TT_ALL_KMERS ||
                               (is_m_kmer && opt::training_target == TT_METHYLATED_KMERS) ||
                               (!is_m_kmer && opt::training_target == TT_UNMETHYLATED_KMERS);
            if (not update_kmer or summaries[ki].events.size() < 100) {
                //continue;
                return;
            }

            ParamMixture mixture;

            // train a mixture model where a minority of k-mers aren't methylated
            
            // unmethylated component
            float um_rate = 0.05f;
            std::string um_kmer = mtrain_alphabet->unmethylate(kmer);
            size_t um_ki = mtrain_alphabet->kmer_rank(um_kmer.c_str(), k);

            mixture.log_weights.push_back(std::log(um_rate));
            mixture.params.push_back(model_iter->second.get_parameters(um_ki));

            mixture.log_weights.push_back(std::log(1 - um_rate));
            mixture.params.push_back(model_iter->second.get_parameters(ki));

            LOG("methyltrain", debug)
                << "INIT__MIX " << model_training_iter->first << '\t' << kmer << "\t["
                << std::fixed << std::setprecision(2) << std::exp(mixture.log_weights[0]) << ' '
                << mixture.params[0].level_mean << ' '
                << mixture.params[0].level_stdv << "]\t["
                << std::exp(mixture.log_weights[1]) << ' '
                << mixture.params[1].level_mean << ' '
                << mixture.params[1].level_stdv << "]" << std::endl;

            ParamMixture trained_mixture = train_gaussian_mixture(summaries[ki].events, mixture);

            LOG("methyltrain", debug)
                << "TRAIN_MIX " << model_training_iter->first << '\t' << kmer << "\t["
                << std::fixed << std::setprecision(2) << std::exp(trained_mixture.log_weights[0]) << ' '
                << trained_mixture.params[0].level_mean << ' '
                << trained_mixture.params[0].level_stdv << "]\t["
                << std::exp(trained_mixture.log_weights[1]) << ' '
                << trained_mixture.params[1].level_mean << ' '
                << trained_mixture.params[1].level_stdv << "]" << std::endl;

            new_pm.states[ki] = trained_mixture.params[1];

            if (model_stdv()) {
                ParamMixture ig_mixture;
                // weights
                ig_mixture.log_weights = trained_mixture.log_weights;
                // states
                ig_mixture.params.emplace_back(model_iter->second.get_parameters(um_ki));
                ig_mixture.params.emplace_back(trained_mixture.params[1]);
                // run training
                auto trained_ig_mixture = train_invgaussian_mixture(summaries[ki].events, ig_mixture);
                LOG("methyltrain", debug)
                    << "IG_INIT__MIX " << model_training_iter->first.c_str() << " " << kmer.c_str() << " ["
                    << std::fixed << std::setprecision(5) << ig_mixture.params[0].sd_mean << " "
                    << ig_mixture.params[1].sd_mean << "]" << std::endl
                    << "IG_TRAIN_MIX " << model_training_iter->first.c_str() << " " << kmer.c_str() << " ["
                    << trained_ig_mixture.params[0].sd_mean << " "
                    << trained_ig_mixture.params[1].sd_mean << "]" << std::endl;
                // update state
                new_pm.states[ki] = trained_ig_mixture.params[1];
                new_pm.states[ki].update_sd_stdv();
                new_pm.states[ki].update_logs();
            } // if model_stdv()
             }, // process_item
             nullptr, // do_after
             opt::num_threads,
             1); // pfor
    } // model_training_iter

    // cleanup records
    for(size_t i = 0; i < records.size(); ++i) {
        bam_destroy1(records[i]);
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    return trained_models;
}

void write_models(ModelMap& models)
{
    // file-of-filenames containing the new models
    std::ofstream fofn_writer(opt::out_fofn);

    // Write the model
    for(auto model_iter = models.begin(); 
             model_iter != models.end(); model_iter++) {

        assert(!model_iter->second.model_filename.empty());
        std::string outname   =  get_model_short_name(model_iter->second.name) + opt::out_suffix;
        std::string modelname =  model_iter->first + opt::out_suffix;
        models[model_iter->first].write( outname, modelname );

        fofn_writer << outname << "\n";
    }

}

int methyltrain_main(int argc, char** argv)
{
    parse_methyltrain_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    ModelMap models = read_models_fofn(opt::models_fofn);
    
    // Set the alphabet for this run to be the alphabet for the first model
    assert(!models.empty());
    mtrain_alphabet = models.begin()->second.pmalphabet;

    const size_t TRAINING_ROUNDS = 5;
    for(size_t round = 0; round < TRAINING_ROUNDS; round++) {
        fprintf(stderr, "Starting round %zu\n", round);
        ModelMap trained_models = train_one_round(models, name_map, round);
        if(opt::write_models) {
            write_models(trained_models);
        }
        models = trained_models;
    }
    return EXIT_SUCCESS;
}

