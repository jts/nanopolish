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
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_khmm_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

//
// Structs
//
struct StateSummary
{
    StateSummary() : n(0), mean_sum(0.f), var_sum(0.f) {}
    size_t n;
    float mean_sum;
    float var_sum;
};

//
// Typedefs
//
typedef std::map<std::string, std::vector<PoreModelStateParams>> ModelMap;
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
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --progress                       print out a progress message\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string models_fofn;
    static std::string region;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:w:m:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "reads",            required_argument, NULL, 'r' },
    { "bam",              required_argument, NULL, 'b' },
    { "genome",           required_argument, NULL, 'g' },
    { "window",           required_argument, NULL, 'w' },
    { "threads",          required_argument, NULL, 't' },
    { "models-fofn",      required_argument, NULL, 'm' },
    { "progress",         no_argument,       NULL, OPT_PROGRESS },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Realign the read in event space
void train_read(const ModelMap& model_map,
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
        std::string curr_model = sr.model_name[strand_idx];
        auto model_iter = model_map.find(curr_model);

        if(model_iter != model_map.end()) {
            sr.replace_pore_model(strand_idx, model_iter->second);
        } else {
            assert(false && "Model not found");
        }

        // Align to the new model
        std::vector<EventAlignment> alignment_output = 
            align_read_to_ref(sr, fai, hdr, record, read_idx, strand_idx, region_start, region_end);

        // Update model observations
        #pragma omp critical
        {
            // Get the training data for this model
            auto& emission_map = training[curr_model];

            for(size_t i = 0; i < alignment_output.size(); ++i) {
                const EventAlignment& ea = alignment_output[i];
                std::string model_kmer = ea.rc ? reverse_complement(ea.ref_kmer) : ea.ref_kmer;
                uint32_t rank = kmer_rank(model_kmer.c_str(), K);
                auto& kmer_summary = emission_map[rank];

                // Here we re-scale the event to the model
                float event_mean = sr.get_drift_corrected_level(ea.event_idx, strand_idx);

                // unshifted mean
                float expected_mean = (sr.pore_model[strand_idx].states[rank].level_mean - 
                                       sr.pore_model[strand_idx].shift) / sr.pore_model[strand_idx].scale;

                sr.get_drift_corrected_level(ea.event_idx, strand_idx);
                event_mean = (event_mean - sr.pore_model[strand_idx].shift) / sr.pore_model[strand_idx].scale;
                printf("%zu\t%s\t%s\t%zu\t%zu\t%zu\t%.2lf\t%c\n", round, curr_model.c_str(), model_kmer.c_str(), ea.ref_position, read_idx, ea.event_idx, event_mean, ea.hmm_state);
                kmer_summary.n += 1;
                kmer_summary.mean_sum += event_mean;
                kmer_summary.var_sum += pow(event_mean - expected_mean, 2.0);
            }
        }
    } // for strands
}

ModelMap read_models_fofn(const std::string& fofn_name)
{
    ModelMap out;
    std::ifstream fofn_reader(fofn_name);
    std::string model_filename;
    while(getline(fofn_reader, model_filename)) {
        printf("reading %s\n", model_filename.c_str());

        std::ifstream model_reader(model_filename);
        std::string model_line;

        std::string model_name;
        std::vector<PoreModelStateParams> states;

        while(getline(model_reader, model_line)) {
            std::stringstream parser(model_line);

            // Extract the model name from the header
            if(model_line.find("#model_file") != std::string::npos) {
                std::string dummy;
                parser >> dummy >> model_name;
            }
            
            // skip the rest of the header
            if(model_line[0] == '#' || model_line.find("kmer") == 0) {
                continue;
            }
            
            std::string kmer;
            PoreModelStateParams params;
            parser >> kmer >> params.level_mean >> params.level_stdv >> params.sd_mean >> params.sd_stdv;
            states.push_back(params);
        }
            
        assert(!model_name.empty());
        out[model_name] = states;
    }
    return out;
}

void parse_methyltrain_options(int argc, char** argv)
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
            case 'm': arg >> opt::models_fofn; break;
            case 'v': opt::verbose++; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << METHYLTRAIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << METHYLTRAIN_VERSION_MESSAGE;
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
    
    if(opt::models_fofn.empty()) {
        std::cerr << SUBPROGRAM ": a --models-fofn file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << METHYLTRAIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}


ModelMap train_one_round(const ModelMap& models, const Fast5Map& name_map, size_t round)
{
    // Initialize the training summary stats for each kmer for each model
    ModelTrainingMap model_training_data;
    for(auto model_iter = models.begin(); model_iter != models.end(); model_iter++) {
        std::vector<StateSummary> summaries(model_iter->second.size()); // one per kmer in the model
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
                    train_read(models, name_map, fai, hdr, record, read_idx, clip_start, clip_end, round, model_training_data);
                }
            }

            num_reads_realigned += num_records_buffered;
            num_records_buffered = 0;

            if(num_reads_realigned > 1000) break;
        }

        if(opt::progress) {
            fprintf(stderr, "Realigned %zu reads in %.1lfs\r", num_reads_realigned, progress.get_elapsed_seconds());
        }
    } while(result >= 0);
    
    assert(num_records_buffered == 0);
 
    // Process the training results
    ModelMap trained_models;

    for(auto model_training_iter = model_training_data.begin(); 
             model_training_iter != model_training_data.end(); model_training_iter++) {
        
        // Initialize trained model from input model
        auto model_iter = models.find(model_training_iter->first);
        assert(model_iter != models.end());
        trained_models[model_training_iter->first] = model_iter->second;
        std::vector<PoreModelStateParams>& new_pm = trained_models[model_training_iter->first];

        // Update means for each kmer
        std::string kmer = "AAAAA";
        const std::vector<StateSummary>& summaries = model_training_iter->second;
        for(size_t ki = 0; ki < summaries.size(); ++ki) {

            if(kmer.find("CG") != std::string::npos) {
                float mu_prime = summaries[ki].mean_sum / summaries[ki].n;
                float var_prime = summaries[ki].var_sum / summaries[ki].n;
                new_pm[ki].level_mean = mu_prime;
                new_pm[ki].level_stdv = sqrt(var_prime);
                fprintf(stderr, "%s %s %.2lf %.2lf\n", model_training_iter->first.c_str(), kmer.c_str(), new_pm[ki].level_mean, new_pm[ki].level_stdv);
            }
            lexicographic_next(kmer);
        }
    }


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

int methyltrain_main(int argc, char** argv)
{
    parse_methyltrain_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    ModelMap models = read_models_fofn(opt::models_fofn);
    
    for(size_t round = 0; round < 10; round++) {
        fprintf(stderr, "Starting round %zu\n", round);
        ModelMap trained_models = train_one_round(models, name_map, round);
        // Write the model
        for(auto model_iter = trained_models.begin(); 
                 model_iter != trained_models.end(); model_iter++) {
        
            std::stringstream outname;
            outname << model_iter->first << ".trained.round" << round + 1;
            std::ofstream writer(outname.str());
            
            const std::vector<PoreModelStateParams>& states = model_iter->second;

            std::string curr_kmer = "AAAAA";
            for(size_t ki = 0; ki < states.size(); ++ki) {
                writer << curr_kmer << "\t" << states[ki].level_mean << "\t" << states[ki].level_stdv << "\n";
                lexicographic_next(curr_kmer);
            }
            writer.close();
        }
        models = trained_models;
    }
    return EXIT_SUCCESS;
}

