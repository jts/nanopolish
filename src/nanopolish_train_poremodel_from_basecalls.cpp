//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_train_poremodel_from_basecalls - train a 
// new pore model from the FAST5 output of a basecaller
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
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_methyltrain.h"
#include "training_core.hpp"
#include "profiler.h"
#include "logger.hpp"

//
// Typedefs
//
typedef std::vector<StateTrainingData> TrainingDataVector;
typedef std::vector<TrainingDataVector> KmerTrainingData;

//
// Getopt
//
#define SUBPROGRAM "train-poremodel-from-basecalls"

static const char *TRAIN_POREMODEL_FROM_BASECALLS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2016 Ontario Institute for Cancer Research\n";

static const char *TRAIN_POREMODEL_FROM_BASECALLS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] input.fofn\n"
"Train a new pore model using the basecalled reads in input.fofn\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string fofn_file;
}

static const char* shortopts = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_train_poremodel_from_basecalls_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << TRAIN_POREMODEL_FROM_BASECALLS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TRAIN_POREMODEL_FROM_BASECALLS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1) {
        std::cerr << SUBPROGRAM ": not enough arguments\n";
        die = true;
    }

    if (argc - optind > 1) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << TRAIN_POREMODEL_FROM_BASECALLS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::fofn_file = argv[optind++];
}

PoreModel initialize_pore_model(const std::vector<KmerTrainingData>& read_training_data,
                                const size_t k)
{
    size_t num_kmers_in_alphabet = gDNAAlphabet.get_num_strings(k);
    
    // Select the read with the most events to form the basis for the pore model distribution
    size_t max_events = 0;
    size_t max_events_index = 0;
    for(size_t rti = 0; rti < read_training_data.size(); ++rti) {
        auto& kmer_training_data = read_training_data[rti];
        size_t total_events = 0;

        for(size_t ki = 0; ki < kmer_training_data.size(); ++ki) {
            total_events += kmer_training_data[ki].size();
        }
        printf("read %zu has %zu events (max: %zu, %zu)\n", rti, total_events, max_events, max_events_index);

        if(total_events > max_events) {
            max_events = total_events;
            max_events_index = rti;
        }
    }

    // Set the initial pore model
    PoreModel pore_model(k);
    pore_model.states.resize(num_kmers_in_alphabet);

    auto& kmer_training_data_for_selected = read_training_data[max_events_index];

    for(size_t ki = 0; ki < kmer_training_data_for_selected.size(); ++ki) {
        std::vector<double> values;
        std::stringstream ss;
        for(size_t ei = 0; ei < kmer_training_data_for_selected[ki].size(); ++ei) {
            values.push_back(kmer_training_data_for_selected[ki][ei].level_mean);
            ss << values.back() << " ";
        }

        // Set the kmer's mean parameter to be the median of the recorded values
        std::sort(values.begin(), values.end());

        size_t n = values.size();
        double median;
        if(n == 0) {
            median = 0.0f;
        } else {
            if(n % 2 == 0) {
                median = (values[n / 2 - 1] + values[n/2]) / 2.0f;
            } else {
                median = values[n/2];
            }

            pore_model.states[ki].level_mean = median;
            pore_model.states[ki].level_stdv = 1.0;
            pore_model.states[ki].sd_mean = 0.0;
            pore_model.states[ki].sd_stdv = 0.0;
            pore_model.states[ki].sd_lambda = 0.0;
            pore_model.states[ki].update_logs();

            printf("k: %zu median: %.2lf values: %s\n", ki, median, ss.str().c_str());
        }
    }

    return pore_model;
}

void alignment_to_training_data(const SquiggleRead* read,
                                const std::vector<EventAlignment>& alignment,
                                const size_t k,
                                size_t read_idx,
                                KmerTrainingData* out_data,
                                FILE* tsv_writer)
{
    for(auto const& a : alignment) {
        size_t kmer_rank = gDNAAlphabet.kmer_rank(a.model_kmer.c_str(), k);
        assert(kmer_rank < out_data->size());
        assert(a.strand_idx == 0);
        assert(a.event_idx < read->events[a.strand_idx].size());

        double level = read->get_fully_scaled_level(a.event_idx, a.strand_idx);
        double stdv = read->events[a.strand_idx][a.event_idx].stdv;

        // If the scale/shift values are off, or the events are erroneous, the scaled events can have negative values
        // causing the training to implode. Filter these here.
        if(level >= 1.0) {
            StateTrainingData std(level, stdv, read->scalings[a.strand_idx].var, read->scalings[a.strand_idx].scale);
            out_data->at(kmer_rank).push_back(std);
        }

        if(tsv_writer) {
            fprintf(tsv_writer, "%zu\t%s\t%.2lf\t%.5lf\n", read_idx, a.model_kmer.c_str(), level, read->events[a.strand_idx][a.event_idx].duration);
        }
    }
}


int train_poremodel_from_basecalls_main(int argc, char** argv)
{
    parse_train_poremodel_from_basecalls_options(argc, argv);
#if 0
    std::ifstream fofn_reader(opt::fofn_file);
    std::string fast5_name;
    
    //Logger::set_level_from_option("debug1");

    // parameters 
    unsigned int basecalled_k = 5; // TODO: infer this
    size_t num_kmers_in_alphabet = gDNAAlphabet.get_num_strings(basecalled_k);
    unsigned int training_strand = T_IDX; // template training for now

    // Read input
    std::vector<SquiggleRead*> reads;
    while(getline(fofn_reader, fast5_name)) {
        fprintf(stderr, "Loading %s\n", fast5_name.c_str());
        SquiggleRead* read = new SquiggleRead(fast5_name, fast5_name, SRF_NO_MODEL);
     
        // initialize the scaling parameters to defaults
        PoreModel& read_pore_model = read->pore_model[training_strand];
        read_pore_model.shift = 0.0;
        read_pore_model.scale = 1.0;
        read_pore_model.drift = 0.0;
        read_pore_model.var = 1.0;
        read_pore_model.scale_sd = 1.0;
        read_pore_model.var_sd = 1.0;

        reads.push_back(read);
    }
    fprintf(stderr, "Loaded %zu reads\n", reads.size());

    // This vector is indexed by read, then kmer, then event
    std::vector<KmerTrainingData> read_training_data;

    FILE* tsv_writer = fopen("train_poremodel_from_basecalls.tsv", "w");
    fprintf(tsv_writer, "read_idx\tkmer\tlevel_mean\tduration\n");

    size_t read_idx = 0;
    for(auto* read : reads) {
        
        assert(false && "This is deprecated until it handles label_shift");
        // extract alignment of events to k-mers
        assert(read->read_type != SRT_2D);
        std::vector<EventAlignment> alignment = 
            read->get_eventalignment_for_1d_basecalls(read->read_sequence,
                                                      read->base_to_event_map,
                                                      basecalled_k,
                                                      training_strand,
                                                      0);

        // convert the alignment into model training data for this read
        KmerTrainingData ktd(num_kmers_in_alphabet);
        alignment_to_training_data(read, 
                                   alignment,
                                   basecalled_k,
                                   read_idx,
                                   &ktd,
                                   tsv_writer);

        read_training_data.push_back(ktd);
        read_idx++;
    }
    
    // Select the read with the most events as the "baseline" read for generating the model
    PoreModel initial_pore_model = initialize_pore_model(read_training_data, basecalled_k);
    PoreModel current_pore_model = initial_pore_model;

    for(size_t iteration = 0; iteration < 10; iteration++) {

        // Determine the k-mers that have been trained
        std::vector<bool> trained_kmers(num_kmers_in_alphabet, false);

        size_t num_trained = 0;
        for(size_t kmer_idx = 0; kmer_idx < num_kmers_in_alphabet; ++kmer_idx) {

            // untrained kmers have a mean of 0.0
            trained_kmers[kmer_idx] = current_pore_model.states[kmer_idx].level_mean > 1.0;
            num_trained += trained_kmers[kmer_idx];
        }

        // Recalibrate the scaling parameters for each read and collect new training data
        KmerTrainingData kmer_training_data(num_kmers_in_alphabet);
        for(size_t read_idx = 0; read_idx < reads.size(); ++read_idx) {

            SquiggleRead* read = reads[read_idx];
            
            // Apply new model to the read
            read->pore_model[training_strand] = current_pore_model;

            // generate alignment
            assert(false && "This is deprecated until it handles label_shift");
            assert(read->read_type != SRT_2D);
            std::vector<EventAlignment> alignment = 
                read->get_eventalignment_for_1d_basecalls(read->read_sequence,
                                                          read->base_to_event_map,
                                                          basecalled_k,
                                                          training_strand,
                                                          0);

            // filter the alignment to only contain k-mers that have a distribution
            std::vector<EventAlignment> filtered_alignment;
            for(size_t i = 0; i < alignment.size(); ++i) {
                size_t kmer_rank = gDNAAlphabet.kmer_rank(alignment[i].model_kmer.c_str(), basecalled_k);
                if(trained_kmers[kmer_rank]) {
                    filtered_alignment.push_back(alignment[i]);
                }
            }

            // recalibrate shift/scale/etc using the filtered alignment
            recalibrate_model(*read, 
                              training_strand,
                              filtered_alignment,
                              &gDNAAlphabet,
                              false, true);
        
            const PoreModel& read_model = read->pore_model[training_strand];
            printf("[recalibration] read %zu events: %zu alignment: %zu shift: %.2lf scale: %.2lf drift: %.4lf var: %.2lf\n", 
                read_idx,
                read->events[training_strand].size(),
                alignment.size(),
                read_model.shift, 
                read_model.scale, 
                read_model.drift, 
                read_model.var);

            // skip reads that aren't behaving well TODO: fix
            if(read_model.scale < 0.9 || read_model.scale > 1.1) {
                continue;
            }

            // collect kmer training data from this read
            alignment_to_training_data(read, 
                                       alignment,
                                       basecalled_k,
                                       read_idx,
                                       &kmer_training_data,
                                       NULL);
        }

        // Write the training data as a tsv file
        std::ofstream training_data_tsv("training_data.tsv");
        StateTrainingData::write_header(training_data_tsv);
        std::string model_kmer(basecalled_k, 'A');
        for(size_t kmer_idx = 0; kmer_idx < num_kmers_in_alphabet; kmer_idx++) {

            for(size_t di = 0; di < kmer_training_data[kmer_idx].size(); di++) {
                kmer_training_data[kmer_idx][di].write_tsv(training_data_tsv, "template.5mer", model_kmer);
            }
            gDNAAlphabet.lexicographic_next(model_kmer);
        }

        // Train new gaussians for each k-mer
        model_kmer = std::string(basecalled_k, 'A');
        PoreModel new_pore_model = current_pore_model;
        for(size_t kmer_idx = 0; kmer_idx < num_kmers_in_alphabet; kmer_idx++) {

            // we use the gaussian mixture machinery but only fit one component in the case
            ParamMixture input_mixture;
            fprintf(stderr, "training %s with %zu events\n", model_kmer.c_str(), kmer_training_data[kmer_idx].size());
            // This is intentially broad and doesn't matter in the one-component case
            PoreModelStateParams initial_params;
            initial_params.level_mean = 200;
            initial_params.level_stdv = 50;
            initial_params.update_logs();

            input_mixture.log_weights.push_back(log(1.0));
            input_mixture.params.push_back(initial_params);
               
            ParamMixture trained_mixture = train_gaussian_mixture(kmer_training_data[kmer_idx], input_mixture);
            new_pore_model.states[kmer_idx] = trained_mixture.params[0];
            new_pore_model.states[kmer_idx].level_stdv = 1.5;
            gDNAAlphabet.lexicographic_next(model_kmer);
        }
        new_pore_model.bake_gaussian_parameters();
        current_pore_model = new_pore_model;
    }
    current_pore_model.write("r9.template.5mer.base.model", "r9.template.5mer.base.model");

    // Deallocate input reads
    for(auto* read : reads) {
        delete read;
    }
#endif
    return 0;
}
