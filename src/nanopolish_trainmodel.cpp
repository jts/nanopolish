//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_trainmodel - train a new pore model from
// the FAST5 output of a basecaller
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

//
// Getopt
//
#define SUBPROGRAM "trainmodel"

static const char *TRAINMODEL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2016 Ontario Institute for Cancer Research\n";

static const char *TRAINMODEL_USAGE_MESSAGE =
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

void parse_trainmodel_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << TRAINMODEL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TRAINMODEL_VERSION_MESSAGE;
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
        std::cout << "\n" << TRAINMODEL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::fofn_file = argv[optind++];
}

int trainmodel_main(int argc, char** argv)
{
    parse_trainmodel_options(argc, argv);

    std::ifstream fofn_reader(opt::fofn_file);
    std::string fast5_name;

    // Read input
    std::vector<SquiggleRead*> reads;
    while(getline(fofn_reader, fast5_name)) {
        fprintf(stderr, "Loading %s\n", fast5_name.c_str());
        reads.push_back(new SquiggleRead(fast5_name, fast5_name));
    }
    fprintf(stderr, "Loaded %zu reads\n", reads.size());

    // 
    unsigned int basecalled_k = 5; // TODO: infer this
    size_t num_kmers = gDNAAlphabet.get_num_strings(basecalled_k);
    unsigned int training_strand = T_IDX; // template training for now

    typedef std::vector<StateTrainingData> TrainingDataVector;
    typedef std::vector<TrainingDataVector> KmerTrainingData;

    // This vector is indexed by read, then kmer, then event
    std::vector<KmerTrainingData> read_training_data;

    size_t read_idx = 0;
    for(auto* read : reads) {
        
        // Initialize a vector-of-vectors to hold training data for this read
        read_training_data.push_back(KmerTrainingData());
        KmerTrainingData& kmer_training_data = read_training_data.back();

        // Initialize the events-by-kmer vector 
        kmer_training_data.resize(num_kmers);
        printf("KTD: %zu\n", kmer_training_data.size());

        const std::string& read_sequence = read->read_sequence;
        size_t n_kmers = read_sequence.size() - basecalled_k + 1;
        for(size_t ki = 0; ki < n_kmers; ++ki) {

            IndexPair event_range_for_kmer = read->base_to_event_map[ki].indices[training_strand];
            
            // skip kmers without events and with multiple events
            if(event_range_for_kmer.start == -1 || 
               event_range_for_kmer.start != event_range_for_kmer.stop) {
                continue;
            }

            std::string kmer = read_sequence.substr(ki, basecalled_k);
            size_t kmer_rank = gDNAAlphabet.kmer_rank(kmer.c_str(), basecalled_k);
            assert(kmer_rank < num_kmers);

            size_t event_idx = event_range_for_kmer.start;
            double level = read->events[training_strand][event_idx].mean;
            double stdv = read->events[training_strand][event_idx].stdv;

            printf("read: %zu kmer: %s ki: %zu lvl: %.2lf dur: %.5lf\n", read_idx, kmer.c_str(), kmer_rank, level, read->events[training_strand][event_idx].duration);
            StateTrainingData std(level, stdv, read->pore_model[training_strand].var);
            kmer_training_data[kmer_rank].push_back(std);
        }

        read_idx++;
        /*
        std::string kmer(basecalled_k, 'A');
        for(size_t ki = 0; ki < num_kmers; ki++) {
            size_t num_events = kmer_training_data[ki].size();
            printf("%s:", kmer.c_str());
            for(size_t ei = 0; ei < num_events; ++ei) {
                printf("%.2lf\t", kmer_training_data[ki][ei].level_mean);
            }
            printf("\n");

            gDNAAlphabet.lexicographic_next(kmer);
        }
        */
    }

    // Select the read with the most events as the "baseline" read for generating the model
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

    // Set the initial pore model based on the read with the most events
    PoreModel pore_model(basecalled_k);
    pore_model.states.resize(num_kmers);
    pore_model.scaled_states.resize(num_kmers);
    pore_model.scaled_params.resize(num_kmers);

    pore_model.shift = 0.0;
    pore_model.scale = 1.0;
    pore_model.drift = 0.0;
    pore_model.var = 1.0;
    pore_model.scale_sd = 1.0;
    pore_model.var_sd = 1.0;
    
    auto& kmer_training_data_for_selected = read_training_data[max_events_index];

    std::vector<bool> use_kmer(num_kmers, false);
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

            // mark this kmer as valid
            use_kmer[ki] = true;
            pore_model.states[ki].level_mean = median;
            pore_model.states[ki].level_stdv = 1.0;
            printf("k: %zu median: %.2lf values: %s\n", ki, median, ss.str().c_str());
        }
    }
    pore_model.bake_gaussian_parameters();

    // Apply model to each read
    for(auto* read: reads) {
        read->pore_model[training_strand] = pore_model;
    }

    // Recalibrate read
    for(auto* read: reads) {

        // We generate a vector of event-to-kmer mapping to use the recalibration linear solver
        std::vector<EventAlignment> alignment;
        const std::string& read_sequence = read->read_sequence;
        size_t n_kmers = read_sequence.size() - basecalled_k + 1;
        for(size_t ki = 0; ki < n_kmers; ++ki) {

            IndexPair event_range_for_kmer = read->base_to_event_map[ki].indices[training_strand];
            
            // skip kmers without events and with multiple events
            if(event_range_for_kmer.start == -1 || 
               event_range_for_kmer.start != event_range_for_kmer.stop) {
                continue;
            }

            std::string kmer = read_sequence.substr(ki, basecalled_k);
            size_t kmer_rank = gDNAAlphabet.kmer_rank(kmer.c_str(), basecalled_k);
            assert(kmer_rank < num_kmers);

            size_t event_idx = event_range_for_kmer.start;
            
            // Only use this kmer if it is part of the initial model
            if(use_kmer[kmer_rank]) {
                EventAlignment ea;
                // ref data
                ea.ref_name = ""; // not needed
                ea.ref_kmer = kmer;
                ea.ref_position = ki;
                ea.read_idx = -1; // not needed
                ea.strand_idx = training_strand;
                ea.event_idx = event_idx;
                ea.rc = false;
                ea.model_kmer = kmer;
                ea.hmm_state = 'M'; // recalibration code only uses "M" alignments
                alignment.push_back(ea);
            }
        }

        recalibrate_model(*read, 
                          training_strand,
                          alignment,
                          &gDNAAlphabet,
                          false);
        
        const PoreModel& read_model = read->pore_model[training_strand];
        printf("[recalibration] events: %zu alignment: %zu shift: %.2lf scale: %.2lf drift: %.4lf var: %.2lf\n", 
            read->events[training_strand].size(),
            alignment.size(),
            read_model.shift, 
            read_model.scale, 
            read_model.drift, 
            read_model.var);

    }
    
    // Deallocate input reads
    for(auto* read : reads) {
        delete read;
    }

    return 0;
}
