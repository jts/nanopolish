//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_dump_initial_alignment.cpp -- write out
// the event-to-basecall alignment.
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
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_read_db.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

// Tell KSEQ what functions to use to open/read files
KSEQ_INIT(gzFile, gzread)

//
// Getopt
//
#define SUBPROGRAM "dump-initial-alignment"

static const char *DUMPALIGNMENT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *DUMPALIGNMENT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa\n"
"Align nanopore events to reference k-mers\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static int num_threads = 1;
}

static const char* shortopts = "r:t:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS };

static const struct option longopts[] = {
    { "verbose",             no_argument,       NULL, 'v' },
    { "reads",               required_argument, NULL, 'r' },
    { "threads",             required_argument, NULL, 't' },
    { "help",                no_argument,       NULL, OPT_HELP },
    { "version",             no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_dumpalignment_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case OPT_HELP:
                std::cout << DUMPALIGNMENT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << DUMPALIGNMENT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
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

    if (die)
    {
        std::cout << "\n" << DUMPALIGNMENT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int dumpalignment_main(int argc, char** argv)
{
    parse_dumpalignment_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // Open readers
    FILE* read_fp = fopen(opt::reads_file.c_str(), "r");
    if(read_fp == NULL) {
        fprintf(stderr, "error: could not open %s for read\n", opt::reads_file.c_str());
        exit(EXIT_FAILURE);
    }

    gzFile gz_read_fp = gzdopen(fileno(read_fp), "r");
    if(gz_read_fp == NULL) {
        fprintf(stderr, "error: could not open %s using gzdopen\n", opt::reads_file.c_str());
        exit(EXIT_FAILURE);
    }

    // read input sequences, add to DB and convert to fasta
    int ret = 0;
    kseq_t* seq = kseq_init(gz_read_fp);
    while((ret = kseq_read(seq)) >= 0) {
        //fprintf(stdout, "#alignment for %s\n", seq->name.s);
        SquiggleRead sr(seq->name.s, read_db);

        size_t strand_idx = 0;
        std::vector<int> event_to_base_map(sr.events[strand_idx].size(), -1);
        for(size_t i = 0; i < sr.base_to_event_map.size(); ++i) {
            IndexPair ip = sr.base_to_event_map[i].indices[strand_idx];
            for(int j = ip.start; j != -1 && j <= ip.stop; ++j) {
                event_to_base_map[j] = i;
            }
        }

        int k = 6;
        int prev_idx = 0;
        for(size_t i = 0; i < sr.events[strand_idx].size(); ++i) {
            std::pair<size_t, size_t> sample_range = sr.get_event_sample_idx(strand_idx, i);
            int base_idx = prev_idx;
            std::string kmer = "NNNNNN";
            if(event_to_base_map[i] >= 0) {
                base_idx = event_to_base_map[i];
                kmer = sr.read_sequence.substr(base_idx, k);
            }

            fprintf(stdout, "%zu\t%d\t%d\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%s\n", i, base_idx, 0, sr.get_unscaled_level(i, strand_idx),
                                                                             sr.get_stdv(i, strand_idx), (float)sample_range.first,
                                                                             (float)(sample_range.second - sample_range.first), kmer.c_str());
        }
    }

    return EXIT_SUCCESS;
}
