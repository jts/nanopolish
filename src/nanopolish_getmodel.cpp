//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_getmodel.h - write the pore model for a read
// to stdout
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
#include <fast5.hpp>
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "profiler.h"

//
// Getopt
//
#define SUBPROGRAM "getmodel"

static const char *GETMODEL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *GETMODEL_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] read.fast5\n"
"Write the pore models for the given read to stdout\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string input_file;
}

static const char* shortopts = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_getmodel_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GETMODEL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GETMODEL_VERSION_MESSAGE;
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
        std::cout << "\n" << GETMODEL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::input_file = argv[optind++];
}

int getmodel_main(int argc, char** argv)
{
    parse_getmodel_options(argc, argv);

    fast5::File f(opt::input_file);

    printf("strand\tkmer\tmodel_mean\tmodel_stdv\n");

    for(size_t si = 0; si < 2; ++si) {
        auto g_l = f.get_basecall_strand_group_list(si);
        for (const auto& g : g_l) {
            if (not f.have_basecall_model(si, g)) continue;
            PoreModel pm(&f, si, g);
            assert(pm.get_num_states() == gDNAAlphabet.get_num_strings(pm.k));
            std::string kmer(pm.k, 'A');
            for(size_t ki = 0; ki < pm.get_num_states(); ++ki) {
                PoreModelStateParams params = pm.get_parameters(ki);
                printf("%c\t%s\t%.2lf\t%.2lf\n",
                       (si == 0? 't' : 'c'),
                       kmer.c_str(),
                       params.level_mean,
                       params.level_stdv);
                gDNAAlphabet.lexicographic_next(kmer); // advance kmer
            }
            break;
        }
    }
    return 0;
}
