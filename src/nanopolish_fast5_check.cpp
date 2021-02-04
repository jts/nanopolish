//---------------------------------------------------------
// Copyright 2020 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//
// nanopolish fast5-check - check fast5 files for common
// I/O problems
//
//---------------------------------------------------------
//

//
// Getopt
//
#define SUBPROGRAM "fast5-check"

#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>

#include <fast5.hpp>
#include "nanopolish_fast5_check.h"
#include "nanopolish_common.h"
#include "nanopolish_read_db.h"
#include "fs_support.hpp"
#include "nanopolish_fast5_io.h"

static const char *FAST5_CHECK_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2020 Ontario Institute for Cancer Research\n";

static const char *FAST5_CHECK_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] -r reads.fastq\n"
"Check whether the fast5 files are indexed correctly and readable by nanopolish\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display version\n"
"  -r, --reads                          file containing the basecalled reads\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose = 0;
    static std::string reads_file;
}

static const char* shortopts = "vr:";

enum {
    OPT_HELP = 1,
    OPT_VERSION,
    OPT_LOG_LEVEL,
};

static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { "version",                   no_argument,       NULL, OPT_VERSION },
    { "log-level",                 required_argument, NULL, OPT_LOG_LEVEL },
    { "verbose",                   no_argument,       NULL, 'v' },
    { "reads",                     required_argument, NULL, 'r' },
    { NULL, 0, NULL, 0 }
};

void parse_fast5_check_options(int argc, char** argv)
{
    bool die = false;
    std::vector< std::string> log_level;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case OPT_HELP:
                std::cout << FAST5_CHECK_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FAST5_CHECK_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                log_level.push_back(arg.str());
                break;
            case 'v': opt::verbose++; break;
            case 'r': arg >> opt::reads_file; break;
        }
    }

    if (argc - optind < 0) {
        std::cerr << SUBPROGRAM ": not enough arguments\n";
        die = true;
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << FAST5_CHECK_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

void check_read(fast5_file& f5_fh, const std::string& read_name)
{
    fast5_raw_scaling scaling = fast5_get_channel_params(f5_fh, read_name);
    if(scaling.digitisation != scaling.digitisation) {
        fprintf(stdout, "\t[read] ERROR: could not read scaling for %s\n", read_name.c_str());
    }
    raw_table rt = fast5_get_raw_samples(f5_fh, read_name, scaling);
    if(rt.n <= 0) {
        fprintf(stdout, "\t[read] ERROR: could not read raw samples for %s\n", read_name.c_str());
    } else {
        fprintf(stdout, "\t[read] OK: found %zu raw samples for %s\n", rt.n, read_name.c_str());
    }
    free(rt.raw);
    rt.raw = NULL;
}

int fast5_check_main(int argc, char** argv)
{
    parse_fast5_check_options(argc, argv);

    // Attempt to load the read_db
    ReadDB read_db;
    read_db.load(opt::reads_file);

    std::vector<std::string> fast5_files = read_db.get_unique_fast5s();
    fprintf(stdout, "The readdb file contains %zu fast5 files\n", fast5_files.size());

    for(size_t i = 0; i < fast5_files.size(); ++i) {

        fast5_file f5_fh = fast5_open(fast5_files[i]);
        if(fast5_is_open(f5_fh)) {
            fprintf(stdout, "[fast5] OK: opened %s\n", fast5_files[i].c_str());

            // check the individual reads in the file
            std::vector<std::string> reads = fast5_get_multi_read_groups(f5_fh);
            for(size_t j = 0; j < reads.size(); ++j) {
                std::string read_name = reads[j].substr(5);
                check_read(f5_fh, read_name);
            }
        } else {
            fprintf(stdout, "[fast5] ERROR: failed to open %s\n", fast5_files[i].c_str());
        }

        fast5_close(f5_fh);
    }
    return 0;
}
