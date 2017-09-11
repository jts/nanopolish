//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//
// nanopolish index - build an index from FASTA/FASTQ file
// and the associated signal-level data
//
//---------------------------------------------------------
//

//
// Getopt
//
#define SUBPROGRAM "index"

#include <iostream>
#include <sstream>
#include <getopt.h>

#include "nanopolish_index.h"
#include "nanopolish_common.h"
#include "nanopolish_read_db.h"
#include "fs_support.hpp"
#include "logger.hpp"
#include "fast5.hpp"
#include "profiler.h"

static const char *INDEX_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *INDEX_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] -d nanopore_raw_file_directory reads.fastq\n"
"Build an index mapping from basecalled reads to the signals measured by the sequencer\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display version\n"
"  -v, --verbose                        display verbose output\n"
"  -d, --directory                      path to the directory containing the raw ONT signal files\n"
"  -f, --fast5-fofn                     file containing the paths to each fast5 file for the run\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose = 0;
    static std::string raw_file_directory;
    static std::string fast5_fofn;
    static std::string reads_file;
}
static std::ostream* os_p;

void index_file(ReadDB& read_db, const std::string& fn)
{
    PROFILE_FUNC("index_file")
    fast5::File* fp = new fast5::File(fn);
    if(fp->is_open()) {
        fast5::Raw_Samples_Params params = fp->get_raw_samples_params();
        std::string read_id = params.read_id;
        read_db.add_signal_path(read_id, fn);
    }
    delete fp;
} // process_file

void index_path(ReadDB& read_db, const std::string& path)
{
    fprintf(stderr, "Indexing %s\n", path.c_str());
    if (is_directory(path)) {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list) {
            if(fn == "." or fn == "..") {
                continue;
            }

            std::string full_fn = path + "/" + fn;
            if(is_directory(full_fn)) {
                // recurse
                index_path(read_db, full_fn);
                read_db.print_stats();
            } else if (full_fn.find(".fast5") != -1 && fast5::File::is_valid_file(full_fn)) {
                index_file(read_db, full_fn);
            }
        }
    }
} // process_path

void process_fast5_fofn(ReadDB& read_db, const std::string& fast5_fofn)
{
    //
    std::ifstream in_file(fast5_fofn.c_str());
    if(in_file.bad()) {
        fprintf(stderr, "error: could not file %s\n", fast5_fofn.c_str());
        exit(EXIT_FAILURE);
    }

    // read
    size_t count = 0;
    std::string filename;
    while(getline(in_file, filename)) {
        if(count++ % 1000 == 0) {
            fprintf(stderr, "Loaded %zu from fofn\n", count);
        }
        index_file(read_db, filename);
    }
}

static const char* shortopts = "vd:f:";

enum {
    OPT_HELP = 1,
    OPT_VERSION,
    OPT_LOG_LEVEL,
};

static const struct option longopts[] = {
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "log-level",          required_argument, NULL, OPT_LOG_LEVEL },
    { "verbose",            no_argument,       NULL, 'v' },
    { "directory",          required_argument, NULL, 'd' },
    { "fast5-fofn",         required_argument, NULL, 'f' },
    { NULL, 0, NULL, 0 }
};

void parse_index_options(int argc, char** argv)
{
    bool die = false;
    std::vector< std::string> log_level;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case OPT_HELP:
                std::cout << INDEX_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << INDEX_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                log_level.push_back(arg.str());
                break;
            case 'v': opt::verbose++; break;
            case 'd': arg >> opt::raw_file_directory; break;
            case 'f': arg >> opt::fast5_fofn; break;
        }
    }

    // set log levels
    auto default_level = (int)logger::warning + opt::verbose;
    logger::Logger::set_default_level(default_level);
    logger::Logger::set_levels_from_options(log_level, &std::clog);

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
        std::cout << "\n" << INDEX_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::reads_file = argv[optind++];
}

int index_main(int argc, char** argv)
{
    parse_index_options(argc, argv);
    
    // import read names, and possibly fast5 paths, from the fasta/fastq file
    ReadDB read_db;
    read_db.build(opt::reads_file);

    bool all_reads_have_paths = read_db.check_signal_paths();

    // if the input fastq did not contain a complete set of paths
    // use the fofn/directory provided to augment the index
    if(!all_reads_have_paths) {
        if(!opt::raw_file_directory.empty()) {
            index_path(read_db, opt::raw_file_directory);
        }

        if(!opt::fast5_fofn.empty()) {
            process_fast5_fofn(read_db, opt::fast5_fofn);
        }
    }

    read_db.print_stats();
    read_db.save();
    return 0;
}
