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
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <slow5/slow5.h>
#include "nanopolish_index.h"
#include "nanopolish_common.h"
#include "nanopolish_read_db.h"
#include "fs_support.hpp"
#include "logger.hpp"
#include "profiler.h"
#include "nanopolish_fast5_io.h"

static const char *INDEX_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *INDEX_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] -d nanopore_raw_file_directory reads.fastq\n"
"       " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --slow5 reads.blow5 reads.fastq\n"
"Build an index mapping from basecalled reads to the signals measured by the sequencer\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display version\n"
"  -v, --verbose                        display verbose output\n"
"      --slow5 FILE                     s/blow5 file (-d, -s, -f options below are only for fast5)\n"
"  -d, --directory                      path to the directory containing the raw ONT signal files. This option can be given multiple times.\n"
"  -s, --sequencing-summary             the sequencing summary file from albacore, providing this option will make indexing much faster\n"
"  -f, --summary-fofn                   file containing the paths to the sequencing summary files (one per line)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose = 0;
    static std::vector<std::string> raw_file_directories;
    static std::string reads_file;
    static std::vector<std::string> sequencing_summary_files;
    static std::string sequencing_summary_fofn;
    static std::string slow5_file;
}
static std::ostream* os_p;

//
void index_file_from_map(ReadDB& read_db, const std::string& fn, const std::multimap<std::string, std::string>& fast5_to_read_name_map)
{
    PROFILE_FUNC("index_file_from_map")

    // Check if this fast5 file is in the map
    size_t last_dir_pos = fn.find_last_of("/");
    std::string fast5_basename = last_dir_pos == std::string::npos ? fn : fn.substr(last_dir_pos + 1);

    auto range = fast5_to_read_name_map.equal_range(fast5_basename);
    for(auto iter = range.first; iter != range.second; ++iter) {
        if(read_db.has_read(iter->second)) {
            read_db.add_signal_path(iter->second.c_str(), fn);
        }
    }
}

//
void index_file_from_fast5(ReadDB& read_db, const std::string& fn)
{
    PROFILE_FUNC("index_file_from_fast5")

    fast5_file f5_file = fast5_open(fn);
    if(!fast5_is_open(f5_file)) {
        fprintf(stderr, "could not open fast5 file: %s\n", fn.c_str());
    }

    if(f5_file.is_multi_fast5) {
        std::vector<std::string> read_groups = fast5_get_multi_read_groups(f5_file);
        std::string prefix = "read_";
        for(size_t group_idx = 0; group_idx < read_groups.size(); ++group_idx) {
            std::string group_name = read_groups[group_idx];
            if(group_name.find(prefix) == 0) {
                std::string read_id = group_name.substr(prefix.size());
                read_db.add_signal_path(read_id, fn);
            }
        }
    } else {
        std::string read_id = fast5_get_read_id_single_fast5(f5_file);
        if(read_id != "") {
            read_db.add_signal_path(read_id, fn);
        }
    }
    fast5_close(f5_file);
}

//
void index_path(ReadDB& read_db, const std::string& path, const std::multimap<std::string, std::string>& fast5_to_read_name_map)
{
    fprintf(stderr, "[readdb] indexing %s\n", path.c_str());
    if (is_directory(path)) {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list) {
            if(fn == "." or fn == "..") {
                continue;
            }

            std::string full_fn = path + "/" + fn;
            bool is_fast5 = full_fn.find(".fast5") != -1;
            bool in_map = fast5_to_read_name_map.find(fn) != fast5_to_read_name_map.end();

            // JTS 04/19: is_directory is painfully slow so we first check if the file is in the name map
            // if it is, it is definitely not a directory so we can skip the system call
            if(!in_map && is_directory(full_fn)) {
                // recurse
                index_path(read_db, full_fn, fast5_to_read_name_map);
            } else if (is_fast5) {
                if(in_map) {
                    index_file_from_map(read_db, full_fn, fast5_to_read_name_map);
                } else {
                    index_file_from_fast5(read_db, full_fn);
                }
            }
        }
    }
}

// read sequencing summary files from the fofn and add them to the list
void process_summary_fofn()
{
    if(opt::sequencing_summary_fofn.empty()) {
        return;
    }

    // open
    std::ifstream in_file(opt::sequencing_summary_fofn.c_str());
    if(in_file.bad()) {
        fprintf(stderr, "error: could not file %s\n", opt::sequencing_summary_fofn.c_str());
        exit(EXIT_FAILURE);
    }

    // read
    std::string filename;
    while(getline(in_file, filename)) {
        opt::sequencing_summary_files.push_back(filename);
    }
}

//
void exit_bad_header(const std::string& str, const std::string& filename)
{
    fprintf(stderr, "Could not find %s column in the header of %s\n", str.c_str(), filename.c_str());
    exit(EXIT_FAILURE);
}

//
void parse_sequencing_summary(const std::string& filename, std::multimap<std::string, std::string>& out_map)
{
    // open
    std::ifstream in_file(filename.c_str());
    if(!in_file.good()) {
        fprintf(stderr, "error: could not read file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }

    // read header to get the column index of the read and file name
    std::string header;
    getline(in_file, header);
    std::vector<std::string> fields = split(header, '\t');

    const std::string READ_NAME_STR = "read_id";
    size_t filename_idx = -1;
    size_t read_name_idx = -1;

    for(size_t i = 0; i < fields.size(); ++i) {
        if(fields[i] == READ_NAME_STR) {
            read_name_idx = i;
        }

        // 19/11/05: support live basecalling summary files
        if(fields[i] == "filename" || fields[i] == "filename_fast5") {
            filename_idx = i;
        }
    }

    if(filename_idx == -1) {
        exit_bad_header("fast5 filename", filename);
    }

    if(read_name_idx == -1) {
        exit_bad_header(READ_NAME_STR, filename);
    }

    // read records and add to map
    std::string line;
    while(getline(in_file, line)) {
        fields = split(line, '\t');
        std::string fast5_filename = fields[filename_idx];
        std::string read_name = fields[read_name_idx];
        out_map.insert(std::make_pair(fast5_filename, read_name));
    }
}

static const char* shortopts = "vd:f:s:";

enum {
    OPT_HELP = 1,
    OPT_VERSION,
    OPT_LOG_LEVEL,
    OPT_SLOW5_FILE,
};

static const struct option longopts[] = {
    { "help",                      no_argument,       NULL, OPT_HELP },
    { "version",                   no_argument,       NULL, OPT_VERSION },
    { "log-level",                 required_argument, NULL, OPT_LOG_LEVEL },
    { "slow5",                     required_argument, NULL, OPT_SLOW5_FILE},
    { "verbose",                   no_argument,       NULL, 'v' },
    { "directory",                 required_argument, NULL, 'd' },
    { "sequencing-summary-file",   required_argument, NULL, 's' },
    { "summary-fofn",              required_argument, NULL, 'f' },
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
            case  OPT_SLOW5_FILE : arg >> opt::slow5_file; break;
            case 'v': opt::verbose++; break;
            case 's': opt::sequencing_summary_files.push_back(arg.str()); break;
            case 'd': opt::raw_file_directories.push_back(arg.str()); break;
            case 'f': arg >> opt::sequencing_summary_fofn; break;
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

    if(!opt::slow5_file.empty()){
        if (!opt::sequencing_summary_fofn.empty()) {
            std::cerr << "Option --summary-fofn will be ignored in slow5 mode\n";
        }
        if (!opt::sequencing_summary_files.empty()) {
            std::cerr << "Option --sequencing-summary-file will be ignored in slow5 mode\n";
        }
        if (!opt::raw_file_directories.empty()) {
            std::cerr << "Option -d will be ignored in slow5 mode\n";
        }
    }

    opt::reads_file = argv[optind++];
}

void clean_fast5_map(std::multimap<std::string, std::string>& mmap)
{
    std::map<std::string, int> fast5_read_count;
    for(const auto& iter : mmap) {
        fast5_read_count[iter.first]++;
    }

    // calculate the mode of the read counts per fast5
    std::map<size_t, size_t> read_count_distribution;
    for(const auto& iter : fast5_read_count) {
        read_count_distribution[iter.second]++;
    }

    size_t mode = 0;
    size_t mode_count = 0;
    for(const auto& iter : read_count_distribution) {
        if(iter.second > mode_count) {
            mode = iter.first;
            mode_count = iter.second;
        }
    }

    // this is the default number of reads per fast5 and we expect the summary file to support that
    // this code is to detect whether the user selected a different number of reads (ARTIC recommends 1000)
    // so the summary file still works
    int EXPECTED_ENTRIES_PER_FAST5 = 4000;
    if(mode_count > 20) {
        EXPECTED_ENTRIES_PER_FAST5 = mode;
    }

    std::vector<std::string> invalid_entries;
    for(const auto& iter : fast5_read_count) {
        if(iter.second > EXPECTED_ENTRIES_PER_FAST5) {
            //fprintf(stderr, "warning: %s has %d entries in the summary and will be indexed the slow way\n", iter.first.c_str(), iter.second);
            invalid_entries.push_back(iter.first);
        }
    }

    if(invalid_entries.size() > 0) {
        fprintf(stderr, "warning: detected invalid summary file entries for %zu of %zu fast5 files\n", invalid_entries.size(), fast5_read_count.size());
        fprintf(stderr, "These files will be indexed without using the summary file, which is slow.\n");
    }

    for(const auto& fast5_name : invalid_entries) {
        mmap.erase(fast5_name);
    }
}

int index_main(int argc, char** argv)
{
    parse_index_options(argc, argv);

    if(opt::slow5_file.empty()) {
        // Read a map from fast5 files to read name from the sequencing summary files (if any)
        process_summary_fofn();
        std::multimap<std::string, std::string> fast5_to_read_name_map;
        for(const auto& ss_filename : opt::sequencing_summary_files) {
            if(opt::verbose > 2) {
                fprintf(stderr, "summary: %s\n", ss_filename.c_str());
            }
            parse_sequencing_summary(ss_filename, fast5_to_read_name_map);
        }

        // Detect non-unique fast5 file names in the summary file
        // This occurs when a file with the same name (abc_0.fast5) appears in both fast5_pass
        // and fast5_fail. This will be fixed by ONT but in the meantime we check for
        // fast5 files that have a non-standard number of reads (>4000) and remove them
        // from the map. These fast5s will be indexed the slow way.
        clean_fast5_map(fast5_to_read_name_map);

        // import read names, and possibly fast5 paths, from the fasta/fastq file
        ReadDB read_db;
        read_db.build(opt::reads_file);
        bool all_reads_have_paths = read_db.check_signal_paths();

        // if the input fastq did not contain a complete set of paths
        // use the fofn/directory provided to augment the index
        if(!all_reads_have_paths) {

            for(const auto& dir_name : opt::raw_file_directories) {
                index_path(read_db, dir_name, fast5_to_read_name_map);
            }
        }

        size_t num_with_path = read_db.get_num_reads_with_path();
        size_t num_reads = read_db.get_num_reads();
        if(num_with_path == 0) {
            fprintf(stderr,"No fast5 files found\n");
            exit(EXIT_FAILURE);
        } else {
            read_db.print_stats();
            read_db.save();
            if (num_with_path < num_reads * 0.99 ) {
                fprintf(stderr,"fast5 files could not be located for %ld reads\n",num_reads-num_with_path);
            }
        }
    }
    else{
        slow5_file_t *sp = slow5_open(opt::slow5_file.c_str(), "r");
        if(sp == NULL){
            fprintf(stderr,"Error in opening slowfile %s\n", opt::slow5_file.c_str());
            exit(EXIT_FAILURE);
        }
        int ret = 0;
        ret = slow5_idx_create(sp);
        if(ret < 0) {
            fprintf(stderr, "Error in creating index for slow5 file %s\n", opt::slow5_file.c_str());
            exit(EXIT_FAILURE);
        }
        slow5_close(sp);

        ReadDB read_db;
        read_db.set_slow5_mode(true);
        read_db.build(opt::reads_file);
        read_db.add_signal_path("*", opt::slow5_file.c_str());
        read_db.save();

    }
    return 0;
}
