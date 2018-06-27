//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_vcf2fasta - write a new genome sequence
// by introducing variants from a set of vcf files
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
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
#include "nanopolish_common.h"
#include "nanopolish_variant.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_haplotype.h"

//
// Getopt
//
#define SUBPROGRAM "vcf2fasta"

static const char *VCF2FASTA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2018 Ontario Institute for Cancer Research\n";

static const char *VCF2FASTA_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] segment1.vcf segment2.vcf ...\n"
"Write a new genome sequence by introducing variants from the input files\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -g, --genome=FILE                    the input genome is in FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::vector<std::string> input_vcf_files;
    static std::string genome_file;
}

static const char* shortopts = "g:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { "genome",      required_argument, NULL, 'g' },
    { NULL, 0, NULL, 0 }
};

void parse_vcf2fasta_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'g': arg >> opt::genome_file; break;
            case OPT_HELP:
                std::cout << VCF2FASTA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << VCF2FASTA_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1) {
        std::cerr << SUBPROGRAM ": not enough arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << VCF2FASTA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    for(; optind < argc; ++optind) {
        opt::input_vcf_files.push_back(argv[optind]);
    }
}

int vcf2fasta_main(int argc, char** argv)
{
    parse_vcf2fasta_options(argc, argv);

    // Read genome file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // Read VCF files and gather variants for each contig and the polishing window coordinates
    std::map<std::string, std::vector<Variant>> variants_by_contig;
    std::map<std::string, std::vector<std::pair<int, int>>> windows_by_contig;

    for(const auto& filename : opt::input_vcf_files) {

        std::string window_str;
        std::vector<Variant> out;
        std::ifstream infile(filename);
        std::string line;
        while(getline(infile, line)) {

            // parse header
            if(line[0] == '#') {

                // check for window coordinates
                if(line.find("nanopolish_window") != std::string::npos) {
                    std::vector<std::string> fields = split(line, '=');
                    assert(fields.size() == 2);
                    window_str = fields[1];
                } 
            } else {
                Variant v(line);
                variants_by_contig[v.ref_name].push_back(v);
            }
        }

        if(window_str.empty()) {
            fprintf(stderr, "error: could not detect polishing window from input file %s\n", filename.c_str());
            exit(EXIT_FAILURE);
        }

        std::string window_contig;
        int window_start, window_end;
        parse_region_string(window_str, window_contig, window_start, window_end);
        windows_by_contig[window_contig].push_back(std::make_pair(window_start, window_end));
    }

    size_t n_contigs = faidx_nseq(fai);

    for(size_t contig_idx = 0; contig_idx < n_contigs; ++contig_idx) {
        std::string contig = faidx_iseq(fai, contig_idx);
        int contig_length = faidx_seq_len(fai, contig.c_str());

        // Confirm that all windows on this contig have been polished
        bool window_check_ok = true;
        auto& windows = windows_by_contig[contig];
        
        std::sort(windows.begin(), windows.end());
        if(windows[0].first != 0) {
            fprintf(stderr, "error: first %d bases are not covered by a polished window for contig %s.\n", windows[0].first, contig.c_str());
            window_check_ok = false;
        }

        for(size_t window_idx = 1; window_idx < windows.size(); ++window_idx) {
            int prev_start = windows[window_idx - 1].first;
            int prev_end = windows[window_idx - 1].second;
            int curr_start = windows[window_idx].first;
            int curr_end = windows[window_idx].second;
            if(curr_start > prev_end) {
                fprintf(stderr, "error: adjacent polishing windows do not overlap (%d-%d and %d-%d)\n", prev_start, prev_end, curr_start, curr_end);
                window_check_ok = false;
            }
        }

        int end_gap = contig_length - windows.back().second;
        if(end_gap > 500) {
            fprintf(stderr, "error: last %d bases are not covered by a polished window for contig %s.\n", end_gap, contig.c_str());
            window_check_ok = false;
        }

        if(!window_check_ok) {
            fprintf(stderr, "error: one or more polishing windows are missing. Please check that all nanopolish variants --consensus jobs ran to completion\n");
            exit(EXIT_FAILURE);
        }
        
        int length;
        char* seq = fai_fetch(fai, contig.c_str(), &length);
        if(length < 0) {
            fprintf(stderr, "error: could not fetch contig %s\n", contig.c_str());
            exit(EXIT_FAILURE);
        }
        std::string seq_str(seq);
        free(seq);
        seq = NULL;


        Haplotype derived_haplotype(contig, 0, seq_str);
        auto& variants = variants_by_contig[contig];
        std::sort(variants.begin(), variants.end(), sortByPosition);

        size_t applied_variants = 0;
        for(size_t variant_idx = 0; variant_idx < variants.size(); ++variant_idx) {
            applied_variants += derived_haplotype.apply_variant(variants[variant_idx]);
        }

        fprintf(stderr, "contig: %s %zu/%zu successfully added\n", contig.c_str(), applied_variants, variants.size());
        fprintf(stdout, ">%s\n%s\n", contig.c_str(), derived_haplotype.get_sequence().c_str());
    }

    return 0;
}
