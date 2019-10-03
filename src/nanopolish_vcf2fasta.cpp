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
"Usage: " PACKAGE_NAME " " SUBPROGRAM " -g draft.fa segment1.vcf segment2.vcf ...\n"
"Write a new genome sequence by introducing variants from the input files\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -g, --genome=FILE                    the input genome is in FILE\n"
"  -f, --fofn=FILE                      read the list of VCF files to use from FILE\n"
"      --skip-checks                    skip the sanity checks\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::vector<std::string> input_vcf_files;
    static std::string vcf_fofn;
    static std::string genome_file;
    static bool skip_checks = false;
}

static const char* shortopts = "g:f:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_SKIP_CHECKS };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "skip-checks",   no_argument,       NULL, OPT_SKIP_CHECKS },
    { "genome",        required_argument, NULL, 'g' },
    { "fofn",          required_argument, NULL, 'f' },
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
            case 'f': arg >> opt::vcf_fofn; break;
            case OPT_SKIP_CHECKS: opt::skip_checks = true; break;
            case OPT_HELP:
                std::cout << VCF2FASTA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << VCF2FASTA_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": -g/--genome file is required\n";
        die = true;
    }

    if (argc - optind < 1 && opt::vcf_fofn.empty()) {
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

    // add files from the fofn
    if(!opt::vcf_fofn.empty()) {
        std::ifstream infile(opt::vcf_fofn);
        std::string line;
        while(getline(infile, line)) {
            opt::input_vcf_files.push_back(line);
        }
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
                std::string window_key = "nanopolish_window=";
                size_t key_pos = line.find(window_key);
                if(key_pos != std::string::npos) {
                    window_str = line.substr(key_pos + window_key.size());
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

        if(!opt::skip_checks) {
            if(windows.empty()) {
                fprintf(stderr, "error: no polishing windows found for %s\n", contig.c_str());
                exit(EXIT_FAILURE);
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

            // check the first and last windows of the contig were polished
            if(windows[0].first != 0) {
                fprintf(stderr, "error: first %d bases are not covered by a polished window for contig %s.\n", windows[0].first, contig.c_str());
                window_check_ok = false;
            }

            int end_gap = contig_length - windows.back().second;
            if(end_gap > 500) {
                fprintf(stderr, "error: last %d bases are not covered by a polished window for contig %s.\n", end_gap, contig.c_str());
                window_check_ok = false;
            }
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

        auto& variants = variants_by_contig[contig];
        std::sort(variants.begin(), variants.end(), sortByPosition);

        // remove duplicate variants
        VariantKeyEqualityComp vkec;
        auto last = std::unique(variants.begin(), variants.end(), vkec);
        variants.erase(last, variants.end());

        assert(variants.size() < (1 << 30));
        uint32_t deleted_tag = 1 << 30;
        uint32_t variant_tag = 1 << 31;

        // make a vector holding either a literal character or an index to the variant that needs to be applied
        std::vector<uint32_t> consensus_record(length);
        for(size_t i = 0; i < length; ++i) {
            consensus_record[i] = seq[i];
        }

        size_t num_skipped = 0;
        size_t num_subs = 0;
        size_t num_insertions = 0;
        size_t num_deletions = 0;

        // update the consensus record according to the variants for this contig
        size_t applied_variants = 0;
        for(size_t variant_idx = 0; variant_idx < variants.size(); ++variant_idx) {
            const Variant& v = variants[variant_idx];

            // check if the variant record matches the reference sequence
            bool matches_ref = true;
            for(size_t i = 0; i < v.ref_seq.length(); ++i) {
                matches_ref = matches_ref && v.ref_seq[i] == consensus_record[v.ref_position + i];
            }

            if(!matches_ref) {
                num_skipped += 1;
                continue;
            }

            // mark the first base of the reference sequence as a variant and set the index
            consensus_record[v.ref_position] = variant_tag | variant_idx;

            // mark the subsequent bases of the reference as deleted
            for(size_t i = 1; i < v.ref_seq.length(); ++i) {
                consensus_record[v.ref_position + i] = deleted_tag;
            }

            num_subs += v.ref_seq.length() == v.alt_seq.length();
            num_insertions += v.ref_seq.length() < v.alt_seq.length();
            num_deletions += v.ref_seq.length() > v.alt_seq.length();
        }

        // write out the consensus record
        std::string out;
        out.reserve(length);
        for(size_t i = 0; i < length; ++i) {
            uint32_t r = consensus_record[i];
            if(r & variant_tag) {
                out.append(variants[r & ~variant_tag].alt_seq);
            } else if(r & ~deleted_tag) {
                out.append(1, r);
            } else {
                assert(r & deleted_tag);
            }
        }

        fprintf(stderr, "[vcf2fasta] rewrote contig %s with %zu subs, %zu ins, %zu dels (%zu skipped)\n", contig.c_str(), num_subs, num_insertions, num_deletions, num_skipped);
        fprintf(stdout, ">%s\n%s\n", contig.c_str(), out.c_str());

        free(seq);
        seq = NULL;
    }

    return 0;
}
