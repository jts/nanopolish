//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_polya_estimator.cpp -- estimate the length
// of poly-A tails for each read
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
#include "htslib/faidx.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_read_db.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_polya_estimator.h"
#include "nanopolish_raw_loader.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

using namespace std::placeholders;

//
// Getopt
//
#define SUBPROGRAM "polya"

static const char *POLYA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *POLYA_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Estimate the length of the poly-A tail on direct RNA reads\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:w:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "reads",            required_argument, NULL, 'r' },
    { "bam",              required_argument, NULL, 'b' },
    { "genome",           required_argument, NULL, 'g' },
    { "window",           required_argument, NULL, 'w' },
    { "threads",          required_argument, NULL, 't' },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_polya_options(int argc, char** argv)
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
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << POLYA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << POLYA_VERSION_MESSAGE;
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

    if (die)
    {
        std::cout << "\n" << POLYA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

//
void estimate_polya_for_single_read(const ReadDB& read_db,
                                    const faidx_t* fai,
                                    FILE* out_fp,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);


    // Get the length of the suffix of the read that has been softclipped
    size_t n_cigar = record->core.n_cigar;
    uint32_t suffix_cigar = bam_get_cigar(record)[n_cigar - 1];
    uint32_t suffix_clip = bam_cigar_oplen(suffix_cigar);
    fprintf(out_fp, "read: %s suffix clip: %zu\n", read_name.c_str(), suffix_clip);

    // load read
    SquiggleRead sr(read_name, read_db);

    // Partition read into an adapter (clipped suffix) and valid transcript part (everything not clipped)
    const std::string& read_sequence = sr.read_sequence;
    std::string transcript = read_sequence.substr(0, read_sequence.length() - suffix_clip);
    std::string suffix = read_sequence.substr(read_sequence.length() - suffix_clip);
    std::string polya_placeholder = "AAAAAAAAAA";
    std::string sequenced_transcript = transcript + polya_placeholder + suffix;

    // Align events to polya + suffix
    std::vector<AlignedPair> alignment = adaptive_banded_simple_event_align(sr, sequenced_transcript);
    for(size_t i = 0; i < alignment.size(); ++i) {
        fprintf(out_fp, "%zu %zu %zu %s\n", i, alignment[i].ref_pos, alignment[i].read_pos, sequenced_transcript.substr(alignment[i].ref_pos, 6).c_str());
    }
}

//
int polya_main(int argc, char** argv)
{
    parse_polya_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(estimate_polya_for_single_read, std::ref(read_db), std::ref(fai), stdout, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    return EXIT_SUCCESS;
}
