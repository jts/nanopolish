//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_consensus.cpp -- entry point to consensus functions
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
#include "htslib/htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_khmm_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "profiler.h"

//
// Getopt
//
#define SUBPROGRAM "eventalign"

static const char *EVENTALIGN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *EVENTALIGN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Compute a new consensus sequence for an assembly using a signal-level HMM\n"
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
    static int num_threads = 1;
}

static const char* shortopts = "r:b:g:t:w:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "reads",       required_argument, NULL, 'r' },
    { "bam",         required_argument, NULL, 'b' },
    { "genome",      required_argument, NULL, 'g' },
    { "window",      required_argument, NULL, 'w' },
    { "threads",     required_argument, NULL, 't' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Realign the read in event space
void realign_read(const Fast5Map& name_map, const bam1_t* record)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    std::string fast5_path = name_map.get_path(read_name);
    fprintf(stderr, "Realigning %s\n", read_name.c_str());

    // load read
    SquiggleRead sr(read_name, fast5_path);
}

void parse_eventalign_options(int argc, char** argv)
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
                std::cout << EVENTALIGN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EVENTALIGN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 0) {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } else if (argc - optind > 0) {
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
        std::cout << "\n" << EVENTALIGN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int eventalign_main(int argc, char** argv)
{
    parse_eventalign_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);

    /*
    // Parse the region string
    // Replace ":" and "-" with spaces to make it parseable with stringstream
    std::replace(opt::window.begin(), opt::window.end(), ':', ' ');
    std::replace(opt::window.begin(), opt::window.end(), '-', ' ');

    const int WINDOW_LENGTH = 10000;
    const int WINDOW_OVERLAP = 200;

    std::stringstream parser(opt::window);
    std::string contig;
    int start_window_id;
    int end_window_id;
    
    parser >> contig >> start_window_id >> end_window_id;
    */

    // Open the BAM and iterate over reads

    // load bam file
    htsFile* bam_fh = sam_open(opt::bam_file.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = opt::bam_file + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    //int contig_id = bam_name2id(hdr, contig_name.c_str());
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    // Initialize iteration
    bam1_t* record = bam_init1();
    
    int result;
    while((result = sam_read1(bam_fh, hdr, record)) >= 0) {

        realign_read(name_map, record);

        /*
        // parse alignments to reference
        std::vector<int> read_bases_for_anchors = 
            match_read_to_reference_anchors(record, start, end, stride);

        // Convert the read base positions into event indices for both strands
        HMMReadAnchorSet event_anchors;
        event_anchors.strand_anchors[T_IDX].resize(read_bases_for_anchors.size());
        event_anchors.strand_anchors[C_IDX].resize(read_bases_for_anchors.size());

        bool do_base_rc = bam_is_rev(record);
        bool template_rc = do_base_rc;
        bool complement_rc = !do_base_rc;

        for(size_t ai = 0; ai < read_bases_for_anchors.size(); ++ai) {

            int read_kidx = read_bases_for_anchors[ai];

            // read not aligned to this reference position
            if(read_kidx == -1) {
                continue;
            }

            if(do_base_rc)
                read_kidx = sr.flip_k_strand(read_kidx);

            int template_idx = sr.get_closest_event_to(read_kidx, T_IDX);
            int complement_idx = sr.get_closest_event_to(read_kidx, C_IDX);
            assert(template_idx != -1 && complement_idx != -1);

            event_anchors.strand_anchors[T_IDX][ai] = { template_idx, template_rc };
            event_anchors.strand_anchors[C_IDX][ai] = { complement_idx, complement_rc };
            
            // If this is not the last anchor, extract the sequence of the read
            // from this anchor to the next anchor as an alternative assembly
            if(ai < read_bases_for_anchors.size() - 1) {
                int start_kidx = read_bases_for_anchors[ai];
                int end_kidx = read_bases_for_anchors[ai + 1];
                int max_kidx = sr.read_sequence.size() - K;

                // flip
                if(do_base_rc) {
                    start_kidx = sr.flip_k_strand(start_kidx);
                    end_kidx = sr.flip_k_strand(end_kidx);
                    
                    // swap
                    int tmp = end_kidx;
                    end_kidx = start_kidx;
                    start_kidx = tmp;
                }

                // clamp values within range
                start_kidx = start_kidx >= 0 ? start_kidx : 0;
                end_kidx = end_kidx <= max_kidx ? end_kidx : max_kidx;
                
                std::string s = sr.read_sequence.substr(start_kidx, end_kidx - start_kidx + K);

                if(do_base_rc) {
                    s = reverse_complement(s);
                }

                if(ai >= read_substrings.size())
                    read_substrings.resize(ai + 1);

                read_substrings[ai].push_back(s);
            }
        }
        */
    }

    // cleanup
    //sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    bam_destroy1(record);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
}
