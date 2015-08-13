//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_call_variants -- find variants wrt a reference
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
#include <queue>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_khmm_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "profiler.h"
#include "progress.h"
#include "stdaln.h"

// Macros
#define max3(x,y,z) std::max(std::max(x,y), z)

// Flags to turn on/off debugging information

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1
//#define DEBUG_TRANSITION 1
//#define DEBUG_PATH_SELECTION 1
//#define DEBUG_SINGLE_SEGMENT 1
//#define DEBUG_SHOW_TOP_TWO 1
//#define DEBUG_SEGMENT_ID 193
//#define DEBUG_BENCHMARK 1

//
// Getopt
//
#define SUBPROGRAM "consensus"

static const char *CONSENSUS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CONSENSUS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Compute a new consensus sequence for an assembly using a signal-level HMM\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"      --snps                           only call SNPs\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -e, --event-bam=FILE                 the events aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -o, --outfile=FILE                   write result to FILE [default: stdout]\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string event_bam_file;
    static std::string genome_file;
    static std::string output_file;
    static std::string window;
    static int snps_only = 0;
    static int show_progress = 0;
    static int num_threads = 1;
}

static const char* shortopts = "r:b:g:t:w:o:e:m:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VCF, OPT_PROGRESS, OPT_SNPS_ONLY };

static const struct option longopts[] = {
    { "verbose",           no_argument,       NULL, 'v' },
    { "reads",             required_argument, NULL, 'r' },
    { "bam",               required_argument, NULL, 'b' },
    { "event-bam",         required_argument, NULL, 'e' },
    { "genome",            required_argument, NULL, 'g' },
    { "window",            required_argument, NULL, 'w' },
    { "outfile",           required_argument, NULL, 'o' },
    { "threads",           required_argument, NULL, 't' },
    { "min-read-evidence", required_argument, NULL, 'm' },
    { "snps",              no_argument,       NULL, OPT_SNPS_ONLY },
    { "progress",          no_argument,       NULL, OPT_PROGRESS },
    { "help",              no_argument,       NULL, OPT_HELP },
    { "version",           no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int get_contig_length(const std::string& contig)
{
    faidx_t *fai = fai_load(opt::genome_file.c_str());
    int len = faidx_seq_len(fai, contig.c_str());
    fai_destroy(fai);
    return len;
}

double parallel_score(const std::string& sequence, const std::vector<HMMInputData>& event_sequences)
{
    double score = 0.0f;

    #pragma omp parallel for
    for(size_t i = 0; i < event_sequences.size(); ++i) {
        double s = profile_hmm_score(sequence, event_sequences[i]);

        #pragma omp critical
        score += s;
    }
    return score;
}

Haplotype call_variants_for_region(const std::string& contig, int region_start, int region_end)
{
    const int BUFFER = 20;

    if(region_start < BUFFER)
        region_start = BUFFER;

    // load the region, accounting for the buffering
    AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file);
    alignments.load_region(contig, region_start - BUFFER, region_end + BUFFER);
    Haplotype derived_haplotype(contig,
                                alignments.get_region_start(),
                                alignments.get_reference());

    // Step 1. Discover putative variants across the whole region
    std::vector<Variant> candidate_variants = alignments.get_variants_in_region(contig, region_start, region_end, 0.2, 20);

    // Step 2. Add variants to the haplotypes
    size_t calling_span = 15;
    size_t curr_variant_idx = 0;
    while(curr_variant_idx < candidate_variants.size()) {
 
        // Group the variants that are within calling_span bases of each other       
        size_t end_variant_idx = curr_variant_idx + 1;
        while(end_variant_idx < candidate_variants.size()) {
            int distance = candidate_variants[end_variant_idx].ref_position - 
                           candidate_variants[end_variant_idx - 1].ref_position;
            if(distance > calling_span)
                break;
            end_variant_idx++;
        }
    
        size_t num_variants = end_variant_idx - curr_variant_idx;
        int calling_start = candidate_variants[curr_variant_idx].ref_position - calling_span;
        int calling_end = candidate_variants[end_variant_idx - 1].ref_position +
                          candidate_variants[end_variant_idx - 1].ref_seq.length() +
                          calling_span;
        int calling_size = calling_end - calling_start;

        if(opt::verbose > 2) {
            fprintf(stderr, "%zu variants in span [%d %d]\n", num_variants, calling_start, calling_end);
        }
        
        // Only try to call variants if there is a reasonable amount and the window is not too large
        if(num_variants <= 10 && calling_size <= 100) {

            // Subset the haplotype to the region we are calling
            Haplotype calling_haplotype = 
                derived_haplotype.substr_by_reference(calling_start, calling_end);
            
            // Get the events for the calling region
            std::vector<HMMInputData> event_sequences = 
                alignments.get_event_subsequences(contig, calling_start, calling_end);
            
            // Subset the variants
            std::vector<Variant> calling_variants(candidate_variants.begin() + curr_variant_idx, 
                                                  candidate_variants.begin() + end_variant_idx);
            
            // Select the best set of variants
            std::vector<Variant> selected_variants = 
                select_variant_set(calling_variants, calling_haplotype, event_sequences);

            // Apply them to the final haplotype
            for(size_t vi = 0; vi < selected_variants.size(); vi++) {
                derived_haplotype.apply_variant(selected_variants[vi]);

                if(opt::verbose > 1) {
                    selected_variants[vi].write_vcf(stderr);
                }
            }
        }

        // advance to start of next region
        curr_variant_idx = end_variant_idx;
    }

    return derived_haplotype;
}

void parse_call_variants_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'e': arg >> opt::event_bam_file; break;
            case 'w': arg >> opt::window; break;
            case 'o': arg >> opt::output_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_SNPS_ONLY: opt::snps_only = 1; break;
            case OPT_PROGRESS: opt::show_progress = 1; break;
            case OPT_HELP:
                std::cout << CONSENSUS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONSENSUS_VERSION_MESSAGE;
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
        std::cout << "\n" << CONSENSUS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int call_variants_main(int argc, char** argv)
{
    parse_call_variants_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    // Parse the window string
    // Replace ":" and "-" with spaces to make it parseable with stringstream
    std::replace(opt::window.begin(), opt::window.end(), ':', ' ');
    std::replace(opt::window.begin(), opt::window.end(), '-', ' ');

    std::stringstream parser(opt::window);
    std::string contig;
    int start_base;
    int end_base;
    
    parser >> contig >> start_base >> end_base;
    end_base = std::min(end_base, get_contig_length(contig) - 1);

    FILE* out_fp;
    if(!opt::output_file.empty()) {
        out_fp = fopen(opt::output_file.c_str(), "w");
    } else {
        out_fp = stdout;
    }

    Variant::write_vcf_header(out_fp);

    fprintf(stderr, "TODO: train model\n");
    fprintf(stderr, "TODO: filter data\n");

    Haplotype haplotype = call_variants_for_region(contig, start_base, end_base);

    std::vector<Variant> variants = haplotype.get_variants();
    for(size_t vi = 0; vi < variants.size(); vi++) {
        variants[vi].write_vcf(out_fp);
    }

    if(out_fp != stdout) {
        fclose(out_fp);
    }
}
