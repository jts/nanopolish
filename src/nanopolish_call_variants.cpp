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
#include <fstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_pore_model_set.h"
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
#define SUBPROGRAM "variants"

static const char *CONSENSUS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CONSENSUS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Find SNPs using a signal-level HMM\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"      --snps                           only call SNPs\n"
"  -w, --window=STR                     find variants in window STR (format: ctg:start-end)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the reference genome are in bam FILE\n"
"  -e, --event-bam=FILE                 the events aligned to the reference genome are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -o, --outfile=FILE                   write result to FILE [default: stdout]\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -m, --min-candidate-frequency=F      alternative bases in F proporation of aligned reads are candidate variants (default 0.2)\n"
"  -c, --candidates=VCF                 read variant candidates from VCF, rather than discovering them from aligned reads\n"
"      --calculate-all-support          when making a call, also calculate the support of the 3 other possible bases\n"
"      --models-fofn=FILE               read alternative k-mer models from FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string event_bam_file;
    static std::string genome_file;
    static std::string output_file;
    static std::string candidates_file;
    static std::string models_fofn;
    static std::string window;
    static std::string alternative_model_type = DEFAULT_MODEL_TYPE;
    static double min_candidate_frequency = 0.2f;
    static int calculate_all_support = false;
    static int snps_only = 0;
    static int show_progress = 0;
    static int num_threads = 1;
}

static const char* shortopts = "r:b:g:t:w:o:e:m:c:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VCF, OPT_PROGRESS, OPT_SNPS_ONLY, OPT_CALC_ALL_SUPPORT, OPT_MODELS_FOFN };

static const struct option longopts[] = {
    { "verbose",                 no_argument,       NULL, 'v' },
    { "reads",                   required_argument, NULL, 'r' },
    { "bam",                     required_argument, NULL, 'b' },
    { "event-bam",               required_argument, NULL, 'e' },
    { "genome",                  required_argument, NULL, 'g' },
    { "window",                  required_argument, NULL, 'w' },
    { "outfile",                 required_argument, NULL, 'o' },
    { "threads",                 required_argument, NULL, 't' },
    { "min-candidate-frequency", required_argument, NULL, 'm' },
    { "candidates",              required_argument, NULL, 'c' },
    { "models-fofn",             required_argument, NULL, OPT_MODELS_FOFN },
    { "calculate-all-support",   no_argument,       NULL, OPT_CALC_ALL_SUPPORT },
    { "snps",                    no_argument,       NULL, OPT_SNPS_ONLY },
    { "progress",                no_argument,       NULL, OPT_PROGRESS },
    { "help",                    no_argument,       NULL, OPT_HELP },
    { "version",                 no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int get_contig_length(const std::string& contig)
{
    faidx_t *fai = fai_load(opt::genome_file.c_str());
    int len = faidx_seq_len(fai, contig.c_str());
    fai_destroy(fai);
    return len;
}

void annotate_with_all_support(std::vector<Variant>& variants,
                               Haplotype base_haplotype,
                               const std::vector<HMMInputData>& input,
                               const uint32_t alignment_flags)

{
    for(size_t vi = 0; vi < variants.size(); vi++) {
        
        // Generate a haplotype containing every variant in the set except for vi
        Haplotype test_haplotype = base_haplotype;
        for(size_t vj = 0; vj < variants.size(); vj++) {

            // do not apply the variant we are testing
            if(vj == vi) {
                continue;
            }
            test_haplotype.apply_variant(variants[vj]);
        }
        
        // Make a vector of four haplotypes, one per base
        std::vector<Haplotype> curr_haplotypes;
        Variant tmp_variant = variants[vi];
        for(size_t bi = 0; bi < 4; ++bi) {
            tmp_variant.alt_seq = "ACGT"[bi];
            Haplotype tmp = test_haplotype;
            tmp.apply_variant(tmp_variant);
            curr_haplotypes.push_back(tmp);
        }

        // Test all reads against the 4 haplotypes
        std::vector<int> support_count(4, 0);

        for(size_t input_idx = 0; input_idx < input.size(); ++input_idx) {
            double best_score = -INFINITY;
            size_t best_hap_idx = 0;

            // calculate which haplotype this read supports best
            for(size_t hap_idx = 0; hap_idx < curr_haplotypes.size(); ++hap_idx) {
                double score = profile_hmm_score(curr_haplotypes[hap_idx].get_sequence(), input[input_idx], alignment_flags);
                if(score > best_score) {
                    best_score = score;
                    best_hap_idx = hap_idx;
                }
            }
            support_count[best_hap_idx] += 1;
        }

        std::stringstream ss;
        for(size_t bi = 0; bi < 4; ++bi) {
            ss << support_count[bi] / (double)input.size() << (bi != 3 ? "," : "");
        }

        variants[vi].add_info("AllSupportFractions", ss.str());
    }
}

std::vector<Variant> get_variants_from_vcf(const std::string& filename, 
                                           const std::string& contig,
                                           int region_start,
                                           int region_end)
{
    std::vector<Variant> out;
    std::ifstream infile(filename);
    std::string line;
    while(getline(infile, line)) {
        
        // skip headers
        if(line[0] == '#') {
            continue;
        }
        
        // parse variant
        Variant v(line);

        if(v.ref_name == contig &&
           (int)v.ref_position >= region_start &&
           (int)v.ref_position <= region_end) 
        {
            out.push_back(v);
        }
    }
    return out;
}

Haplotype call_variants_for_region(const std::string& contig, int region_start, int region_end)
{
    const int BUFFER = 20;
    uint32_t alignment_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;
    if(region_start < BUFFER)
        region_start = BUFFER;

    // load the region, accounting for the buffering
    AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file);
    
    if(!opt::alternative_model_type.empty()) {
        alignments.set_alternative_model_type(opt::alternative_model_type);
    }

    alignments.load_region(contig, region_start - BUFFER, region_end + BUFFER);
    Haplotype derived_haplotype(contig,
                                alignments.get_region_start(),
                                alignments.get_reference());

    // Step 1. Discover putative variants across the whole region
    std::vector<Variant> candidate_variants;
    if(opt::candidates_file.empty()) {
        candidate_variants = alignments.get_variants_in_region(contig, region_start, region_end, opt::min_candidate_frequency, 20);
    } else {
        candidate_variants = get_variants_from_vcf(opt::candidates_file, contig, region_start, region_end);
    }

    // Step 2. Add variants to the haplotypes
    int calling_span = 10;
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
        if(num_variants <= 15 && calling_size <= 100) {

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
                select_variant_set(calling_variants, calling_haplotype, event_sequences, alignment_flags);

            // optionally annotate each variant with fraction of reads supporting A,C,G,T at this position
            if(opt::calculate_all_support) {
                annotate_with_all_support(selected_variants, calling_haplotype, event_sequences, alignment_flags);
            }

            // Apply them to the final haplotype
            for(size_t vi = 0; vi < selected_variants.size(); vi++) {

                derived_haplotype.apply_variant(selected_variants[vi]);

                if(opt::verbose > 1) {
                    selected_variants[vi].write_vcf(stderr);
                }
            }
        } else {
            fprintf(stderr, "Warning: %zu variants in span, region not called [%d %d]\n", num_variants, calling_start, calling_end);
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
            case 'm': arg >> opt::min_candidate_frequency; break;
            case 'c': arg >> opt::candidates_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_MODELS_FOFN: arg >> opt::models_fofn; break;
            case OPT_CALC_ALL_SUPPORT: opt::calculate_all_support = 1; break;
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

    if(opt::models_fofn.empty()) {
        std::cerr << SUBPROGRAM ": a --models file must be provided\n";
        die = true;
    } else {
        // initialize the model set from the fofn
        PoreModelSet::initialize(opt::models_fofn);
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
    return 0;
}
