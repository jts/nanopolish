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
"      --consensus                      run in consensus calling mode\n"
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
    static std::string alternative_model_type = "reftrained";
    static double min_candidate_frequency = 0.2f;
    static int calculate_all_support = false;
    static int snps_only = 0;
    static int show_progress = 0;
    static int num_threads = 1;
    static int calibrate = 0;
    static int consensus_mode = 0;
}

static const char* shortopts = "r:b:g:t:w:o:e:m:c:v";

enum { OPT_HELP = 1,
       OPT_VERSION,
       OPT_VCF,
       OPT_PROGRESS,
       OPT_SNPS_ONLY,
       OPT_CALC_ALL_SUPPORT,
       OPT_CONSENSUS,
       OPT_MODELS_FOFN };

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
    { "consensus",               no_argument,       NULL, OPT_CONSENSUS },
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

std::vector<Variant> generate_candidate_single_base_edits(const AlignmentDB& alignments,
                                                          int region_start,
                                                          int region_end,
                                                          int calling_span,
                                                          uint32_t alignment_flags)
{
    std::vector<Variant> out_variants;

    std::string contig = alignments.get_region_contig();

    // Add all positively-scoring single-base changes into the candidate set
    for(size_t i = region_start; i < region_end; ++i) {

        std::vector<Variant> sb_variants;

        for(size_t j = 0; j < 4; ++j) {
            // Substitutions
            Variant v;
            v.ref_name = contig;
            v.ref_position = i;
            v.ref_seq = alignments.get_reference_substring(contig, i, i);
            v.alt_seq = "ACGT"[j];

            if(v.ref_seq != v.alt_seq) {
                sb_variants.push_back(v);
            }

            // Insertions
            v.alt_seq = v.ref_seq + "ACGT"[j];
            // ignore insertions of the type "A" -> "AA" as these are redundant
            if(v.alt_seq[1] != v.ref_seq[0]) {
                sb_variants.push_back(v);
            }
        }

        // deletion
        Variant del;
        del.ref_name = contig;
        del.ref_position = i - 1;
        del.ref_seq = alignments.get_reference_substring(contig, i - 1, i);
        del.alt_seq = del.ref_seq[0];

        // ignore deletions of the type "AA" -> "A" as these are redundant
        if(del.alt_seq[0] != del.ref_seq[1]) {
            sb_variants.push_back(del);
        }

        // Screen for positively scoring variants
        int calling_start = sb_variants[0].ref_position - calling_span;
        int calling_end = sb_variants[0].ref_position + calling_span;

        Haplotype test_haplotype(contig,
                                 calling_start,
                                 alignments.get_reference_substring(contig, calling_start, calling_end));

        std::vector<HMMInputData> event_sequences =
            alignments.get_event_subsequences(contig, calling_start, calling_end);

        std::vector<Variant> passed_variants = select_positive_scoring_variants(sb_variants, test_haplotype, event_sequences, alignment_flags);
        for(const auto& v : passed_variants) {
            out_variants.push_back(v);
            v.write_vcf(stderr);
        }
    }
    return out_variants;
}

std::vector<Variant> expand_variants(const AlignmentDB& alignments,
                                     const std::vector<Variant>& candidate_variants,
                                     int region_start,
                                     int region_end,
                                     int calling_span,
                                     uint32_t alignment_flags)
{
    std::vector<Variant> out_variants;

    std::string contig = alignments.get_region_contig();

    for(size_t vi = 0; vi < candidate_variants.size(); ++vi) {
        const Variant& in_variant = candidate_variants[vi];

        // add the variant unmodified
        out_variants.push_back(in_variant);

        // don't do anything with substitutions
        if(in_variant.ref_seq.size() == 1 && in_variant.alt_seq.size() == 1) {
            continue;
        }

        // deletion
        Variant v = candidate_variants[vi];
        v.ref_seq = alignments.get_reference_substring(v.ref_name, v.ref_position, v.ref_position + v.ref_seq.size());
        assert(v.ref_seq != candidate_variants[vi].ref_seq);
        assert(v.ref_seq.substr(0, candidate_variants[vi].ref_seq.size()) == candidate_variants[vi].ref_seq);
        out_variants.push_back(v);

        // insertion
        for(size_t j = 0; j < 4; ++j) {
            v = candidate_variants[vi];
            v.alt_seq.append(1, "ACGT"[j]);
            out_variants.push_back(v);
        }
    }
    return out_variants;
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

Haplotype call_haplotype_from_candidates(const AlignmentDB& alignments,
                                         const std::vector<Variant>& candidate_variants,
                                         int calling_span,
                                         uint32_t alignment_flags)
{
    Haplotype derived_haplotype(alignments.get_region_contig(), alignments.get_region_start(), alignments.get_reference());

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
                alignments.get_event_subsequences(alignments.get_region_contig(), calling_start, calling_end);

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


Haplotype call_variants_for_region(const std::string& contig, int region_start, int region_end)
{
    int calling_span = 10;
    const int BUFFER = 20;
    uint32_t alignment_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

    // load the region, accounting for the buffering
    if(region_start < BUFFER)
        region_start = BUFFER;
    AlignmentDB alignments(opt::reads_file, opt::genome_file, opt::bam_file, opt::event_bam_file, opt::calibrate);

    if(!opt::alternative_model_type.empty()) {
        alignments.set_alternative_model_type(opt::alternative_model_type);
    }

    alignments.load_region(contig, region_start - BUFFER, region_end + BUFFER);

    // Step 1. Discover putative variants across the whole region
    std::vector<Variant> candidate_variants;
    if(opt::candidates_file.empty()) {
        candidate_variants = alignments.get_variants_in_region(contig, region_start, region_end, opt::min_candidate_frequency, 20);
    } else {
        candidate_variants = get_variants_from_vcf(opt::candidates_file, contig, region_start, region_end);
    }

    if(opt::consensus_mode) {

        // generate single-base edits that have a positive haplotype score
        std::vector<Variant> single_base_edits = generate_candidate_single_base_edits(alignments, region_start, region_end, calling_span, alignment_flags);

        // insert these into the candidate set
        candidate_variants.insert(candidate_variants.end(), single_base_edits.begin(), single_base_edits.end());

        // deduplicate variants
        std::set<Variant, VariantKeyComp> dedup_set(candidate_variants.begin(), candidate_variants.end());
        candidate_variants.clear();
        candidate_variants.insert(candidate_variants.end(), dedup_set.begin(), dedup_set.end());
        std::sort(candidate_variants.begin(), candidate_variants.end(), sortByPosition);
    }

    // Step 2. Call variants

    // in consensus mode we iterate until a maximum number of rounds is reached
    // or the variant set converges
    size_t round = 0;
    size_t MAX_ROUNDS = 5;
    Haplotype called_haplotype(alignments.get_region_contig(),
                               alignments.get_region_start(),
                               alignments.get_reference());

    while(opt::consensus_mode && round++ < MAX_ROUNDS) {
        fprintf(stderr, "Round %zu\n", round);
        called_haplotype = call_haplotype_from_candidates(alignments,
                                                          candidate_variants,
                                                          calling_span,
                                                          alignment_flags);

        if(opt::consensus_mode) {
            // Expand the called variant set by adding nearby variants
            std::vector<Variant> called_variants = called_haplotype.get_variants();
            candidate_variants = expand_variants(alignments,
                                                 called_variants,
                                                 region_start,
                                                 region_end,
                                                 calling_span,
                                                 alignment_flags);
        }
    }

    return called_haplotype;
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
            case OPT_CONSENSUS: opt::consensus_mode = 1; break;
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

    if(opt::consensus_mode) {
        fprintf(stderr, "Haplotype: %s\n", haplotype.get_sequence().c_str());
    }

    if(out_fp != stdout) {
        fclose(out_fp);
    }
    return 0;
}
