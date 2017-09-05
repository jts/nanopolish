//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_phase_reads -- phase variants onto reads
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <map>
#include <omp.h>
#include <getopt.h>
#include <cstddef>
#include "htslib/faidx.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_bam_utils.h"
#include "nanopolish_index.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"
#include "logger.hpp"

using namespace std::placeholders;

//
// structs
//

//
// Getopt
//
#define SUBPROGRAM "phase-reads"

static const char *PHASE_READS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2017 Ontario Institute for Cancer Research\n";

static const char *PHASE_READS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa variants.vcf\n"
"Train a duration model\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -w, --window=STR                     only phase reads in the window STR (format: ctg:start-end)\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --progress                       print out a progress message\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string variants_file;
    static std::string region;
    
    static unsigned progress = 0;
    static unsigned num_threads = 1;
    static unsigned batch_size = 128;
    static int min_flanking_sequence = 30;
}

static const char* shortopts = "r:b:g:t:w:v";

enum { OPT_HELP = 1,
       OPT_VERSION,
       OPT_PROGRESS,
       OPT_LOG_LEVEL
     };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "reads",              required_argument, NULL, 'r' },
    { "bam",                required_argument, NULL, 'b' },
    { "genome",             required_argument, NULL, 'g' },
    { "threads",            required_argument, NULL, 't' },
    { "window",             required_argument, NULL, 'w' },
    { "progress",           no_argument,       NULL, OPT_PROGRESS },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "log-level",          required_argument, NULL, OPT_LOG_LEVEL },
    { NULL, 0, NULL, 0 }
};

void parse_phase_reads_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'w': arg >> opt::region; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << PHASE_READS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PHASE_READS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                logger::Logger::set_level_from_option(arg.str());
                break;
        }
    }

    if(argc - optind > 0) {
        opt::variants_file = argv[optind++];
    } else {
        fprintf(stderr, "Error, variants file is missing\n");
        die = true;
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

    if (die) {
        std::cout << "\n" << PHASE_READS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

void phase_single_read(const ReadDB& read_db,
                       const faidx_t* fai,
                       const std::vector<Variant>& variants,
                       samFile* sam_fp,
                       const bam_hdr_t* hdr,
                       const bam1_t* record,
                       size_t read_idx,
                       int region_start,
                       int region_end)
{
    const double MAX_Q_SCORE = 30;
    const double BAM_Q_OFFSET = 0;
    int tid = omp_get_thread_num();
    uint32_t alignment_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;
    
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);

    // load read
    SquiggleRead sr(read_name, read_db);
    
    std::string ref_name = hdr->target_name[record->core.tid];
    int alignment_start_pos = record->core.pos;
    int alignment_end_pos = bam_endpos(record);
    
    // Search the variant collection for the index of the first/last variants to phase
    Variant lower_search;
    lower_search.ref_name = ref_name;
    lower_search.ref_position = alignment_start_pos;
    auto lower_iter = std::lower_bound(variants.begin(), variants.end(), lower_search, sortByPosition);

    Variant upper_search;
    upper_search.ref_name = ref_name;
    upper_search.ref_position = alignment_end_pos;
    auto upper_iter = std::upper_bound(variants.begin(), variants.end(), upper_search, sortByPosition);

    fprintf(stderr, "%s %s:%u-%u %zu\n", read_name.c_str(), ref_name.c_str(), alignment_start_pos, alignment_end_pos, upper_iter - lower_iter);

    // no variants to phase?
    if(lower_iter == variants.end()) {
        return;
    }

    int fetched_len;
    std::string reference_seq = get_reference_region_ts(fai, 
                                                        ref_name.c_str(), 
                                                        alignment_start_pos, 
                                                        alignment_end_pos, 
                                                        &fetched_len);

    std::string read_outseq = reference_seq;
    std::string read_outqual(reference_seq.length(), MAX_Q_SCORE + BAM_Q_OFFSET);

    Haplotype reference_haplotype(ref_name, alignment_start_pos, reference_seq);
    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {

        // skip if 1D reads and this is the wrong strand
        if(!sr.has_events_for_strand(strand_idx)) {
            continue;
        }

        // only phase using template strand
        if(strand_idx != 0) {
            continue;
        }

        SequenceAlignmentRecord seq_align_record(record);
        EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);

        // 
        for(; lower_iter < upper_iter; ++lower_iter) {

            const Variant& v = *lower_iter;

            if(!v.is_snp()) {
                continue;
            }

            int calling_start = v.ref_position - opt::min_flanking_sequence;
            int calling_end = v.ref_position + opt::min_flanking_sequence;

            HMMInputData data;
            data.read = event_align_record.sr;
            data.strand = event_align_record.strand;
            data.rc = event_align_record.rc;
            data.event_stride = event_align_record.stride;
            
            int e1,e2;
            bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events, 
                                                            calling_start,
                                                            calling_end,
                                                            e1, 
                                                            e2);

            // The events of this read do not span the calling window, skip
            if(!bounded || fabs(e2 - e1) / (calling_start - calling_end) > MAX_EVENT_TO_BP_RATIO) {
                continue;
            }

            data.event_start_idx = e1;
            data.event_stop_idx = e2;

            Haplotype calling_haplotype =
                reference_haplotype.substr_by_reference(calling_start, calling_end);
        
            double ref_score = profile_hmm_score(calling_haplotype.get_sequence(), data, alignment_flags);
            bool good_haplotype = calling_haplotype.apply_variant(v);
            if(good_haplotype) {
                double alt_score = profile_hmm_score(calling_haplotype.get_sequence(), data, alignment_flags);
                double log_sum = add_logs(alt_score, ref_score);
                double log_p_ref = ref_score - log_sum;
                double log_p_alt = alt_score - log_sum;
                char call;
                double log_p_wrong;
                if(alt_score > ref_score) {
                     call = v.alt_seq[0];
                     log_p_wrong = log_p_ref;
                } else {
                    call = v.ref_seq[0];
                    log_p_wrong = log_p_alt;
                }

                double q_score = -10 * log_p_wrong / log(10);
                q_score = std::min(MAX_Q_SCORE, q_score);
                char q_char = (int)q_score + 33;
                //fprintf(stderr, "\t%s score: %.2lf %.2lf %c p_wrong: %.3lf Q: %d QC: %c\n", v.key().c_str(), ref_score, alt_score, call, log_p_wrong, (int)q_score, q_char);

                int out_position = v.ref_position - alignment_start_pos;
                assert(read_outseq[out_position] == v.ref_seq[0]);
                read_outseq[out_position] = call;
                read_outqual[out_position] = q_char;
            }
        }

        // Construct the output bam record
        bam1_t* out_record = bam_init1();
        
        // basic stats
        out_record->core.tid = record->core.tid;
        out_record->core.pos = alignment_start_pos;
        out_record->core.qual = record->core.qual;
        out_record->core.flag = record->core.flag & BAM_FSUPPLEMENTARY;
        
        // no read pairs
        out_record->core.mtid = -1;
        out_record->core.mpos = -1;
        out_record->core.isize = 0;

        std::vector<uint32_t> cigar;
        uint32_t cigar_op = read_outseq.size() << BAM_CIGAR_SHIFT | BAM_CMATCH;
        cigar.push_back(cigar_op);
        write_bam_vardata(out_record, read_name, cigar, read_outseq, read_outqual);

        #pragma omp critical
        {
            sam_write1(sam_fp, hdr, out_record);
        }
        bam_destroy1(out_record); // automatically frees malloc'd segment

    } // for strand
}

int phase_reads_main(int argc, char** argv)
{
    parse_phase_reads_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());
  
    std::vector<Variant> variants;  
    if(!opt::region.empty()) {
        std::string contig;
        int start_base;
        int end_base;
        parse_region_string(opt::region, contig, start_base, end_base);

        // Read the variants for this region
        variants = read_variants_for_region(opt::variants_file, contig, start_base, end_base);
    } else {
         variants = read_variants_from_file(opt::variants_file);
    }

    // Sort variants by reference coordinate
    std::sort(variants.begin(), variants.end(), sortByPosition);
    
    // remove hom reference
    auto new_end = std::remove_if(variants.begin(), variants.end(), [](Variant v) { return v.genotype == "0/0"; });
    variants.erase( new_end, variants.end());
    
    samFile* sam_out = sam_open("-", "w");

    // the BamProcessor framework calls the input function with the 
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(phase_single_read, std::ref(read_db), std::ref(fai), std::ref(variants), sam_out, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    
    // Copy the bam header to std
    sam_hdr_write(sam_out, processor.get_bam_header());
    
    processor.parallel_run(f);
    
    fai_destroy(fai);
    sam_close(sam_out);
    
    return EXIT_SUCCESS;
}
