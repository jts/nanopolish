//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_call_methylation -- identify methylated bases
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
#include <fstream>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_read_db.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

using namespace std::placeholders;

//
// Structs
//
struct OutputHandles
{
    FILE* site_writer;
};

struct ScoredSite
{
    ScoredSite() 
    { 
        ll_unmethylated[0] = 0;
        ll_unmethylated[1] = 0;
        ll_methylated[0] = 0;
        ll_methylated[1] = 0;
        strands_scored = 0;
    }

    std::string chromosome;
    int start_position;
    int end_position;
    int n_cpg;
    std::string sequence;

    // scores per strand
    double ll_unmethylated[2];
    double ll_methylated[2];
    int strands_scored;

    //
    static bool sort_by_position(const ScoredSite& a, const ScoredSite& b) { return a.start_position < b.start_position; }

};

//
Alphabet* mtest_alphabet = &gMCpGAlphabet;

//
// Getopt
//
#define SUBPROGRAM "call-methylation"

static const char *CALL_METHYLATION_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CALL_METHYLATION_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Classify nucleotides as methylated or not.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --progress                       print out a progress message\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string models_fofn;
    static std::string region;
    static std::string cpg_methylation_model_type = "reftrained";
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:w:m:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "reads",            required_argument, NULL, 'r' },
    { "bam",              required_argument, NULL, 'b' },
    { "genome",           required_argument, NULL, 'g' },
    { "window",           required_argument, NULL, 'w' },
    { "threads",          required_argument, NULL, 't' },
    { "models-fofn",      required_argument, NULL, 'm' },
    { "progress",         no_argument,       NULL, OPT_PROGRESS },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Test CpG sites in this read for methylation
void calculate_methylation_for_read(const OutputHandles& handles,
                                    const ReadDB& read_db,
                                    const faidx_t* fai,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    SquiggleRead sr(read_name, read_db);

    // An output map from reference positions to scored CpG sites
    std::map<int, ScoredSite> site_score_map;

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        if(!sr.has_events_for_strand(strand_idx)) {
            continue;
        }

        size_t k = sr.get_model_k(strand_idx);

        // check if there is a cpg model for this strand
        if(!PoreModelSet::has_model(sr.get_model_kit_name(strand_idx),
                                    "cpg",
                                    sr.get_model_strand_name(strand_idx),
                                    k))
        {
            continue;
        }

        // Build the event-to-reference map for this read from the bam record
        SequenceAlignmentRecord seq_align_record(record);
        EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);

        std::vector<double> site_scores;
        std::vector<int> site_starts;
        std::vector<int> site_ends;
        std::vector<int> site_count;


        std::string contig = hdr->target_name[record->core.tid];
        int ref_start_pos = record->core.pos;
        int ref_end_pos =  bam_endpos(record);

        // Extract the reference sequence for this region
        int fetched_len = 0;
        assert(ref_end_pos >= ref_start_pos);
        std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos, 
                                                      ref_end_pos, &fetched_len);
        
        // Remove non-ACGT bases from this reference segment
        ref_seq = gDNAAlphabet.disambiguate(ref_seq);

        // Scan the sequence for CpGs
        std::vector<int> cpg_sites;
        assert(ref_seq.size() != 0);
        for(size_t i = 0; i < ref_seq.size() - 1; ++i) {
            if(ref_seq[i] == 'C' && ref_seq[i+1] == 'G') {
                cpg_sites.push_back(i);
            }
        }

        // Batch the CpGs together into groups that are separated by some minimum distance
        int min_separation = 10;
        std::vector<std::pair<int, int>> groups;
        
        size_t curr_idx = 0;
        while(curr_idx < cpg_sites.size()) {
            // Find the endpoint of this group of sites
            size_t end_idx = curr_idx + 1;
            while(end_idx < cpg_sites.size()) {
                if(cpg_sites[end_idx] - cpg_sites[end_idx - 1] > min_separation)
                    break;
                end_idx += 1; 
            }
            groups.push_back(std::make_pair(curr_idx, end_idx));
            curr_idx = end_idx;
        }

        for(size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {

            size_t start_idx = groups[group_idx].first;
            size_t end_idx = groups[group_idx].second;

            // the coordinates on the reference substring for this group of sites
            int sub_start_pos = cpg_sites[start_idx] - min_separation;
            int sub_end_pos = cpg_sites[end_idx - 1] + min_separation;
            int span = cpg_sites[end_idx - 1] - cpg_sites[start_idx];

            // skip if too close to the start of the read alignment or
            // if the reference range is too large to efficiently call
            if(sub_start_pos <= min_separation || span > 200) {
                continue;
            }

            std::string subseq = ref_seq.substr(sub_start_pos, sub_end_pos - sub_start_pos + 1);
            std::string rc_subseq = mtest_alphabet->reverse_complement(subseq);

            int calling_start = sub_start_pos + ref_start_pos;
            int calling_end = sub_end_pos + ref_start_pos;

            // using the reference-to-event map, look up the event indices for this segment
            int e1,e2;
            bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events, 
                                                            calling_start,
                                                            calling_end,
                                                            e1, 
                                                            e2);

            double ratio = fabs(e2 - e1) / (calling_start - calling_end); 
            
            // Only process this region if the the read is aligned within the boundaries
            // and the span between the start/end is not unusually short
            if(!bounded || abs(e2 - e1) <= 10 || ratio > MAX_EVENT_TO_BP_RATIO) {
                continue;
            }

            uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

            // Set up event data
            HMMInputData data;
            data.read = &sr;
            data.pore_model = sr.get_model(strand_idx, "cpg");
            data.strand = strand_idx;
            data.rc = event_align_record.rc;
            data.event_start_idx = e1;
            data.event_stop_idx = e2;
            data.event_stride = data.event_start_idx <= data.event_stop_idx ? 1 : -1;
         
            // Calculate the likelihood of the unmethylated sequence
            HMMInputSequence unmethylated(subseq, rc_subseq, mtest_alphabet);
            double unmethylated_score = profile_hmm_score(unmethylated, data, hmm_flags);

            // Methylate all CpGs in the sequence and score again
            std::string mcpg_subseq = mtest_alphabet->methylate(subseq);
            std::string rc_mcpg_subseq = mtest_alphabet->reverse_complement(mcpg_subseq);
            
            // Calculate the likelihood of the methylated sequence
            HMMInputSequence methylated(mcpg_subseq, rc_mcpg_subseq, mtest_alphabet);
            double methylated_score = profile_hmm_score(methylated, data, hmm_flags);

            // Aggregate score
            int start_position = cpg_sites[start_idx] + ref_start_pos;
            auto iter = site_score_map.find(start_position);
            if(iter == site_score_map.end()) {
                // insert new score into the map
                ScoredSite ss;
                ss.chromosome = contig;
                ss.start_position = start_position;
                ss.end_position = cpg_sites[end_idx - 1] + ref_start_pos;
                ss.n_cpg = end_idx - start_idx;

                // extract the CpG site(s) with a k-mers worth of surrounding context
                size_t site_output_start = cpg_sites[start_idx] - k + 1;
                size_t site_output_end =  cpg_sites[end_idx - 1] + k;
                ss.sequence = ref_seq.substr(site_output_start, site_output_end - site_output_start);
            
                // insert into the map    
                iter = site_score_map.insert(std::make_pair(start_position, ss)).first;
            }
            
            // set strand-specific score
            // upon output below the strand scores will be summed
            iter->second.ll_unmethylated[strand_idx] = unmethylated_score;
            iter->second.ll_methylated[strand_idx] = methylated_score;
            iter->second.strands_scored += 1;
        } // for group
    } // for strands
    
    #pragma omp critical(call_methylation_write)
    {
        // write all sites for this read
        for(auto iter = site_score_map.begin(); iter != site_score_map.end(); ++iter) {

            const ScoredSite& ss = iter->second;
            double sum_ll_m = ss.ll_methylated[0] + ss.ll_methylated[1];
            double sum_ll_u = ss.ll_unmethylated[0] + ss.ll_unmethylated[1];
            double diff = sum_ll_m - sum_ll_u;

            fprintf(handles.site_writer, "%s\t%d\t%d\t", ss.chromosome.c_str(), ss.start_position, ss.end_position);
            fprintf(handles.site_writer, "%s\t%.2lf\t", sr.read_name.c_str(), diff);
            fprintf(handles.site_writer, "%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
            fprintf(handles.site_writer, "%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());
        }
    }
}

void parse_call_methylation_options(int argc, char** argv)
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
            case 'm': arg >> opt::models_fofn; break;
            case 'w': arg >> opt::region; break;
            case 'v': opt::verbose++; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << CALL_METHYLATION_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CALL_METHYLATION_VERSION_MESSAGE;
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

    if(!opt::models_fofn.empty()) {
        // initialize the model set from the fofn
        PoreModelSet::initialize(opt::models_fofn);
    }

    if (die)
    {
        std::cout << "\n" << CALL_METHYLATION_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int call_methylation_main(int argc, char** argv)
{
    parse_call_methylation_options(argc, argv);
    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    // Initialize writers
    OutputHandles handles;
    handles.site_writer = stdout;
    
    // Write header
    fprintf(handles.site_writer, "chromosome\tstart\tend\tread_name\t"
                                 "log_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\t"
                                 "num_calling_strands\tnum_cpgs\tsequence\n");

    // the BamProcessor framework calls the input function with the 
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(calculate_methylation_for_read, std::ref(handles), std::ref(read_db), fai, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.parallel_run(f);

    // cleanup
    if(handles.site_writer != stdout) {
        fclose(handles.site_writer);
    }

    fai_destroy(fai);

    return EXIT_SUCCESS;
}

