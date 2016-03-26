//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_methyltest -- test CpG sites for methylation
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
#include "nanopolish_fast5_map.h"
#include "nanopolish_methyltrain.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

//
// Structs
//
struct OutputHandles
{
    FILE* site_writer;
    FILE* read_writer;
    FILE* strand_writer;
};

struct ScoredSite
{
    ScoredSite() 
    { 
        ll_unmethylated[0] = 0;
        ll_unmethylated[1] = 0;
        ll_methylated[0] = 0;
        ll_methylated[1] = 0;
    }

    std::string chromosome;
    int start_position;
    int end_position;
    int n_cpg;
    std::string sequence;

    // scores per strand
    double ll_unmethylated[2];
    double ll_methylated[2];

    //
    static bool sort_by_position(const ScoredSite& a, const ScoredSite& b) { return a.start_position < b.start_position; }

};

//
Alphabet* mtest_alphabet = &gMCpGAlphabet;

//
// Getopt
//
#define SUBPROGRAM "methyltest"

static const char *METHYLTEST_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *METHYLTEST_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Test CpG sites for methylation\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -m, --models-fofn=FILE               read the models from the FOFN\n"
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
void calculate_methylation_for_read(const ModelMap& model_map,
               const Fast5Map& name_map, 
               const faidx_t* fai, 
               const bam_hdr_t* hdr, 
               const bam1_t* record, 
               size_t read_idx,
               const OutputHandles& handles)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    std::string fast5_path = name_map.get_path(read_name);
    SquiggleRead sr(read_name, fast5_path);

    // An output map from reference positions to scored CpG sites
    std::map<int, ScoredSite> site_score_map;

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        std::vector<double> site_scores;
        std::vector<int> site_starts;
        std::vector<int> site_ends;
        std::vector<int> site_count;

        // replace the baked-in pore model with the methylation model
        // (including unmethylated kmers) for this strand
        std::string curr_model = sr.pore_model[strand_idx].name;

        std::string methyl_model = curr_model + ".ecoli_er2925.pcr_MSssI.timp.021216.alphabet_cpg.model";
        auto model_iter = model_map.find(methyl_model);

        if(model_iter != model_map.end()) {
            sr.pore_model[strand_idx].update_states( model_iter->second );
        } else {
            fprintf(stderr, "Error, methylated model %s not found\n", methyl_model.c_str());
            exit(EXIT_FAILURE);
        }
        
        size_t k = sr.pore_model[strand_idx].k;

        // Align in event space using the new model
        EventAlignmentParameters params;
        params.sr = &sr;
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
        params.strand_idx = strand_idx;
        params.read_idx = read_idx;
        params.alphabet = mtest_alphabet;

        std::vector<EventAlignment> alignment_output = align_read_to_ref(params);
        if(alignment_output.empty())
            continue;
        std::string contig = alignment_output.front().ref_name.c_str();
        
        // Convert the EventAlignment to a map between reference positions and events
        std::vector<AlignedPair> event_aligned_pairs;
        for(size_t i = 0; i < alignment_output.size(); ++i) {

            AlignedPair ap = { alignment_output[i].ref_position,
                               alignment_output[i].event_idx };
            event_aligned_pairs.push_back(ap);
        }

        int ref_start_pos = event_aligned_pairs.front().ref_pos;
        int ref_end_pos = event_aligned_pairs.back().ref_pos;

        // Extract the reference sequence for this region
        int fetched_len = 0;
        assert(ref_end_pos >= ref_start_pos);
        std::string ref_seq = get_reference_region_ts(params.fai, contig.c_str(), ref_start_pos, 
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
        size_t curr_idx = 0;
        while(curr_idx < cpg_sites.size()) {
            
            // Find the endpoint of this group of sites
            size_t end_idx = curr_idx + 1;
            while(end_idx < cpg_sites.size()) {
                if(cpg_sites[end_idx] - cpg_sites[end_idx - 1] > min_separation)
                    break;
                end_idx += 1; 
            }

            // the coordinates on the reference substring for this group of sites
            int sub_start_pos = cpg_sites[curr_idx] - min_separation;
            int sub_end_pos = cpg_sites[end_idx - 1] + min_separation;

            if(sub_start_pos > min_separation && cpg_sites[end_idx - 1] - cpg_sites[curr_idx] < 200) {
    
                std::string subseq = ref_seq.substr(sub_start_pos, sub_end_pos - sub_start_pos + 1);
                std::string rc_subseq = mtest_alphabet->reverse_complement(subseq);

                // using the reference-to-event map, look up the event indices for this segment
                AlignedPairRefLBComp lb_comp;
                AlignedPairConstIter start_iter = std::lower_bound(event_aligned_pairs.begin(), event_aligned_pairs.end(),
                                                                   sub_start_pos + ref_start_pos, lb_comp);

                AlignedPairConstIter stop_iter = std::lower_bound(event_aligned_pairs.begin(), event_aligned_pairs.end(),
                                                                  sub_end_pos + ref_start_pos, lb_comp);
                
                // Only process this region if the the read is aligned within the boundaries
                // and the span between the start/end is not unusually short
                if(start_iter != event_aligned_pairs.end() && stop_iter != event_aligned_pairs.end() &&
                    abs(start_iter->read_pos - stop_iter->read_pos) > 10) 
                {
                    
                    uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

                    // Set up event data
                    HMMInputData data;
                    data.read = &sr;
                    data.anchor_index = -1; // unused
                    data.strand = strand_idx;
                    data.rc = alignment_output.front().rc;
                    data.event_start_idx = start_iter->read_pos;
                    data.event_stop_idx = stop_iter->read_pos;
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
                    int start_position = cpg_sites[curr_idx] + ref_start_pos;
                    auto iter = site_score_map.find(start_position);
                    if(iter == site_score_map.end()) {
                        // insert new score into the map
                        ScoredSite ss;
                        ss.chromosome = contig;
                        ss.start_position = start_position;
                        ss.end_position = cpg_sites[end_idx - 1] + ref_start_pos;
                        ss.n_cpg = end_idx - curr_idx;

                        // extract the CpG site(s) with a k-mers worth of surrounding context
                        size_t site_output_start = cpg_sites[curr_idx] - k + 1;
                        size_t site_output_end =  cpg_sites[end_idx - 1] + k;
                        ss.sequence = ref_seq.substr(site_output_start, site_output_end - site_output_start);
                    
                        // insert into the map    
                        iter = site_score_map.insert(std::make_pair(start_position, ss)).first;
                    }
                    
                    // set strand-specific score
                    // upon output below the strand scores will be summed
                    iter->second.ll_unmethylated[strand_idx] = unmethylated_score;
                    iter->second.ll_methylated[strand_idx] = methylated_score;
                }
            }

            curr_idx = end_idx;
        }
    } // for strands
    
    #pragma omp critical(methyltest_write)
    {
        // these variables are sums over all sites within a read
        double ll_ratio_sum_strand[2] = { 0.0f, 0.0f };
        double ll_ratio_sum_both = 0;
        size_t num_positive = 0;

        // write all sites for this read
        for(auto iter = site_score_map.begin(); iter != site_score_map.end(); ++iter) {

            const ScoredSite& ss = iter->second;

            double sum_ll_m = ss.ll_methylated[0] + ss.ll_methylated[1];
            double sum_ll_u = ss.ll_unmethylated[0] + ss.ll_unmethylated[1];

            double diff = sum_ll_m - sum_ll_u;
            num_positive += diff > 0;

            fprintf(handles.site_writer, "%s\t%d\t%d\t", ss.chromosome.c_str(), ss.start_position, ss.end_position);
            fprintf(handles.site_writer, "ReadIdx=%zu;", read_idx);
            fprintf(handles.site_writer, "LogLikMeth=%.2lf;LogLikUnmeth=%.2lf;LogLikRatio=%.2lf;", sum_ll_m, sum_ll_u, diff);
            fprintf(handles.site_writer, "LogLikMethByStrand=%.2lf,%.2lf;", ss.ll_methylated[0], ss.ll_methylated[1]);
            fprintf(handles.site_writer, "LogLikUnmethByStrand=%.2lf,%.2lf;", ss.ll_unmethylated[0], ss.ll_unmethylated[1]);
            fprintf(handles.site_writer, "NumCpGs=%d;Sequence=%s\n", ss.n_cpg, ss.sequence.c_str());

            ll_ratio_sum_strand[0] += ss.ll_methylated[0] - ss.ll_unmethylated[0];
            ll_ratio_sum_strand[1] += ss.ll_methylated[1] - ss.ll_unmethylated[1];
            ll_ratio_sum_both += diff;
        }
        std::string complement_model = sr.pore_model[C_IDX].name;
        fprintf(handles.read_writer, "%s\t%.2lf\t%zu\t%s\tNumPositive=%zu\n", fast5_path.c_str(), ll_ratio_sum_both, site_score_map.size(), complement_model.c_str(), num_positive);
    
        for(size_t si = 0; si < NUM_STRANDS; ++si) {
            std::string model = sr.pore_model[si].name;
            fprintf(handles.strand_writer, "%s\t%.2lf\t%zu\t%s\n", fast5_path.c_str(), ll_ratio_sum_strand[si], site_score_map.size(), model.c_str());
        }
    }
}

void parse_methyltest_options(int argc, char** argv)
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
            case 'v': opt::verbose++; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << METHYLTEST_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << METHYLTEST_VERSION_MESSAGE;
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
    
    if(opt::models_fofn.empty()) {
        std::cerr << SUBPROGRAM ": a --models-fofn file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << METHYLTEST_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int methyltest_main(int argc, char** argv)
{
    parse_methyltest_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    ModelMap models = read_models_fofn(opt::models_fofn, mtest_alphabet);
    
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
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    hts_itr_t* itr;

    // If processing a region of the genome, only emit events aligned to this window
    int clip_start = -1;
    int clip_end = -1;

    if(opt::region.empty()) {
        // TODO: is this valid?
        itr = sam_itr_queryi(bam_idx, HTS_IDX_START, 0, 0);
    } else {

        fprintf(stderr, "Region: %s\n", opt::region.c_str());
        itr = sam_itr_querys(bam_idx, hdr, opt::region.c_str());
        hts_parse_reg(opt::region.c_str(), &clip_start, &clip_end);
    }

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    // Initialize writers
    OutputHandles handles;
    handles.site_writer = fopen(std::string(opt::bam_file + ".methyltest.sites.bed").c_str(), "w");
    handles.read_writer = fopen(std::string(opt::bam_file + ".methyltest.reads.tsv").c_str(), "w");
    handles.strand_writer = fopen(std::string(opt::bam_file + ".methyltest.strand.tsv").c_str(), "w");

    // Write a header to the reads.tsv file
    fprintf(handles.read_writer, "name\tsum_ll_ratio\tn_cpg\tcomplement_model\ttags\n");
    
    // strand header
    fprintf(handles.strand_writer, "name\tsum_ll_ratio\tn_cpg\tmodel\n");


    // Initialize iteration
    std::vector<bam1_t*> records(opt::batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_processed = 0;
    size_t num_records_buffered = 0;
    Progress progress("[methyltest]");

    do {
        assert(num_records_buffered < records.size());
        
        // read a record into the next slot in the buffer
        result = sam_itr_next(bam_fh, itr, records[num_records_buffered]);
        num_records_buffered += result >= 0;

        // realign if we've hit the max buffer size or reached the end of file
        if(num_records_buffered == records.size() || result < 0) {
            
            #pragma omp parallel for
            for(size_t i = 0; i < num_records_buffered; ++i) {
                bam1_t* record = records[i];
                size_t read_idx = num_reads_processed + i;
                if( (record->core.flag & BAM_FUNMAP) == 0) {
                    calculate_methylation_for_read(models, name_map, fai, hdr, record, read_idx, handles);
                }
            }

            num_reads_processed += num_records_buffered;
            num_records_buffered = 0;

        }
    } while(result >= 0);
    
    assert(num_records_buffered == 0);
    progress.end();

    // cleanup records
    for(size_t i = 0; i < records.size(); ++i) {
        bam_destroy1(records[i]);
    }

    // cleanup
    fclose(handles.site_writer);
    fclose(handles.read_writer);
    fclose(handles.strand_writer);

    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    
    return EXIT_SUCCESS;
}

