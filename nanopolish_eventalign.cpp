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
#include "nanopolish_iupac.h"
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
"Align nanopore events to reference k-mers\n"
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

// Modify the aligned_pairs vector to ensure the highest read position
// does not exceed max_kmer
void trim_aligned_pairs_to_kmer(std::vector<AlignedPair>& aligned_pairs, int max_kmer_idx)
{
    int idx = aligned_pairs.size() - 1;
    while(idx >= 0 && aligned_pairs[idx].read_pos > max_kmer_idx)
        idx -= 1;

    if(idx < 0)
        aligned_pairs.clear(); // no valid data
    else
        aligned_pairs.resize(idx + 1);
}

// Modify the aligned_pairs vector to ensure there are no alignments
// outside of the given reference coordinates
void trim_aligned_pairs_to_ref_region(std::vector<AlignedPair>& aligned_pairs, int ref_start, int ref_end)
{
    std::vector<AlignedPair> trimmed;
    for(size_t i = 0; i < aligned_pairs.size(); ++i) {
        if(aligned_pairs[i].ref_pos >= ref_start && 
           aligned_pairs[i].ref_pos <= ref_end) {
            trimmed.push_back(aligned_pairs[i]);
        }
    }
    
    aligned_pairs.swap(trimmed);
}

// Returns the index into the aligned_pairs vector that has the highest ref_pos
// that is not greater than ref_pos_max. It starts the search at pair_idx
int get_end_pair(const std::vector<AlignedPair>& aligned_pairs, int ref_pos_max, int pair_idx)
{
    while(pair_idx < aligned_pairs.size()) {
        if(aligned_pairs[pair_idx].ref_pos > ref_pos_max)
            return pair_idx - 1;
        pair_idx += 1;
    }
    
    return aligned_pairs.size() - 1;
}

struct EventAlignment
{
    // ref data
    std::string ref_name;
    std::string ref_kmer;
    int ref_position;

    // event data
    size_t read_idx;
    int strand_idx;
    int event_idx;
    bool rc;
};

void emit_header(FILE* fp)
{
    fprintf(fp, "%s\t%s\t%s\t%s\t%s\t", "contig", "position", "reference_kmer", "read_index", "strand");
    fprintf(fp, "%s\t%s\t%s\t", "event_index", "event_level_mean", "event_length");
    fprintf(fp, "%s\t%s\t%s\t%s\n", "model_kmer", "model_mean", "model_stdv", "model_name");

}

void emit_event_alignment(FILE* fp,
                          const SquiggleRead& sr,
                          const EventAlignment& ea)
{
    // basic information
    fprintf(fp, "%s\t%d\t%s\t%zu\t%c\t", ea.ref_name.c_str(), ea.ref_position, ea.ref_kmer.c_str(), ea.read_idx, "tc"[ea.strand_idx]);

    // event information
    float event_mean = sr.get_drift_corrected_level(ea.event_idx, ea.strand_idx);
    float event_duration = sr.get_duration(ea.event_idx, ea.strand_idx);
    fprintf(fp, "%d\t%.2lf\t%.3lf\t", ea.event_idx, event_mean, event_duration);

    // model information
    std::string model_kmer = ea.rc ? reverse_complement(ea.ref_kmer) : ea.ref_kmer;

    uint32_t rank = kmer_rank(model_kmer.c_str(), K);
    GaussianParameters model = sr.pore_model[ea.strand_idx].get_scaled_parameters(rank);
    fprintf(fp, "%s\t%.2lf\t%.2lf\t%s\n", model_kmer.c_str(), model.mean, model.stdv, sr.model_name[ea.strand_idx].c_str());
}

// Realign the read in event space
void realign_read(FILE* fp, 
                  const Fast5Map& name_map, 
                  const faidx_t* fai, 
                  const bam_hdr_t* hdr, 
                  const bam1_t* record, 
                  size_t read_idx,
                  int region_start,
                  int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    std::string fast5_path = name_map.get_path(read_name);

    // load read
    SquiggleRead sr(read_name, fast5_path);
    
    if(opt::verbose > 0) {
        fprintf(stderr, "Realigning %s [%zu %zu]\n", 
                read_name.c_str(), sr.events[0].size(), sr.events[1].size());
    }

    // Extract the reference subsequence for the entire alignment
    int fetched_len = 0;
    int ref_offset = record->core.pos;
    std::string ref_name(hdr->target_name[record->core.tid]);

    char* cref_seq = faidx_fetch_seq(fai, ref_name.c_str(), ref_offset, bam_endpos(record), &fetched_len);
    std::string ref_seq(cref_seq);
    free(cref_seq);

    // If the reference sequence contains ambiguity codes, switch them to the lexicographically lowest base
    ref_seq = IUPAC::disambiguate_to_lowest(ref_seq);

    if(ref_offset == 0)
        return;

    // Make a vector of aligned (ref_pos, read_pos) pairs
    std::vector<AlignedPair> aligned_pairs = get_aligned_pairs(record);

    if(region_start != -1 && region_end != -1) {
        trim_aligned_pairs_to_ref_region(aligned_pairs, region_start, region_end);
    }

    // Trim the aligned pairs to be within the range of the maximum kmer index
    int max_kmer_idx = sr.read_sequence.size() - K;
    trim_aligned_pairs_to_kmer(aligned_pairs, max_kmer_idx);

    if(aligned_pairs.empty())
        return;

    bool do_base_rc = bam_is_rev(record);
    bool rc_flags[2] = { do_base_rc, !do_base_rc }; // indexed by strand
    const int align_stride = 100; // approximately how many reference bases to align to at once
    const int output_stride = 50; // approximately how many event alignments to output at once

    for(int strand_idx = 0; strand_idx < 2; ++strand_idx) {

        // get the event range of the read to re-align
        int read_kidx_start = aligned_pairs.front().read_pos;
        int read_kidx_end = aligned_pairs.back().read_pos;
        
        if(do_base_rc) {
            read_kidx_start = sr.flip_k_strand(read_kidx_start);
            read_kidx_end = sr.flip_k_strand(read_kidx_end);
        }
        
        assert(read_kidx_start >= 0);
        assert(read_kidx_end >= 0);

        int first_event = sr.get_closest_event_to(read_kidx_start, strand_idx);
        int last_event = sr.get_closest_event_to(read_kidx_end, strand_idx);
        bool forward = first_event < last_event;

        int last_event_output = -1;
        int curr_start_event = first_event;
        int curr_start_ref = aligned_pairs.front().ref_pos;
        int curr_pair_idx = 0;

        while( (forward && curr_start_event < last_event) ||
               (!forward && curr_start_event > last_event)) {

            // Get the index of the aligned pair approximately align_stride away
            int end_pair_idx = get_end_pair(aligned_pairs, curr_start_ref + align_stride, curr_pair_idx);
        
            int curr_end_ref = aligned_pairs[end_pair_idx].ref_pos;
            int curr_end_read = aligned_pairs[end_pair_idx].read_pos;

            if(do_base_rc) {
                curr_end_read = sr.flip_k_strand(curr_end_read);
            }
            assert(curr_end_read >= 0);

            std::string ref_subseq = ref_seq.substr(curr_start_ref - ref_offset, curr_end_ref - curr_start_ref + 1);
            
            // Nothing to align to
            if(ref_subseq.length() < K)
                break;

            // Set up HMM input
            HMMInputData input;
            input.read = &sr;
            input.anchor_index = 0; // not used here
            input.event_start_idx = curr_start_event;
            input.event_stop_idx = sr.get_closest_event_to(curr_end_read, strand_idx);

            // A limitation of the segment-by-segment alignment is that we can't jump
            // over very large deletions wrt to the reference. The effect of this
            // is that we can get segments that have very few alignable events. We
            // just stop processing them for now
            if(abs(input.event_start_idx - input.event_stop_idx) < 2)
                break;

            input.strand = strand_idx;
            input.event_stride = input.event_start_idx < input.event_stop_idx ? 1 : -1;
            input.rc = rc_flags[strand_idx];
            
            std::vector<AlignmentState> event_alignment = profile_hmm_align(ref_subseq, input);
            //print_alignment("test", 0, 0, ref_subseq, input, event_alignment);
            
            // Output alignment
            size_t num_output = 0;
            size_t event_align_idx = 0;

            // If we aligned to the last event, output everything and stop
            bool last_section = end_pair_idx == aligned_pairs.size() - 1;

            int last_event_output = 0;
            int last_ref_kmer_output = 0;

            for(; event_align_idx < event_alignment.size() && (num_output < output_stride || last_section); event_align_idx++) {

                AlignmentState& as = event_alignment[event_align_idx];
                if(as.state != 'K' && as.event_idx != curr_start_event) {

                    //printf("Outputting event %d aligned to k: %d pos: %d\n", as.event_idx, as.kmer_idx, curr_start_ref + as.kmer_idx);
                    EventAlignment ea;
                    
                    // ref
                    ea.ref_name = ref_name;
                    ea.ref_position = curr_start_ref + as.kmer_idx;
                    ea.ref_kmer = ref_seq.substr(ea.ref_position - ref_offset, K);

                    // event
                    ea.read_idx = read_idx;
                    ea.strand_idx = strand_idx;
                    ea.event_idx = as.event_idx;
                    ea.rc = input.rc;

                    emit_event_alignment(fp, sr, ea);

                    // update
                    last_event_output = as.event_idx;
                    last_ref_kmer_output = curr_start_ref + as.kmer_idx;
                    num_output += 1;
                }
            }

            // Advance the pair iterator to the ref base
            curr_start_event = last_event_output;
            curr_start_ref = last_ref_kmer_output;
            curr_pair_idx = get_end_pair(aligned_pairs, curr_start_ref, curr_pair_idx);
        }
    }
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
        std::cout << "\n" << EVENTALIGN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int eventalign_main(int argc, char** argv)
{
    parse_eventalign_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);

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

    // Initialize iteration
    bam1_t* record = bam_init1();
    
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

    // Write the header
    emit_header(stdout);

    int result;
    size_t read_idx = 0;
    while((result = sam_itr_next(bam_fh, itr, record)) >= 0) {

        realign_read(stdout, name_map, fai, hdr, record, read_idx, clip_start, clip_end);
        read_idx += 1;
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    bam_destroy1(record);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
}
