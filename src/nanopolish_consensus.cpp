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
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "nanopolish_hmm_input_sequence.h"
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
#define SUBPROGRAM "consensus"

static const char *CONSENSUS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CONSENSUS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] -w contig:start-end --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Compute a new consensus sequence for an assembly using a signal-level HMM\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -o, --outfile=FILE                   write result to FILE [default: stdout]\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --models-fofn=FILE               read alternative k-mer models from FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string output_file;
    static std::string window;
    static std::string models_fofn;
    static int show_progress = 0;
    static int num_threads = 1;
}

static const char* shortopts = "r:b:g:t:w:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VCF, OPT_PROGRESS, OPT_MODELS_FOFN };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "reads",       required_argument, NULL, 'r' },
    { "bam",         required_argument, NULL, 'b' },
    { "genome",      required_argument, NULL, 'g' },
    { "window",      required_argument, NULL, 'w' },
    { "outfile",     required_argument, NULL, 'o' },
    { "threads",     required_argument, NULL, 't' },
    { "models-fofn", required_argument, NULL, OPT_MODELS_FOFN },
    { "progress",    no_argument,       NULL, OPT_PROGRESS },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

std::vector<HMMInputData> get_input_for_columns(HMMRealignmentInput& window,
                                                const HMMAnchoredColumn& start_column,
                                                const HMMAnchoredColumn& end_column)
{
    assert(start_column.anchors.size() == end_column.anchors.size());

    std::vector<HMMInputData> input;
    for(uint32_t rsi = 0; rsi < start_column.anchors.size(); ++rsi) {

        HMMStrandAnchor start_sa = start_column.anchors[rsi];
        HMMStrandAnchor end_sa = end_column.anchors[rsi];

        // sanity checks
        // This read strand does not have events at both anchors
        if(start_sa.event_idx == -1 || end_sa.event_idx == -1)
            continue;

        // require a minimum number of events
        uint32_t n_events = abs(start_sa.event_idx - end_sa.event_idx);
        if(n_events < 20 || n_events > 500)
            continue;

        if(start_sa.rc != end_sa.rc)
            continue;


        HMMInputData data;

        uint32_t read_idx = rsi / 2;
        assert(read_idx < window.reads.size());
        data.anchor_index = rsi;
        data.read = window.reads[read_idx].get();
        data.strand = rsi % 2;
        data.event_start_idx = start_sa.event_idx;
        data.event_stop_idx = end_sa.event_idx;
        if(data.event_start_idx < data.event_stop_idx)
            data.event_stride = 1;
        else
            data.event_stride = -1;
        data.rc = start_sa.rc;

        input.push_back(data);
    }
    return input;
}

// Handy wrappers for scoring/debugging functions
// The consensus algorithms call into these so we can switch
// scoring functions without writing a bunch of code
double score_sequence(const std::string& sequence, const HMMInputData& data)
{
    return profile_hmm_score(sequence, data);
}

void update_training_with_segment(const HMMInputSequence& sequence, const HMMInputData& data)
{
    std::vector<HMMAlignmentState> alignment = profile_hmm_align(sequence, data);
    data.read->parameters[data.strand].add_training_from_alignment(sequence, data, alignment);
}

struct PathCons
{
    // default constructor
    PathCons(const std::string& s) : path(s), score(0.0f), sum_rank(0) {}

    std::string path;
    
    double score;
    size_t sum_rank;
    size_t num_improved;
    size_t num_scored;
    
    std::string mutdesc;
    
};
typedef std::vector<PathCons> PathConsVector;

bool sortPathConsScoreDesc(const PathCons& a, const PathCons& b)
{
    return a.score > b.score;
}

bool sortPathConsRankAsc(const PathCons& a, const PathCons& b)
{
    return a.sum_rank < b.sum_rank;
}

bool sortPathConsRankDesc(const PathCons& a, const PathCons& b)
{
    return a.sum_rank > b.sum_rank;
}

struct IndexedPathScore
{
    double score;
    uint32_t path_index;
};

bool sortIndexedPathScoreDesc(const IndexedPathScore& a, const IndexedPathScore& b)
{
    return a.score > b.score;
}

// This scores each path using the HMM and 
// sorts the paths into ascending order by score
void score_paths(PathConsVector& paths, const std::vector<HMMInputData>& input)
{
    PROFILE_FUNC("score_paths")
    size_t CULL_RATE = 5;
    double CULL_MIN_SCORE = -30.0f;
    double CULL_MIN_IMPROVED_FRACTION = 0.2f;

    // cache the initial sequence
    std::string first = paths[0].path;
    
    PathConsVector dedup_paths;

    // initialize and deduplicate paths to avoid redundant computation
    std::set<std::string> path_string_set;
    for(size_t pi = 0; pi < paths.size(); ++pi) {

        if(path_string_set.find(paths[pi].path) == path_string_set.end()) {
            paths[pi].score = 0;
            paths[pi].sum_rank = 0;
            paths[pi].num_improved = 0;
            paths[pi].num_scored = 0;
            dedup_paths.push_back(paths[pi]);
            path_string_set.insert(paths[pi].path);
        }
    }
    paths.clear();
    paths.swap(dedup_paths);
    

    // Score all reads
    for(uint32_t ri = 0; ri < input.size(); ++ri) {

        if(opt::verbose > 2) {
            fprintf(stderr, "Scoring %d\n", ri);
        }

        //const HMMInputData& data = input[ri];
        std::vector<IndexedPathScore> result(paths.size());

        // Score all paths
        #pragma omp parallel for
        for(size_t pi = 0; pi < paths.size(); ++pi) {
            double curr = score_sequence(paths[pi].path, input[ri]);
            result[pi].score = curr;
            result[pi].path_index = pi;
        }

        // Save score of first path
        double first_path_score = result[0].score;

        // Sort result by score
        std::stable_sort(result.begin(), result.end(), sortIndexedPathScoreDesc);

        for(size_t pri = 0; pri < result.size(); ++pri) {
            size_t pi = result[pri].path_index;

            paths[pi].score += (result[pri].score - first_path_score);
            uint32_t rank_score = pri;
            paths[pi].sum_rank += rank_score;
            paths[pi].num_improved += (result[pri].score > first_path_score);
            paths[pi].num_scored += 1;
        }

        // Cull paths
        if(ri > 0 && ri % CULL_RATE == 0) {
            PathConsVector retained_paths;
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                
                // We keep a path if any of these conditions are met:
                //  1) it is the original unmodified sequence
                //  2) its score is greater than CULL_MIN_SCORE
                //  3) the fraction of reads that score better on this
                //     path compared to the original sequence is greater
                //     than CULL_MIN_IMPROVED_FRACTION
                double f = (double)paths[pi].num_improved / (double)paths[pi].num_scored;
                if(pi == 0 || paths[pi].score > CULL_MIN_SCORE || f >= CULL_MIN_IMPROVED_FRACTION) {
                    retained_paths.push_back(paths[pi]);
                }
            }
            paths.swap(retained_paths);
        }
    }

    // select new sequence
    //std::stable_sort(paths.begin(), paths.end(), sortPathConsRankAsc);
    std::stable_sort(paths.begin(), paths.end(), sortPathConsScoreDesc);

#if DEBUG_PATH_SELECTION
    double MIN_FIT = INFINITY;
    for(size_t pi = 0; pi < paths.size(); ++pi) {

        // Calculate the length of the matching prefix with the initial sequence
        const std::string& s = paths[pi].path;

        char initial = s == first ? 'I' : ' ';

        fprintf(stderr, "%zu\t%s\t%.1lf\t%zu %c %s", pi, paths[pi].path.c_str(), paths[pi].score, paths[pi].sum_rank, initial, paths[pi].mutdesc.c_str());
        // If this is the truth path or the best path, show the scores for all reads
        if(pi <= 1 || initial == 'I') {
            for(uint32_t ri = 0; ri < input.size(); ++ri) {
                const HMMInputData& data = input[ri];
                double curr = score_sequence(paths[pi].path, input[ri]);
                fprintf(stderr, "%.2lf ", curr);
            }
        }
        fprintf(stderr, "\n");
    }
#endif

}

void extend_paths(PathConsVector& paths, int maxk = 2)
{
    // Insert all possible extensions into the path sequence
    // for k in 1 to maxk
    PathConsVector new_paths;

    for(int k = 1; k <= maxk; ++k) {

        for(unsigned pi = 0; pi < paths.size(); ++pi) {
    
            std::string first(k, 'A');
            std::string extension = first;

            do {
                std::string current = paths[pi].path;
                std::string ns = current.insert(current.size() - 5, extension);
                PathCons ps(ns);
                new_paths.push_back(ps);
                gDNAAlphabet.lexicographic_next(extension);
            } while(extension != first);
        }
    }

    paths.swap(new_paths);
}

PathConsVector generate_mutations(const std::string& sequence, const uint32_t k)
{
    PathConsVector mutations;

    // Add the unmutated sequence
    {
        PathCons pc(sequence);
        mutations.push_back(pc);
    }

    // Mutate every base except for in the first/last k-mer
    for(size_t si = k; si < sequence.size() - k; ++si) {
        
        // All subs
        for(size_t bi = 0; bi < 4; bi++) {
            char b = "ACGT"[bi];
            if(sequence[si] == b)
                continue;
            PathCons pc(sequence);
            pc.path[si] = b;
            std::stringstream ss;
            ss << "sub-" << si << "-" << b;
            pc.mutdesc = ss.str();
            mutations.push_back(pc);
        }

        // 1bp del at this position
        {
            PathCons pc(sequence);
            pc.path.erase(si, 1);
            
            std::stringstream ss;
            ss << "del-" << si;
            pc.mutdesc = ss.str();
            
            mutations.push_back(pc);
        }

        // All 1bp ins before this position
        for(size_t bi = 0; bi < 4; bi++) {
            char b = "ACGT"[bi];
            PathCons pc(sequence);
            pc.path.insert(si, 1, b);
            
            std::stringstream ss;
            ss << "ins-" << si << "-" << b;
            pc.mutdesc = ss.str();
            
            mutations.push_back(pc);
        }
    }

    return mutations;
}

// Run the mutation algorithm to generate an improved consensus sequence
std::string run_mutation(const std::string& base, const std::vector<HMMInputData>& input)
{
    PROFILE_FUNC("run_mutation")
    std::string result = base;

    // assume models for all the reads have the same k
    assert(!input.empty());
    const uint32_t k = input[0].read->pore_model[input[0].strand].k;

    int iteration = 0;
    while(iteration++ < 10) {

        // Generate possible sequences
        PathConsVector paths = generate_mutations(result, k);

        // score them in the HMM
        score_paths(paths, input);

        // check if no improvement was made
        if(paths[0].path == result)
            break;
        result = paths[0].path;
    }

    return result;
}

void generate_alt_paths(PathConsVector& paths, const std::string& base, const std::vector<std::string>& alts, 
                        const uint32_t k)
{
    // Generate alternatives
    for(uint32_t ai = 0; ai < alts.size(); ++ai) {
        const std::string& alt = alts[ai];

        if(alt.size() < k)
            continue;

        kLCSResult result = kLCS(base, alt, k);

#ifdef DEBUG_ALT_GENERATION
        printf("Match to alt %s\n", alt.c_str());
        for(size_t mi = 0; mi < result.size(); ++mi) {
            std::string extend = "";
            if(mi < result.size() - 1 && result[mi].j + 1 != result[mi + 1].j) {
                extend = alt.substr(result[mi].j, result[mi + 1].j - result[mi].j + k);
            }
            printf("\t%zu %zu %s %s\n", result[mi].i, result[mi].j, base.substr(result[mi].i, k).c_str(), extend.c_str());
        }
#endif

        uint32_t match_idx = 0;
        while(match_idx < result.size()) {
            uint32_t last_idx = result.size() - 1;

            // advance the match to the next point of divergence
            while(match_idx != last_idx && 
                  result[match_idx].i == result[match_idx + 1].i - 1 &&
                  result[match_idx].j == result[match_idx + 1].j - 1) {
                match_idx++;
            }
            // no more divergences to process
            if(match_idx == last_idx)
                break;

            uint32_t bl = result[match_idx + 1].i - result[match_idx].i;
            uint32_t rl = result[match_idx + 1].j - result[match_idx].j;

            std::string base_subseq = base.substr(result[match_idx].i, bl);
            std::string alt_subseq = alt.substr(result[match_idx].j, rl);
            
            // Perform the splice
            PathCons new_path(base);
            new_path.path.replace(result[match_idx].i, bl, alt_subseq);
            paths.push_back(new_path);
            
            match_idx += 1;
        }
    }
}

// Run the block substitution algorithm to generate an improved consensus sequence
std::string run_block_substitution(const std::string& base,
                                   const std::vector<HMMInputData>& input,
                                   const std::vector<std::string>& alts)
{
    std::string result = base;

    // assume models for all the reads have the same k
    assert(!input.empty());
    const uint32_t k = input[0].read->pore_model[input[0].strand].k;

    uint32_t max_rounds = 6;
    uint32_t round = 0;
    while(round++ < max_rounds) {
        
        PathConsVector paths;
        PathCons initial_path(result);
        paths.push_back(initial_path);
        
        generate_alt_paths(paths, result, alts, k);
        score_paths(paths, input);

        if(paths[0].path == result)
            break;
        result = paths[0].path;
    }
    return result;
}

//
// Outlier filtering
//
void filter_outlier_data(std::vector<HMMInputData>& input, const std::string& sequence)
{
    std::vector<HMMInputData> out_rs;
    for(uint32_t ri = 0; ri < input.size(); ++ri) {
        const HMMInputData& rs = input[ri];

        double curr = score_sequence(sequence, rs);
        double n_events = rs.event_start_idx - rs.event_stop_idx + 1.0f;
        double lp_per_event = curr / n_events;

        if(opt::verbose >= 1) {
            fprintf(stderr, "OUTLIER_FILTER %d %.2lf %.2lf %.2lf\n", ri, curr, n_events, lp_per_event);
        }
        
        // R9 thresholds
        double threshold = model_stdv() ? 8.0f : 4.0f; // TODO: check
        if(fabs(lp_per_event) < threshold) {
            out_rs.push_back(rs);
        }
    }
    input.swap(out_rs);
}

std::string join_sequences_at_kmer(const std::string& a, const std::string& b, const uint32_t k)
{
    // this is a special case to make the calling code cleaner
    if(a.empty())
        return b;

    // These sequences must have a k-mer match at the start/end
    std::string a_last_kmer = a.substr(a.size() - k);
    std::string b_last_kmer = b.substr(0, k);
    assert(a_last_kmer == b_last_kmer);
    return a + b.substr(k);
}

void run_splice_segment(HMMRealignmentInput& window, uint32_t segment_id, const uint32_t k)
{
    // The structure of the data looks like this:

    // --------------------------------------------------------
    // S                       M                              E
    // where is the start column, M is the middle column and E
    // is the end column. We want to call a new consensus from S
    // to E. We do this by generating the base sequence from S to E
    // and then applying all of the alternatives indicated by the
    // start and middle column. We score these alternatives using
    // the read strands spanning from S to E. After a new consensus
    // has been selected, we re-calculate the alignments of events to
    // the middle anchor.

    // Get the segments
    assert(segment_id + 2 < window.anchored_columns.size());
    HMMAnchoredColumn& start_column = window.anchored_columns[segment_id];
    HMMAnchoredColumn& middle_column = window.anchored_columns[segment_id + 1];
    HMMAnchoredColumn& end_column = window.anchored_columns[segment_id + 2];

    std::string s_m_base = start_column.base_sequence;
    std::string m_e_base = middle_column.base_sequence;

    // The collection of alternative sequences
    std::vector<std::string> alts;

    for(uint32_t ai = 0; ai < start_column.alt_sequences.size(); ++ai) {
        alts.push_back(start_column.alt_sequences[ai]);
    }

    // set up the input data for the HMM
    std::vector<HMMInputData> data = get_input_for_columns(window, start_column, end_column);

    if(opt::verbose > 0) {
        fprintf(stderr, "correcting segment %u with %zu reads\n", segment_id, data.size());
    }
    
    // The current consensus sequence
    std::string original = join_sequences_at_kmer(s_m_base, m_e_base, k);
    std::string base = original;

    // filter out poor quality reads
    filter_outlier_data(data, base);

    // Only attempt correction if there are any reads here
    if(!data.empty()) {
        
        std::string bs_result = run_block_substitution(base, data, alts);
        std::string mut_result = run_mutation(bs_result, data);
        base = mut_result;
    }

    if(opt::verbose > 0) {
        fprintf(stderr, "ORIGINAL[%d] %s\n", segment_id, original.c_str());
        fprintf(stderr, "RESULT[%d]   %s\n", segment_id, base.c_str());
    }
        
    // Update the sequences for the start and middle segments
    // by cutting the new consensus in the middle
    // We maintain the k-mer match invariant by requiring the
    // sequences to overlap by k-bp
    assert(base.length() >= k);
    uint32_t midpoint_kmer = (base.length() - k + 1) / 2;

    std::string s_m_fixed = base.substr(0, midpoint_kmer + k);
    std::string m_e_fixed = base.substr(midpoint_kmer);

    assert(s_m_fixed.substr(s_m_fixed.size() - k) == m_e_fixed.substr(0, k));

    start_column.base_sequence = s_m_fixed;
    middle_column.base_sequence = m_e_fixed;

    // Update the event indices in the first column to match 
    for(uint32_t ri = 0; ri < data.size(); ++ri) {

        // Realign to the consensus sequence
        std::vector<HMMAlignmentState> decodes = profile_hmm_align(base, data[ri]);

        // Get the closest event aligned to the target kmer
        int32_t min_k_dist = base.length();
        uint32_t event_idx = 0;
        for(uint32_t di = 0; di < decodes.size(); ++di) {
            int32_t dist = decodes[di].kmer_idx - midpoint_kmer;
            if(dist <= min_k_dist) {
                min_k_dist = dist;
                event_idx = decodes[di].event_idx;
            }
        }

        middle_column.anchors[data[ri].anchor_index].event_idx = event_idx;
    }
}

// update the training data on the current segment
void train_segment(HMMRealignmentInput& window, uint32_t segment_id)
{
    // Get the segments
    assert(segment_id + 2 < window.anchored_columns.size());
    HMMAnchoredColumn& start_column = window.anchored_columns[segment_id];
    HMMAnchoredColumn& middle_column = window.anchored_columns[segment_id + 1];
    HMMAnchoredColumn& end_column = window.anchored_columns[segment_id + 2];

    std::string s_m_base = start_column.base_sequence;
    std::string m_e_base = middle_column.base_sequence;

    // Set up the the input data for the HMM
    std::vector<HMMInputData> input = get_input_for_columns(window, start_column, end_column);

    // no training can be performed if there are no reads for this segment
    if(input.empty()) {
        return;
    }

    // assume models for all the reads have the same k
    const uint32_t k = input[0].read->pore_model[input[0].strand].k;

    std::string segment_sequence = join_sequences_at_kmer(s_m_base, m_e_base, k);
     
    for(uint32_t ri = 0; ri < input.size(); ++ri) {
        std::vector<HMMAlignmentState> decodes = profile_hmm_align(segment_sequence, input[ri]);
        update_training_with_segment(segment_sequence, input[ri]);
    }
}

void train(HMMRealignmentInput& window)
{
    // train on current consensus
    uint32_t num_segments = window.anchored_columns.size();
    for(uint32_t segment_id = 0; segment_id < num_segments - 2; ++segment_id) {
        train_segment(window, segment_id);
    }

    // Update model parameters
    for(uint32_t ri = 0; ri < window.reads.size(); ++ri) {
        window.reads[ri]->parameters[0].train();
        window.reads[ri]->parameters[1].train();
    }
}

std::string call_consensus_for_window(const Fast5Map& name_map, const std::string& contig, int start_base, int end_base)
{
    const int minor_segment_stride = 50;
    HMMRealignmentInput window = build_input_for_region(opt::bam_file,
                                                        opt::genome_file,
                                                        name_map,
                                                        contig,
                                                        start_base,
                                                        end_base,
                                                        minor_segment_stride);
    uint32_t num_segments = window.anchored_columns.size();

    // If there are not reads or not enough segments do not try to call a consensus sequence
    if(window.reads.empty() || num_segments < 3) {
        // No data for this window, just return the original sequence as the consensus
        assert(!window.original_sequence.empty());
        return window.original_sequence;
    }
    
    if(opt::verbose > 0) {
        fprintf(stderr, "correcting window %s:%d-%d with %zu reads\n", contig.c_str(), start_base, end_base, window.reads.size());
    }

    //
    // Train the HMM
    //
    WARN_ONCE("Debug: using default transition parameters"); 
    //train(window);

    // assume models for all the reads have the same k
    const uint32_t k = window.reads[0]->pore_model[T_IDX].k;

    //
    // Compute the new consensus sequence
    //
    std::string reference = "";
    std::string consensus = "";

    uint32_t start_segment_id = 0;

    // Copy the base segments before they are updated
    // by the consensus algorithm
    std::vector<std::string> ref_segments;
    for(uint32_t segment_id = 0; segment_id < num_segments; ++segment_id) {
        ref_segments.push_back(window.anchored_columns[segment_id].base_sequence);
    }

    // Initialize progress status
    std::stringstream message;
    message << "[consensus] " << contig << ":" << start_base << "-" << end_base;
    Progress progress(message.str());

    for(uint32_t segment_id = start_segment_id; segment_id < num_segments - 2; ++segment_id) {

        // update progress
        if(opt::show_progress) {
            progress.print((float)segment_id / (num_segments - 2));
        }

        // run the consensus algorithm for this segment
        run_splice_segment(window, segment_id, k);

        // run_splice_segment updates the base_sequence of the current anchor, grab it and append
        std::string base = window.anchored_columns[segment_id].base_sequence;

        // append the new sequences in, respecting the K overlap
        reference = join_sequences_at_kmer(reference, ref_segments[segment_id], k);
        consensus = join_sequences_at_kmer(consensus, base, k);

        if(opt::verbose > 0) {
            fprintf(stderr, "UNCORRECT[%d]: %s\n", segment_id, reference.c_str());
            fprintf(stderr, "CONSENSUS[%d]: %s\n", segment_id, consensus.c_str());
        }
    }

    // Append segment that ends at the last anchor
    reference = join_sequences_at_kmer(reference, ref_segments[num_segments - 2], k);
    const std::string& last_segment = 
        window.anchored_columns[num_segments - 2].base_sequence;
    consensus = join_sequences_at_kmer(consensus, last_segment, k);

    if(opt::show_progress) {
        progress.end();
    }

    return consensus;
}

void parse_consensus_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'w': arg >> opt::window; break;
            case 'o': arg >> opt::output_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'v': opt::verbose++; break;
            case OPT_MODELS_FOFN: arg >> opt::models_fofn; break;
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

    if(opt::window.empty()) {
        std::cerr << SUBPROGRAM ": the -w (window) parameter must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << CONSENSUS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int consensus_main(int argc, char** argv)
{
    parse_consensus_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    
    // Parse the window string
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

    FILE* out_fp = NULL;

    if(!opt::output_file.empty()) {
        out_fp = fopen(opt::output_file.c_str(), "w");
    } else {
        out_fp = stdout;
    }

    for(int window_id = start_window_id; window_id < end_window_id; ++window_id) {
        int start_base = window_id * WINDOW_LENGTH;
        int end_base = start_base + WINDOW_LENGTH + WINDOW_OVERLAP;
        
        std::string window_consensus = call_consensus_for_window(name_map, contig, start_base, end_base);
        fprintf(out_fp, ">%s:%d\n%s\n", contig.c_str(), window_id, window_consensus.c_str());
    }

    if(out_fp != stdout) {
        fclose(out_fp);
    }
    return 0;
}
