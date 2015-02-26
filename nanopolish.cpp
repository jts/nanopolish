//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish.cpp -- entry point to consensus functions
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
#include "nanopolish_poremodel.h"
#include "nanopolish_interface.h"
#include "nanopolish_khmm_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_khmm.h"
#include "nanopolish_profile_hmm.h"
#include "profiler.h"

// Macros
#define max3(x,y,z) std::max(std::max(x,y), z)

// Constants

// strands
const uint8_t T_IDX = 0;
const uint8_t C_IDX = 1;
const uint8_t NUM_STRANDS = 2;

const static double LOG_KMER_INSERTION = log(0.1);
const static double P_RANDOM_SKIP = 0.05;
const static double EVENT_DETECTION_THRESHOLD = 1.0f;

// Flags to turn on/off debugging information

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1
//#define DEBUG_TRANSITION 1
#define DEBUG_PATH_SELECTION 1
//#define DEBUG_SINGLE_SEGMENT 1
//#define DEBUG_SHOW_TOP_TWO 1
//#define DEBUG_SEGMENT_ID 97

struct HMMReadAnchor
{
    int32_t event_idx;
    bool rc; // with respect to consensus
};

struct HMMAnchoredColumn
{
    std::vector<HMMReadAnchor> anchors;
    std::string base_sequence;
    std::vector<std::string> alt_sequences;
};

// A global vector used to store data we've received from the python code
struct HmmConsData
{
    int num_threads;

    //
    std::vector<SquiggleRead> reads;
    std::vector<HMMAnchoredColumn> anchored_columns;
    
    //
    std::string consensus_result;
};

HmmConsData g_data;
bool g_initialized = false;

extern "C"
void initialize(int num_threads)
{
    g_initialized = true;
    g_data.num_threads = num_threads;
}

extern "C"
void clear_data()
{
    g_data.reads.clear();
    g_data.anchored_columns.clear();
    g_data.consensus_result.clear();
}

extern "C"
const char* get_consensus_result()
{
    return g_data.consensus_result.c_str();
}

extern "C"
void add_read(CSquiggleReadInterface params)
{
    g_data.reads.push_back(SquiggleRead());

    SquiggleRead& sr = g_data.reads.back();
    sr.read_id = g_data.reads.size() - 1;

    for(uint32_t i = 0; i < NUM_STRANDS; ++i) {
        // Initialize pore model   
        sr.pore_model[i].scale = params.pore_model[i].scale;
        sr.pore_model[i].shift = params.pore_model[i].shift;
        sr.pore_model[i].drift = params.pore_model[i].drift;
        sr.pore_model[i].var = params.pore_model[i].var;
        sr.pore_model[i].scale_sd = params.pore_model[i].scale_sd;
        sr.pore_model[i].var_sd = params.pore_model[i].var_sd;
        
        assert(params.pore_model[i].n_states == 1024);
        for(uint32_t j = 0; j < params.pore_model[i].n_states; ++j) {
            
            sr.pore_model[i].state[j].level_mean = params.pore_model[i].level_mean[j];
            sr.pore_model[i].state[j].level_stdv = params.pore_model[i].level_stdv[j];
            
            sr.pore_model[i].state[j].sd_mean = params.pore_model[i].sd_mean[j];
            sr.pore_model[i].state[j].sd_stdv = params.pore_model[i].sd_stdv[j];
         }
    
        // Initialize events
        sr.events[i].n_events = params.events[i].n_events;
        sr.events[i].level = params.events[i].level;
        sr.events[i].stdv = params.events[i].stdv;
        sr.events[i].time = params.events[i].time;

        /*
        printf("Model[%zu] scale: %lf shift: %lf %lf %lf\n", i, sr.pore_model[i].scale, 
                                                                 sr.pore_model[i].shift,
                                                                 sr.pore_model[i].state[0].level_mean, 
                                                                 sr.pore_model[i].state[0].level_stdv);
    
        printf("First 100 events of %d\n", sr.events[i].n_events);
        for(int j = 0; j < 100; ++j)
            printf("%d: %lf\n", j, sr.events[i].level[j]);
        */
    }

    // Initialize hmm parameters for both strands of the read
    khmm_parameters_initialize(sr.parameters[0]);
    khmm_parameters_initialize(sr.parameters[1]);
}

// This is called by python to tell us we want to start a new anchored column
extern "C"
void start_anchored_column()
{
    HMMAnchoredColumn ac;
    g_data.anchored_columns.push_back(ac);
}

extern "C"
void add_read_anchor(CReadAnchorInterface in_ra)
{
    assert(!g_data.anchored_columns.empty());

    HMMReadAnchor ra = { in_ra.event_idx, in_ra.rc };
    g_data.anchored_columns.back().anchors.push_back(ra);
}

extern "C"
void add_base_sequence(char* str)
{
    assert(!g_data.anchored_columns.empty());
    g_data.anchored_columns.back().base_sequence = str;
}

extern "C"
void add_alt_sequence(char* str)
{
    assert(!g_data.anchored_columns.empty());
    g_data.anchored_columns.back().alt_sequences.push_back(str);
}

// This is called by python to tell us we want to start a new anchored column
extern "C"
void end_anchored_column()
{
    // Validate that we received two read anchors per read
    assert(g_data.anchored_columns.back().anchors.size() == g_data.reads.size() * 2);
}

std::vector<HMMConsReadState> get_read_states_for_columns(const HMMAnchoredColumn& start_column,  
                                                          const HMMAnchoredColumn& end_column)
{
    assert(start_column.anchors.size() == end_column.anchors.size());

    std::vector<HMMConsReadState> read_states;
    for(uint32_t rsi = 0; rsi < start_column.anchors.size(); ++rsi) {

        HMMReadAnchor start_ra = start_column.anchors[rsi];
        HMMReadAnchor end_ra = end_column.anchors[rsi];

        // This read strand does not have events at both anchors
        if(start_ra.event_idx == -1 || end_ra.event_idx == -1)
            continue;

        HMMConsReadState crs;

        uint32_t read_idx = rsi / 2;
        assert(read_idx < g_data.reads.size());
        crs.anchor_index = rsi;
        crs.read = &g_data.reads[read_idx];
        crs.strand = rsi % 2;
        crs.event_start_idx = start_ra.event_idx;
        crs.event_stop_idx = end_ra.event_idx;
        if(crs.event_start_idx < crs.event_stop_idx)
            crs.stride = 1;
        else
            crs.stride = -1;
        assert(start_ra.rc == end_ra.rc);
        crs.rc = start_ra.rc;

        read_states.push_back(crs);
    }
    return read_states;
}

// Handy wrappers for scoring/debugging functions
// The consensus algorithms call into these so we can switch
// scoring functinos without writing a bunch of code
double score_sequence(const std::string& sequence, const HMMConsReadState& state)
{
    //return score_skip_merge(sequence, state);
    //return score_khmm_model_postmerge(sequence, state);
    //return khmm_score(sequence, state, AP_GLOBAL);
    return profile_hmm_score(sequence, state);
    //return score_emission_dp(sequence, state);
}


std::vector<AlignmentState> hmm_align(const std::string& sequence, const HMMConsReadState& state)
{
    return profile_hmm_align(sequence, state);
//    return khmm_posterior_decode(sequence, state);
}

void debug_sequence(const std::string& name, uint32_t seq_id, uint32_t read_id, const std::string& sequence, const HMMConsReadState& state)
{
    std::vector<AlignmentState> alignment = hmm_align(sequence, state);
    print_alignment(name, seq_id, read_id, sequence, state, alignment);
}

void update_training_with_segment(const std::string& sequence, const HMMConsReadState& state)
{
    profile_hmm_update_training(sequence, state);
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
void score_paths(PathConsVector& paths, const std::vector<HMMConsReadState>& read_states)
{
    PROFILE_FUNC("score_paths")
    double MIN_FIT = INFINITY;
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
    for(uint32_t ri = 0; ri < read_states.size(); ++ri) {
        printf("Scoring %d\n", ri);

        const HMMConsReadState& read_state = read_states[ri];
        const KHMMParameters& parameters = read_state.read->parameters[read_state.strand];
 
        if( fabs(parameters.fit_quality) > MIN_FIT)
            continue;

        std::vector<IndexedPathScore> result(paths.size());

        // Score all paths
        omp_set_num_threads(g_data.num_threads);
        #pragma omp parallel for
        for(size_t pi = 0; pi < paths.size(); ++pi) {
            double curr = score_sequence(paths[pi].path, read_states[ri]);
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
    for(size_t pi = 0; pi < paths.size(); ++pi) {

        // Calculate the length of the matching prefix with the initial sequence
        const std::string& s = paths[pi].path;

        char initial = s == first ? 'I' : ' ';

        printf("%zu\t%s\t%.1lf\t%zu %c %s", pi, paths[pi].path.c_str(), paths[pi].score, paths[pi].sum_rank, initial, paths[pi].mutdesc.c_str());
        // If this is the truth path or the best path, show the scores for all reads
        if(pi <= 1 || initial == 'I') {
            for(uint32_t ri = 0; ri < read_states.size(); ++ri) {
                const HMMConsReadState& read_state = read_states[ri];
                const KHMMParameters& parameters = read_state.read->parameters[read_state.strand];
                if( fabs(parameters.fit_quality) > MIN_FIT)
                    continue;

                double curr = score_sequence(paths[pi].path, read_states[ri]);
                printf("%.1lf,%.2lf ", parameters.fit_quality, curr);
            }
        }
        printf("\n");
    }
#endif

}

void extend_paths(PathConsVector& paths, int maxk = 2)
{
    // Insert all possible extensions into the path sequence
    // for k in 1 to maxk
    PathConsVector new_paths;

    for(int k = 1; k <= maxk; ++k) {

        for(int pi = 0; pi < paths.size(); ++pi) {
    
            std::string first(k, 'A');
            std::string extension = first;

            do {
                std::string current = paths[pi].path;
                std::string ns = current.insert(current.size() - 5, extension);
                PathCons ps(ns);
                new_paths.push_back(ps);
                lexicographic_next(extension);
            } while(extension != first);
        }
    }

    paths.swap(new_paths);
}

PathConsVector generate_mutations(const std::string& sequence)
{
    PathConsVector mutations;

    // Add the unmutated sequence
    {
        PathCons pc(sequence);
        mutations.push_back(pc);
    }

    // Mutate every base except for in the first/last k-mer
    for(size_t si = K; si < sequence.size() - K; ++si) {
        
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

void run_mutation(std::string& base, const std::vector<HMMConsReadState>& read_states, std::string& second_best)
{
    PROFILE_FUNC("run_mutation")
    int iteration = 0;
    while(iteration++ < 10) {

        // Generate possible sequences
        PathConsVector paths = generate_mutations(base);

        score_paths(paths, read_states);

        second_best = paths[1].path;
        // check if no improvement was made
        if(paths[0].path == base)
            break;
        base = paths[0].path;
    }
}

void generate_alt_paths(PathConsVector& paths, const std::string& base, const std::vector<std::string>& alts)
{
    // Generate alternatives
    for(uint32_t ai = 0; ai < alts.size(); ++ai) {
        const std::string& alt = alts[ai];
        kLCSResult result = kLCS(base, alt, K);

#ifdef DEBUG_ALT_GENERATION
        printf("Match to alt %s\n", alt.c_str());
        for(size_t mi = 0; mi < result.size(); ++mi) {
            std::string extend = "";
            if(mi < result.size() - 1 && result[mi].j + 1 != result[mi + 1].j) {
                extend = alt.substr(result[mi].j, result[mi + 1].j - result[mi].j + K);
            }
            printf("\t%zu %zu %s %s\n", result[mi].i, result[mi].j, base.substr(result[mi].i, K).c_str(), extend.c_str());
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

//
// Outlier filtering
//
void filter_outlier_read_states(std::vector<HMMConsReadState>& read_states, const std::string& sequence)
{
    std::vector<HMMConsReadState> out_rs;
    for(uint32_t ri = 0; ri < read_states.size(); ++ri) {
        const HMMConsReadState& rs = read_states[ri];

        double curr = score_sequence(sequence, rs);
        double n_events = abs(rs.event_start_idx - rs.event_stop_idx) + 1.0f;
        double lp_per_event = curr / n_events;
        printf("OUTLIER_FILTER %d %.2lf %.2lf %.2lf\n", ri, curr, n_events, lp_per_event);
        if(fabs(lp_per_event) < 3.5f) {
            out_rs.push_back(rs);
        }
    }
    read_states.swap(out_rs);
}

std::string join_sequences_at_kmer(const std::string& a, const std::string& b)
{
    // These sequences must have a k-mer match at the start/end
    std::string a_last_kmer = a.substr(a.size() - K);
    std::string b_last_kmer = b.substr(0, K);
    assert(a_last_kmer == b_last_kmer);
    return a + b.substr(K);
}

void run_splice_segment(uint32_t segment_id)
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

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
    assert(segment_id + 2 < g_data.anchored_columns.size());
    HMMAnchoredColumn& start_column = g_data.anchored_columns[segment_id];
    HMMAnchoredColumn& middle_column = g_data.anchored_columns[segment_id + 1];
    HMMAnchoredColumn& end_column = g_data.anchored_columns[segment_id + 2];

    std::string s_m_base = start_column.base_sequence;
    std::string m_e_base = middle_column.base_sequence;

    // The current consensus sequence
    std::string original = join_sequences_at_kmer(s_m_base, m_e_base);
    std::string base = original;
    
    // The collection of alternative sequences
    std::vector<std::string> alts;

    for(uint32_t ai = 0; ai < start_column.alt_sequences.size(); ++ai) {
        alts.push_back(start_column.alt_sequences[ai]);
    }

    // Set up the HMMReadStates, which are used to calculate
    // the probability of the data given a possible consensus sequence
    std::vector<HMMConsReadState> read_states = get_read_states_for_columns(start_column, end_column);

    //
    filter_outlier_read_states(read_states, base);

    // Only attempt correction if there are any reads here
    if(!read_states.empty()) {
        uint32_t num_rounds = 6;
        uint32_t round = 0;
        while(round++ < num_rounds) {
            
            PathConsVector paths;
            PathCons base_path(base);
            paths.push_back(base_path);
            
            generate_alt_paths(paths, base, alts);
            score_paths(paths, read_states);

            if(paths[0].path == base)
                break;
            base = paths[0].path;
        }

        std::string second_best;
        run_mutation(base, read_states, second_best);

#if DEBUG_SHOW_TOP_TWO
        assert(!second_best.empty());
        for(uint32_t ri = 0; ri < read_states.size(); ++ri) {
            debug_sequence("best", segment_id, ri, base, read_states[ri]);
            debug_sequence("second", segment_id, ri, second_best, read_states[ri]);
        }
#endif
    }

    printf("ORIGINAL[%d] %s\n", segment_id, original.c_str());
    printf("RESULT[%d]   %s\n", segment_id, base.c_str());

    // Update the sequences for the start and middle segments
    // by cutting the new consensus in the middle
    // We maintain the k-mer match invariant by requiring the
    // sequences to overlap by 5bp
    assert(base.length() >= K);
    uint32_t midpoint_kmer = (base.length() - K + 1) / 2;

    std::string s_m_fixed = base.substr(0, midpoint_kmer + K);
    std::string m_e_fixed = base.substr(midpoint_kmer);

    assert(s_m_fixed.substr(s_m_fixed.size() - K) == m_e_fixed.substr(0, K));

    start_column.base_sequence = s_m_fixed;
    middle_column.base_sequence = m_e_fixed;

    // Update the event indices in the first column to match 
    for(uint32_t ri = 0; ri < read_states.size(); ++ri) {

        // Realign to the consensus sequence
        std::vector<AlignmentState> decodes = hmm_align(base, read_states[ri]);

        // Get the closest event aligned to the target kmer
        int32_t min_k_dist = base.length();
        uint32_t event_idx = 0;
        for(uint32_t di = 0; di < decodes.size(); ++di) {
            int32_t dist = abs(decodes[di].kmer_idx - midpoint_kmer);
            if(dist <= min_k_dist) {
                min_k_dist = dist;
                event_idx = decodes[di].event_idx;
            }
        }

        middle_column.anchors[read_states[ri].anchor_index].event_idx = event_idx;
    }
}

extern "C"
void run_splice()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }
 
    std::string uncorrected = "";
    std::string consensus = "";

    uint32_t start_segment_id = 0;
#ifdef DEBUG_SINGLE_SEGMENT
    start_segment_id = DEBUG_SEGMENT_ID;
#endif

    uint32_t num_segments = g_data.anchored_columns.size();
    for(uint32_t segment_id = start_segment_id; segment_id < num_segments - 2; ++segment_id) {

        // Track the original sequence for reference
        if(uncorrected.empty()) {
            uncorrected = g_data.anchored_columns[segment_id].base_sequence;
        } else {
            uncorrected.append(g_data.anchored_columns[segment_id].base_sequence.substr(K));
        }

        // run the consensus algorithm for this segment
        run_splice_segment(segment_id);

        // run_splice_segment updates the base_sequence of the current anchor, grab it and append
        std::string base = g_data.anchored_columns[segment_id].base_sequence;

        if(consensus.empty()) {
            consensus = base;
        } else {
            // The first 5 bases of the incoming sequence must match
            // the last 5 bases of the growing consensus
            // run_splice_segment must ensure this
            assert(consensus.substr(consensus.size() - K) == base.substr(0, K));
            consensus.append(base.substr(K));
        }

        printf("UNCORRECT[%d]: %s\n", segment_id, uncorrected.c_str());
        printf("CONSENSUS[%d]: %s\n", segment_id, consensus.c_str());
#ifdef DEBUG_SINGLE_SEGMENT
        break;
#endif
    }

    g_data.consensus_result = consensus;
}

// update the training data on the current segment
extern "C"
void train_segment(uint32_t segment_id)
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    // Get the segments
    assert(segment_id + 2 < g_data.anchored_columns.size());
    HMMAnchoredColumn& start_column = g_data.anchored_columns[segment_id];
    HMMAnchoredColumn& middle_column = g_data.anchored_columns[segment_id + 1];
    HMMAnchoredColumn& end_column = g_data.anchored_columns[segment_id + 2];

    std::string s_m_base = start_column.base_sequence;
    std::string m_e_base = middle_column.base_sequence;

    std::string segment_sequence = join_sequences_at_kmer(s_m_base, m_e_base);

    // Set up the HMMReadStates, which are used to calculate
    // the probability of the data given a possible consensus sequence
    std::vector<HMMConsReadState> read_states = get_read_states_for_columns(start_column, end_column);
     
    for(uint32_t ri = 0; ri < read_states.size(); ++ri) {

        std::vector<AlignmentState> decodes = hmm_align(segment_sequence, read_states[ri]);
        update_training_with_segment(segment_sequence, read_states[ri]);
    }
}

extern "C"
void train()
{
    // train on current consensus
    uint32_t num_segments = g_data.anchored_columns.size();
    for(uint32_t segment_id = 0; segment_id < num_segments - 2; ++segment_id) {
        printf("Training segment %d\n", segment_id);
        train_segment(segment_id);
    }

    // Update model parameters
    for(uint32_t ri = 0; ri < g_data.reads.size(); ++ri) {
        khmm_parameters_train(g_data.reads[ri].parameters[0]);
        khmm_parameters_train(g_data.reads[ri].parameters[1]);
    }
}

int main(int argc, char** argv)
{

}
