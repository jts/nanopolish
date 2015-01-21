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
#include "hmmcons_poremodel.h"
#include "hmmcons_interface.h"
#include "Profiler.h"

// Macros
#define max3(x,y,z) std::max(std::max(x,y), z)

// Constants

// strands
const uint8_t T_IDX = 0;
const uint8_t C_IDX = 1;
const uint8_t NUM_STRANDS = 2;

// 
const uint8_t K = 5;

const static double LOG_KMER_INSERTION = log(0.1);
const static double SELF_KMER_TRANSITION = 0.2;
const static double P_RANDOM_SKIP = 0.015; // 0.015 is estimated from the event duration distribution for sample data
const static double EVENT_DETECTION_THRESHOLD = 1.0f;
const static uint32_t KHMM_MAX_JUMP = 5;

static const uint8_t base_rank[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

enum AlignmentPolicy
{
    AP_GLOBAL,
    AP_SEMI_KMER
};

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1
//#define DEBUG_TRANSITION 1
struct CEventSequence
{
    uint32_t n_events;
    const double* level;
    const double* stdv;
    const double* time;
};

struct CSquiggleRead
{
    // unique identifier of the read
    uint32_t read_id;

    // one model for each strand
    CPoreModel pore_model[2];

    // one event sequence for each strand as well
    CEventSequence events[2];
};

struct HMMConsReadState
{
    CSquiggleRead* read;
    uint32_t event_start_idx;
    uint32_t event_stop_idx;
    uint32_t kmer_idx;
    uint8_t strand;
    int8_t stride;
    uint8_t rc;
    std::string alignment;
};

struct ExtensionResult
{
    double b[4];
    std::string best_path;
    double best_path_score;
};

//
//
//
struct HMMParameters
{
    // transition matrix
    const static uint32_t n_states = 4;

    // The transition matrix is described using pseudocounts, within
    // the initialize function it is normalized and log_scaled
    double t[n_states][n_states];
};

void initialize_hmm(HMMParameters& params)
{
    // These values are trained from ONT's 2D alignments

    // transitions from the match state
    params.t[0][0] = 8980.f; // M
    params.t[0][1] = 2147.f; // E
    params.t[0][2] = 1973.f; // K
    params.t[0][3] =  200.f; // terminal
    
    // transitions from the event insertion state
    params.t[1][0] = 1746.f; // M
    params.t[1][1] = 1352.f; // E
    params.t[1][2] =  401.f; // K
    params.t[1][3] =   35.f; // terminal

    // transitions from the k-mer insertion state
    params.t[2][0] = 2373.f; // M
    params.t[2][1] =    0.f; // E
    params.t[2][2] =  519.f; // K
    params.t[2][3] =   30.f; // terminal

    // transitions from the terminal state
    params.t[3][0] = 0.f;  // M
    params.t[3][1] = 0.f;  // E
    params.t[3][2] = 0.f;  // K
    params.t[3][3] = 1.f;  // terminal

    // Row normalize and log scale
    for(uint32_t i = 0; i < params.n_states; ++i) {
        double sum = 0.0f;
        for(uint32_t j = 0; j < params.n_states; ++j)
            sum += params.t[i][j];

        for(uint32_t j = 0; j < params.n_states; ++j)
            params.t[i][j] = log(params.t[i][j] / sum);
    }
}

// A global vector used to store data we've received from the python code
struct HmmConsData
{
    //
    std::vector<CSquiggleRead> reads;

    //
    std::vector<HMMConsReadState> read_states;
    std::vector<std::string> candidate_consensus;
    std::string consensus_result;
    
    //
    HMMParameters hmm_params;
};
HmmConsData g_data;
bool g_initialized = false;

extern "C"
void initialize()
{
    initialize_hmm(g_data.hmm_params);
    g_initialized = true;
}

extern "C"
void clear_state()
{
    g_data.read_states.clear();
    g_data.candidate_consensus.clear();
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
    g_data.reads.push_back(CSquiggleRead());
    CSquiggleRead& sr = g_data.reads.back();
 
    for(uint32_t i = 0; i < NUM_STRANDS; ++i) {
        // Initialize pore model   
        sr.pore_model[i].scale = params.pore_model[i].scale;
        sr.pore_model[i].shift = params.pore_model[i].shift;
        sr.pore_model[i].drift = params.pore_model[i].drift;
        sr.pore_model[i].var = params.pore_model[i].var;
        
        assert(params.pore_model[i].n_states == 1024);
        for(uint32_t j = 0; j < params.pore_model[i].n_states; ++j) {
            sr.pore_model[i].state[j].mean = params.pore_model[i].mean[j];
            sr.pore_model[i].state[j].sd = params.pore_model[i].sd[j];
         }
    
        // Initialize events
        sr.events[i].n_events = params.events[i].n_events;
        sr.events[i].level = params.events[i].level;
        sr.events[i].stdv = params.events[i].stdv;
        sr.events[i].time = params.events[i].time;

        /*
        printf("Model[%zu] scale: %lf shift: %lf %lf %lf\n", i, sr.pore_model[i].scale, 
                                                                 sr.pore_model[i].shift,
                                                                 sr.pore_model[i].state[0].mean, 
                                                                 sr.pore_model[i].state[0].sd);
    
        printf("First 100 events of %d\n", sr.events[i].n_events);
        for(int j = 0; j < 100; ++j)
            printf("%d: %lf\n", j, sr.events[i].level[j]);
        */
    }
}

extern "C"
void add_read_state(CReadStateInterface params)
{
    g_data.read_states.push_back(HMMConsReadState());
    HMMConsReadState& rs = g_data.read_states.back();
    rs.read = &g_data.reads[params.read_idx];
    rs.event_start_idx = params.event_start_idx;
    rs.event_stop_idx = params.event_stop_idx;
    rs.kmer_idx = 0;
    rs.strand = params.strand;
    rs.stride = params.stride;
    rs.rc = params.rc;
}

extern "C"
void add_candidate_consensus(char* str)
{
    g_data.candidate_consensus.push_back(str);
}


//
// HMM matrix
//
struct HMMCell
{
    double M;
    double E;
    double K;
};

struct HMMMatrix
{
    HMMCell* cells;
    uint32_t n_rows;
    uint32_t n_cols;
};

void allocate_matrix(HMMMatrix& matrix, uint32_t n_rows, uint32_t n_cols)
{
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    uint32_t N = matrix.n_rows * matrix.n_cols;
    matrix.cells = (HMMCell*)malloc(N * sizeof(HMMCell));
    memset(matrix.cells, 0, N * sizeof(HMMCell));
}

void free_matrix(HMMMatrix& matrix)
{
    free(matrix.cells);
    matrix.cells = NULL;
}

inline uint32_t cell(const HMMMatrix& matrix, uint32_t row, uint32_t col)
{
    return row * matrix.n_cols + col;
}

void print_matrix(const HMMMatrix& matrix)
{
    for(uint32_t i = 0; i < matrix.n_rows; ++i) {
        for(uint32_t j = 0; j < matrix.n_cols; ++j) {
            uint32_t c = cell(matrix, i, j);
            printf("%.1lf,%.1lf,%.1f\t", matrix.cells[c].M, matrix.cells[c].E, matrix.cells[c].K);
        }
        printf("\n");
    }
}

//
// Double Matrix
//
struct DoubleMatrix
{
    double* cells;
    uint32_t n_rows;
    uint32_t n_cols;
};

void allocate_matrix(DoubleMatrix& matrix, uint32_t n_rows, uint32_t n_cols)
{
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    
    uint32_t N = matrix.n_rows * matrix.n_cols;
    matrix.cells = (double*)malloc(N * sizeof(double));
    memset(matrix.cells, 0, N * sizeof(double));
}

void free_matrix(DoubleMatrix& matrix)
{
    assert(matrix.cells != NULL);
    free(matrix.cells);
    matrix.cells = NULL;
}

inline uint32_t cell(const DoubleMatrix& matrix, uint32_t row, uint32_t col)
{
    return row * matrix.n_cols + col;
}

inline double set(const DoubleMatrix& matrix, uint32_t row, uint32_t col, double v)
{
    uint32_t c = cell(matrix, row, col);
    matrix.cells[c] = v;
}

inline double get(const DoubleMatrix& matrix, uint32_t row, uint32_t col)
{
    uint32_t c = cell(matrix, row, col);
    return matrix.cells[c];
}

void print_matrix(const DoubleMatrix& matrix, bool do_exp = false)
{
    for(uint32_t i = 0; i < matrix.n_rows; ++i) {
        for(uint32_t j = 0; j < matrix.n_cols; ++j) {
            uint32_t c = cell(matrix, i, j);
            double v = matrix.cells[c];
            if(do_exp)
                v = exp(v);
            printf("%.3lf\t", v);
        }
        printf("\n");
    }
}

//
// Kmer Ranks
//
inline uint32_t kmer_rank(const char* str, uint32_t K)
{
    uint32_t rank = 0;
    for(uint32_t i = 0; i < K; ++i)
        rank |= base_rank[str[i]] << 2 * (K - i - 1);
    return rank;
}

inline uint32_t rc_kmer_rank(const char* str, uint32_t K)
{
    uint32_t rank = 0;
    for(int32_t i = K - 1; i >= 0; --i)
        rank |= ((3 - base_rank[str[i]]) << 2 * i);
    return rank;
}


// Increment the input string to be the next sequence in lexicographic order
void lexicographic_next(std::string& str)
{
    int carry = 1;
    int i = str.size() - 1;
    do {
        uint32_t r = base_rank[str[i]] + carry;
        str[i] = "ACGT"[r % 4];
        carry = r / 4;
        i -= 1;
    } while(carry > 0 && i >= 0);
}

// From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
// TODO: replace with a lookup table that can be interpolated
inline double normal_pdf(double x, double m, double s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

inline double log_normal_pdf(double x, double m, double s)
{
    static const double log_inv_sqrt_2pi = log(0.3989422804014327);
    double a = (x - m) / s;
    return log_inv_sqrt_2pi - log(s) + (-0.5f * a * a);
}

// The probability that a standard normal RV is <= x
// from http://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
inline double log_standard_normal_cdf(double x)
{
    return 0.5 * erfc(-x * M_SQRT1_2);
}

// The probability that a normal RV is <= x
inline double log_normal_cdf(double x, double m, double s)
{
    double a = (x - m) / s;
    return log(0.5 * (1 + erf(a * M_SQRT1_2)));
}

inline double log_probability_match(const CSquiggleRead& read,
                                    uint32_t kmer_rank,
                                    uint32_t event_idx, 
                                    uint8_t strand)
{
    const CPoreModel& pm = read.pore_model[strand];

    // Extract event
    double level = read.events[strand].level[event_idx];
    
    // correct level by drift
    double start = read.events[strand].time[0];
    double time = read.events[strand].time[event_idx] - start;
    level -= (time * pm.drift);

    double m = pm.state[kmer_rank].mean * pm.scale + pm.shift;
    double s = pm.state[kmer_rank].sd * pm.var;
    double lp = log_normal_pdf(level, m, s);

#if DEBUG_HMM_EMISSION
    printf("Event[%d] Kmer: %d -- L:%.1lf m: %.1lf s: %.1lf p: %.3lf p_old: %.3lf\n", event_idx, kmer_rank, level, m, s, exp(lp), normal_pdf(level, m, s));
#endif

    return lp;
}

inline double log_probability_event_insert(const CSquiggleRead& read,
                                           uint32_t kmer_rank,
                                           uint32_t event_idx, 
                                           uint8_t strand)
{
    return log_probability_match(read, kmer_rank, event_idx, strand);
}

inline double log_probability_kmer_insert(const CSquiggleRead& read,
                                          uint32_t kmer_rank,
                                          uint32_t event_idx, 
                                          uint8_t strand)

{
    return log_probability_match(read, kmer_rank, event_idx, strand);
}

void fill_khmm_transitions(DoubleMatrix& matrix, const std::string& consensus, const HMMConsReadState& state)
{
    PROFILE_FUNC("fill_khmm_transitions")

    const CPoreModel& pm = state.read->pore_model[state.strand];
    uint32_t n_kmers = consensus.size() - K + 1;
    uint32_t n_states = n_kmers + 2;
    uint32_t terminal_state = n_states - 1;

    assert(matrix.n_rows == n_states && matrix.n_cols == n_states);

    // Initialize the transition matrix to -INFINITY for all states
    for(size_t si = 0; si < n_states; ++si)
        for(size_t sj = 0; sj < n_states; ++sj)
            set(matrix, si, sj, -INFINITY);

    // Start state transitions -- only allowed to go to k_0
    set(matrix, 0, 1, 0.0f);
    
    // TODO: calculate in log space
    for(size_t si = 1; si < n_states - 1; si++) {
        size_t ki = si - 1;
        double sum = 0.0f;

        uint32_t last_valid_state = si + KHMM_MAX_JUMP;
        if(last_valid_state >= terminal_state)
            last_valid_state = terminal_state - 1;

        for(size_t sj = si; sj <= last_valid_state; ++sj) {
            
            size_t kj = sj - 1;
                        
            // transition probability
            double p_i_j = 0.0f;

            if(ki == kj) {
                p_i_j = SELF_KMER_TRANSITION;
            } else {
        
                uint32_t rank_i = !state.rc ? 
                    kmer_rank(consensus.c_str() + ki, K) : 
                    rc_kmer_rank(consensus.c_str() + ki, K);
                
                uint32_t rank_j = !state.rc ? 
                    kmer_rank(consensus.c_str() + kj, K) : 
                    rc_kmer_rank(consensus.c_str() + kj, K);
        
                double level_i = (pm.state[rank_i].mean + pm.shift) * pm.scale;
                double sd_i = pm.state[rank_i].sd * pm.scale;
                
                double level_j = (pm.state[rank_j].mean + pm.shift) * pm.scale;
                double sd_j = pm.state[rank_j].sd * pm.scale;

                double diff_mean = level_i - level_j;
                double diff_sd = sqrt(pow(sd_i, 2.0) + pow(sd_j, 2.0));

                double p_within_threshold = exp(log_normal_cdf(EVENT_DETECTION_THRESHOLD, diff_mean, diff_sd)) - 
                                            exp(log_normal_cdf(-EVENT_DETECTION_THRESHOLD, diff_mean, diff_sd));
                double p_skip = P_RANDOM_SKIP + (1 - P_RANDOM_SKIP) * p_within_threshold;

                p_i_j = (1 - sum) * (1 - p_skip);
#ifdef DEBUG_TRANSITION
                printf("\t\t%zu -> %zu %.2lf %.2lf %.2lf %.2lf p_within_threshold: %.4lf p_skip: %.4lf p: %.2lf\n", ki, kj, level_i, level_j, diff_mean, diff_sd, p_within_threshold, p_skip, p_i_j);
#endif
            }

            sum += p_i_j;
            set(matrix, si, sj, log(p_i_j));
        }
    }

    // Transition to end state -- only the last k-mer can go to the end state
    // TODO: allow the last k-mer to be skipped??
    set(matrix, n_states - 2, n_states - 1, 0.0f);
}

void initialize_forward_khmm(DoubleMatrix& fm)
{
    // initialize forward calculation
    for(uint32_t si = 0; si < fm.n_cols; si++)
        set(fm, 0, si, -INFINITY);
    for(uint32_t ri = 0; ri < fm.n_rows; ri++)
        set(fm, ri, 0, -INFINITY);

    set(fm, 0, 0, 0.0f); // probability 1 in the start state for the null row
}

// Terminate the forward algorithm by calculating
// the probability of transitioning to the end state
// for all columns and a given row
double forward_khmm_terminate(const DoubleMatrix& fm,
                              const DoubleMatrix& tm,
                              uint32_t row)
{
    double sum = -INFINITY;
    uint32_t tcol = fm.n_cols - 1;
    for(uint32_t sk = 0; sk < fm.n_cols - 1; sk++) {
        // transition probability from state k to state l
        double t_kl = get(tm, sk, tcol);
        double fm_k = get(fm, row, sk);
        sum = log(exp(sum) + exp(t_kl + fm_k));
    }
    return sum;
}

double fill_forward_khmm(DoubleMatrix& fm, // forward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMConsReadState& state,
                         uint32_t e_start)
{
    PROFILE_FUNC("fill_forward_khmm")

    // Fill in matrix
    for(uint32_t row = 1; row < fm.n_rows; row++) {
        for(uint32_t sl = 1; sl < fm.n_cols - 1; sl++) {

            // cell indices
            //uint32_t c = cell(matrix, row, col);
            //uint32_t event_i = e_start + (row - 1) * state.stride;
            //uint32_t kmer_idx = k_start + col - 1;

            // Sum over states for previous row
            double sum = -INFINITY;

            // Only look back as far as the first state that can jump here
            uint32_t first_possible_state = 0;
            if(sl >= KHMM_MAX_JUMP)
                first_possible_state = sl - KHMM_MAX_JUMP;

            for(uint32_t sk = first_possible_state; sk <= sl; sk++) {

                // transition probability from state k to state l
                double t_kl = get(tm, sk, sl);
                double fm_k = get(fm, row - 1, sk);
                sum = log(exp(sum) + exp(t_kl + fm_k));
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf f: %.2lf s: %.2lf\n", row, sl, sk, t_kl, fm_k, sum);
#endif
            }

            // Emission probability for event i in state sl
            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t kmer_idx = sl - 1;

            uint32_t rank = !state.rc ? 
                kmer_rank(sequence + kmer_idx, K) : 
                rc_kmer_rank(sequence + kmer_idx, K);
            double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);
            
            set(fm, row, sl, lp_e + sum);

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) ei: %zu ki: %zu\n", row, sl, event_idx, kmer_idx);
            printf("(%d %d) sum: %.2lf lp_e: %.2lf fm: %.2lf\n", row, sl, sum, lp_e, get(fm, row, sl));
#endif
        }
    }

    // terminate by summing the last row and transitioning to end state
    double sum = -INFINITY;
    uint32_t tcol = fm.n_cols - 1;
    uint32_t lrow = fm.n_rows - 1;
    for(uint32_t sk = 0; sk < fm.n_cols - 1; sk++) {

        // transition probability from state k to state l
        double t_kl = get(tm, sk, tcol);
        double fm_k = get(fm, lrow, sk);
        sum = log(exp(sum) + exp(t_kl + fm_k));
    }
    return sum;
}

void initialize_backward_khmm(DoubleMatrix& bm, const DoubleMatrix& tm)
{
    // initialize forward calculation
    uint32_t tcol = tm.n_cols - 1;
    uint32_t row = bm.n_rows - 1;

    for(uint32_t si = 0; si < bm.n_cols; si++)
        set(bm, row, si, get(tm, si, tcol));
}

void fill_backward_khmm(DoubleMatrix& bm, // backward matrix
                         const DoubleMatrix& tm, //transitions
                         const char* sequence,
                         const HMMConsReadState& state,
                         uint32_t e_start)
{
    // Fill in matrix
    for(uint32_t row = bm.n_rows - 2; row > 0; row--) {
        for(uint32_t sk = 1; sk < bm.n_cols - 1; sk++) {

            // Sum over states for next row
            double sum = -INFINITY;
            for(uint32_t sl = 1; sl < bm.n_cols - 1; sl++) {

                // transition probability from state k to state l
                double t_kl = get(tm, sk, sl);
                double bm_l = get(bm, row + 1, sl);

                // Emit E_(i+1) in state sl
                uint32_t event_idx = e_start + row * state.stride; // for i + 1
                uint32_t kmer_idx = sl - 1;
                uint32_t rank = !state.rc ? 
                    kmer_rank(sequence + kmer_idx, K) : 
                    rc_kmer_rank(sequence + kmer_idx, K);
                double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);

                sum = log(exp(sum) + exp(lp_e + t_kl + bm_l));
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf f: %.2lf e: %.2lf s: %.2lf\n", row, sk, sl, t_kl, fm_l, lp_e, sum);
#endif
            }
            
            set(bm, row, sk, sum);

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) bm: %.2lf\n", row, sl, sum, lp_e, get(bm, row, sl));
#endif
        }
    }
}

double score_khmm_model(const std::string& consensus, const HMMConsReadState& state, AlignmentPolicy policy)
{
    uint32_t n_kmers = consensus.size() - K + 1;
    uint32_t n_states = n_kmers + 2; // one start and one end state

    DoubleMatrix tm;
    allocate_matrix(tm, n_states, n_states);

    fill_khmm_transitions(tm, consensus, state);
    
    uint32_t e_start = state.event_start_idx;
    uint32_t e_end = state.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate a matrix to hold the HMM result
    DoubleMatrix fm;
    allocate_matrix(fm, n_rows, n_states);

    initialize_forward_khmm(fm);
    fill_forward_khmm(fm, tm, consensus.c_str(), state, e_start);

    double score = 0.0f;
    if(policy == AP_GLOBAL) {
        // score by the bottom-right cell
        uint32_t last_row = fm.n_rows - 1;
        score = forward_khmm_terminate(fm, tm, last_row);
    } else if(policy == AP_SEMI_KMER) {

        // score by the best cell in the last column
        double best_score = -INFINITY;
        uint32_t best_row = 0;
        for(size_t row = 1; row < fm.n_rows - 1; ++row) {
            double s = forward_khmm_terminate(fm, tm, row);
            if(s > best_score) {
                best_score = s;
                best_row = row;
            }
        }

        score = best_score;
    } else {
        assert(false);
    }

    free_matrix(tm);
    free_matrix(fm);
    return score;
}

void debug_khmm_model(const std::string& name, const std::string& consensus, const HMMConsReadState& state)
{
    uint32_t n_kmers = consensus.size() - K + 1;
    uint32_t n_states = n_kmers + 2; // one start and one end state

    DoubleMatrix tm;
    allocate_matrix(tm, n_states, n_states);

    fill_khmm_transitions(tm, consensus, state);
    
    uint32_t e_start = state.event_start_idx;
    uint32_t e_end = state.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate and compute forward matrix
    DoubleMatrix fm;
    allocate_matrix(fm, n_rows, n_states);

    initialize_forward_khmm(fm);
    double lf = fill_forward_khmm(fm, tm, consensus.c_str(), state, e_start);

    // Allocate and compute backward matrix
    DoubleMatrix bm;
    allocate_matrix(bm, n_rows, n_states);

    initialize_backward_khmm(bm, tm);
    fill_backward_khmm(bm, tm, consensus.c_str(), state, e_start);

    //print_matrix(fm);
    //print_matrix(bm);

    // posterior decode
    std::vector<uint32_t> matches;
    std::vector<double> posteriors;
    std::vector<double> lpfm;
    
    uint32_t row = fm.n_rows - 1;
    uint32_t col = fm.n_cols - 2;
    while(row > 0) {

        // Calculate posterior probability that e_i is matched to k_j
        double max_posterior = -INFINITY;
        uint32_t max_s = 0;

        for(uint32_t si = 1; si < fm.n_cols - 1; ++si) {
            double lp = get(fm, row, si) + get(bm, row, si) - lf;
            if(lp > max_posterior) {
                max_posterior = lp;
                max_s = si;
            }
        }
        
        matches.push_back(max_s - 1); // kmer index
        posteriors.push_back(max_posterior);
        lpfm.push_back(get(fm, row, max_s));
        row -= 1;
    }

    std::reverse(matches.begin(), matches.end());
    std::reverse(posteriors.begin(), posteriors.end());
    std::reverse(lpfm.begin(), lpfm.end());
    
    printf("%s align max value: %.2lf\n", name.c_str(), lf);

    const CPoreModel& pm = state.read->pore_model[state.strand];
    uint32_t prev_ki = -1;

    // summary states
    double sum_emissions = 0.0f;
    uint32_t n_matches = 0;
    uint32_t n_skips = 0;
    uint32_t n_merges = 0;

    for(size_t i = 0; i < matches.size(); ++i) {

        uint32_t ei = e_start + i * state.stride;
        uint32_t ki = matches[i];
        
        // Emit kmer skip symbols
        char s;
        while(prev_ki != -1 && prev_ki + 1 < ki) {
            std::string kmer = consensus.substr(prev_ki + 1, K);
            
            uint32_t rank = !state.rc ? 
                kmer_rank(kmer.c_str(), K) : 
                rc_kmer_rank(kmer.c_str(), K);
        
            double model_m = (pm.state[rank].mean + pm.shift) * pm.scale;
            double model_s = pm.state[rank].sd * pm.scale;
            s = 'K';
            printf("%c\t%d\t%d\t0.0\t%s\t0.0\t%.2lf\t%.2lf\n", s, -1, prev_ki + 1, kmer.c_str(), model_m, lpfm[i]);
            n_skips += 1;
            prev_ki++;
        }

        if(prev_ki != -1 && ki - prev_ki == 0)
            s = 'E';
        else
            s = 'M';
        
        double level = state.read->events[state.strand].level[ei];
        double sd = state.read->events[state.strand].stdv[ei];
        
        uint32_t rank = !state.rc ? 
            kmer_rank(consensus.c_str() + ki, K) : 
            rc_kmer_rank(consensus.c_str() + ki, K);
        
        const CPoreModel& pm = state.read->pore_model[state.strand];
        double model_m = (pm.state[rank].mean + pm.shift) * pm.scale;
        double model_s = pm.state[rank].sd * pm.scale;
        double norm_level = (level - model_m) / model_s;
        std::string kmer = consensus.substr(ki, K);
        printf("%c\t%d\t%d\t%.2lf\t%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n", s, ei, ki, exp(posteriors[i]), kmer.c_str(), level, model_m, norm_level, lpfm[i]);

        // summary
        n_matches += (s == 'M');
        n_merges += (s == 'E');
        if(s != 'K')
            sum_emissions += log_probability_match(*state.read, rank, ei, state.strand);
        prev_ki = ki;
    }

    printf("Summary: sum_emission: %.2lf M: %zu E: %zu K: %zu\n", sum_emissions, n_matches, n_merges, n_skips);
    //print_matrix(fm);
//    printf("%s %.2lf\n", consensus.c_str(), lf);

    free_matrix(tm);
    free_matrix(fm);
    free_matrix(bm);
}

double score_emission_dp(const std::string& sequence, const HMMConsReadState& state)
{
    uint32_t n_kmers = sequence.size() - K + 1;
    uint32_t n_cols = n_kmers + 1;

    uint32_t e_start = state.event_start_idx;
    uint32_t e_end = state.event_stop_idx;
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;

    uint32_t n_rows = n_events + 1;

    // Allocate a matrix to hold the HMM result
    DoubleMatrix m;
    allocate_matrix(m, n_rows, n_cols);

    // Initialize matrix to -INFINITY except for (0,0)
    for(uint32_t row = 0; row < m.n_rows; ++row) {
        for(uint32_t col = 0; col < m.n_cols; ++col) {
            set(m, row, col, -INFINITY);
        }
    }
    set(m, 0, 0, 0.0f);

    // Fill matrix
    for(uint32_t row = 1; row < m.n_rows; ++row) {
        for(uint32_t col = 1; col < m.n_cols; ++col) {
            
            // Get the emission probability for matching this event to this kmer
            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t kmer_idx = col - 1;

            uint32_t rank = !state.rc ? 
                kmer_rank(sequence.c_str() + kmer_idx, K) : 
                rc_kmer_rank(sequence.c_str() + kmer_idx, K);

            double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);
            
            // Get the scores for the 3 possible preceding cells           
            double up = lp_e + get(m, row - 1, col);
            double diag = lp_e + get(m, row - 1, col - 1);
            double left = get(m, row, col - 1);
            double best = max3(up, diag, left);

            set(m, row, col, best);
        }
    }

    double score = get(m, m.n_rows - 1, m.n_cols - 1);
    free_matrix(m);
    return score;
}

// Handy wrappers for scoring/debugging functions
// The consensus algorithms call into these so we can switch
// scoring functinos without writing a bunch of code
double score_sequence(const std::string& sequence, const HMMConsReadState& state)
{
    return score_khmm_model(sequence, state, AP_GLOBAL);
    //return score_emission_dp(sequence, state);
}

void debug_sequence(const std::string& name, const std::string& sequence, const HMMConsReadState& state)
{
    return debug_khmm_model(name, sequence, state);
}


extern "C"
void run_debug()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    for(uint32_t i = 0; i < g_data.read_states.size(); ++i) {
        debug_sequence("input",  "AACAGTCCACTATTGGATGGTAAAGCCAACAGAAATTTTTACGCAAG", g_data.read_states[i]);
        debug_sequence("oldbad", "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAATTTACGCAAG", g_data.read_states[i]);
        debug_sequence("final", "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTACGCAAG", g_data.read_states[i]);
        debug_sequence("truth", "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTTACGCAAG", g_data.read_states[i]);
    }
}

struct PathCons
{
    std::string path;
    double score;
    std::string mutdesc;
    
};
typedef std::vector<PathCons> PathConsVector;

bool sortPathConsAsc(const PathCons& a, const PathCons& b)
{
    return a.score > b.score;
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
                PathCons ps = { ns, 0.0f };
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
        PathCons pc = { sequence, 0.0f };
        mutations.push_back(pc);
    }

    // Mutate every base except for in the first/last k-mer
    for(size_t si = K; si < sequence.size() - K; ++si) {
        
        // All subs
        for(size_t bi = 0; bi < 4; bi++) {
            char b = "ACGT"[bi];
            if(sequence[si] == b)
                continue;
            PathCons pc = { sequence, 0.0f };
            pc.path[si] = b;
            std::stringstream ss;
            ss << "sub-" << si << "-" << b;
            pc.mutdesc = ss.str();
            mutations.push_back(pc);


        }

        // 1bp del at this position
        {
            PathCons pc = { sequence, 0.0f };
            pc.path.erase(si, 1);
            
            std::stringstream ss;
            ss << "del-" << si;
            pc.mutdesc = ss.str();
            
            mutations.push_back(pc);
        }

        // All 1bp ins before this position
        for(size_t bi = 0; bi < 4; bi++) {
            char b = "ACGT"[bi];
            PathCons pc = { sequence, 0.0f };
            pc.path.insert(si, 1, b);
            
            std::stringstream ss;
            ss << "ins-" << si << "-" << b;
            pc.mutdesc = ss.str();
            
            mutations.push_back(pc);
        }
    }

    return mutations;
}

extern "C"
void run_mutation()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    std::string sequence = g_data.candidate_consensus.front();
    printf("Initial consensus: %s\n", sequence.c_str());
    int iteration = 0;
    while(iteration++ < 10) {

        // Generate possible sequences
        PathConsVector paths = generate_mutations(sequence);

        // Score all reads
        for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
            std::vector<double> scores;
            double sum_score = -INFINITY;

            // Score all paths
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
                sum_score = log(exp(sum_score) + exp(curr));
                scores.push_back(curr);
            }
            
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                paths[pi].score += (scores[pi] - scores[0]);
            }
        }
        
        // Cull paths
        std::sort(paths.begin(), paths.end(), sortPathConsAsc);
        std::string truth = "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTTACGCAAGCTAAAGCCCGGCAGATGATTATCTGTGCCGATATGATCAAACCGCGGTTGAATGAAAC";
        
        printf("Iteration %d\n", iteration);
        for(size_t pi = 0; pi < paths.size(); ++pi) {

            // Calculate the length of the matching prefix with the truth
            const std::string& s = paths[pi].path;

            uint32_t plen = 0;
            while(s[plen] == truth[plen] && plen < s.length() && plen < truth.length())
                plen++;

            // Match info
            char match = plen == s.length() ? '*' : ' ';
            
            char initial = s == sequence ? 'I' : ' ';

            printf("%zu\t%s\t%.1lf %d %c %c %s", pi, paths[pi].path.c_str(), paths[pi].score, plen, match, initial, paths[pi].mutdesc.c_str());
            // If this is the truth path or the best path, show the scores for all reads
            if(pi == 0 || match == '*') {
                for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
                    double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
                    printf("%.1lf ", curr);
                }
            }
            printf("\n");
        }
        
        // check if no improvement was made
        if(paths[0].path == sequence)
            break;
        sequence = paths[0].path;
    }
    
    g_data.consensus_result = sequence;
}

PathConsVector generate_rewrites(const std::string& sequence, uint32_t position)
{
    PathConsVector rewrites;
    if(position >= sequence.size() - 5)
        return rewrites;

    // Add the unmutated sequence
    {
        PathCons pc = { sequence, 0.0f };
        rewrites.push_back(pc);
    }

    std::string prefix = sequence.substr(0, position);
    std::string suffix = sequence.substr(position);

    printf("REWRITE -- %zu %s %s\n", position, prefix.c_str(), suffix.c_str());

    // Insert every possible 5-mer at the stated position
    std::string first(K, 'A'); // AAAAA
    std::string kmer = first;

    do {
        for(uint32_t del = 0; del < 5; ++del) {
            std::string ns = prefix + kmer + suffix.substr(del);
            PathCons ps = { ns, 0.0f };
            rewrites.push_back(ps);
        }
        lexicographic_next(kmer);
    } while(kmer != first);

    return rewrites;
}

extern "C"
void run_rewrite()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    std::string sequence = g_data.candidate_consensus.front();
    printf("Initial consensus: %s\n", sequence.c_str());

    uint32_t position = 5;

    while(1) {

        printf("Position %d\n", position);
        
        // Generate possible sequences
        PathConsVector paths = generate_rewrites(sequence, position++);

        if(paths.size() == 0)
            break;

        // Score all reads
        for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
            std::vector<double> scores;
            double sum_score = -INFINITY;

            // Score all paths
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
                sum_score = log(exp(sum_score) + exp(curr));
                scores.push_back(curr);
            }
            
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                paths[pi].score += (scores[pi] - scores[0]);
            }
        }
        
        // select new sequence
        std::sort(paths.begin(), paths.end(), sortPathConsAsc);
        
        for(size_t pi = 0; pi < paths.size(); ++pi) {

            // Calculate the length of the matching prefix with the truth
            const std::string& s = paths[pi].path;

            char initial = s == sequence ? 'I' : ' ';

            printf("%zu\t%s\t%.1lf %c %s", pi, paths[pi].path.c_str(), paths[pi].score, initial, paths[pi].mutdesc.c_str());
            // If this is the truth path or the best path, show the scores for all reads
            if(pi == 0) {
                for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
                    double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
                    printf("%.1lf ", curr);
                }
            }
            printf("\n");
        }
        
        sequence = paths[0].path;
    }
    
    g_data.consensus_result = sequence;
}

extern "C"
void run_selection()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    PathConsVector paths;

    for(uint32_t ci = 0; ci < g_data.candidate_consensus.size(); ++ci) {
        PathCons pc = { g_data.candidate_consensus[ci], 0.0f };
        paths.push_back(pc);

    }

    // Score all reads
    for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
        std::vector<double> scores;
        double sum_score = -INFINITY;

        // Score all paths
        for(size_t pi = 0; pi < paths.size(); ++pi) {
            double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
            debug_sequence(paths[pi].path, paths[pi].path, g_data.read_states[ri]);
            sum_score = log(exp(sum_score) + exp(curr));
            scores.push_back(curr);
        }

        for(size_t pi = 0; pi < paths.size(); ++pi) {
            paths[pi].score += scores[pi];
        }
    }
        
    std::sort(paths.begin(), paths.end(), sortPathConsAsc);
    for(size_t pi = 0; pi < paths.size(); ++pi) {
        printf("%zu\t%s\t%.1lf\n", pi, paths[pi].path.c_str(), paths[pi].score);
    }

    g_data.consensus_result = paths[0].path;
}

extern "C"
void run_consensus()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    std::string consensus = "AACAG";
#if 0
    // Populate initial pathcons vector
    PathConsVector paths;
    PathCons pc = { consensus, 0.0f };
    paths.push_back( pc );

    for(size_t i = 0; i < paths.size(); ++i) {
        printf("%zu %s\n", i, paths[i].path.c_str());
    }
    
    int iteration = 0;
    while(iteration++ < 10) {
        
        extend_paths(paths);

        // Score all reads
        for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
            std::vector<double> scores;
            double sum_score = -INFINITY;

            // Score all paths
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                double curr = score_khmm_model(paths[pi].path, g_data.read_states[ri], AP_SEMI_KMER);
                sum_score = log(exp(sum_score) + exp(curr));
                scores.push_back(curr);
            }
            
            for(size_t pi = 0; pi < paths.size(); ++pi) {
                paths[pi].score += (scores[pi] - sum_score);
            }
        }
        
        // Cull paths
        std::sort(paths.begin(), paths.end(), sortPathConsAsc);
        std::string truth = "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTTACGCAAGCTAAAGCCCGGCAGATGATTATCTGTGCCGATATGATCAAACCGCGGTTGAATGAAAC";
        
        printf("Iteration %d\n", iteration);
        for(size_t pi = 0; pi < paths.size(); ++pi) {

            // Calculate the length of the matching prefix with the truth
            const std::string& s = paths[pi].path;

            uint32_t plen = 0;
            while(s[plen] == truth[plen] && plen < s.length() && plen < truth.length())
                plen++;

            // Match info
            char match = plen == s.length() ? '*' : ' ';
            
            printf("%zu %s %.1lf %d %c", pi, paths[pi].path.c_str(), paths[pi].score, plen, match);
            // If this is the truth path or the best path, show the scores for all reads
            if(pi == 0 || match == '*') {
                for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
                    double curr = score_khmm_model(paths[pi].path, g_data.read_states[ri], AP_SEMI_KMER);
                    printf("%.1lf ", curr);
                }
            }
            printf("\n");
        }

        paths.resize(std::min(paths.size(), (size_t)256));
    }
#endif

    // Populate initial pathcons vector
    PathConsVector paths;
    std::string extension = "AAAAA";
    do {
        PathCons ps = { consensus + extension, 0.0f };
        paths.push_back(ps);
        lexicographic_next(extension);
    } while(extension != "AAAAA");

    uint32_t window = 20;
    int iteration = 0;
    while(iteration++ < 30) {

        // Extend paths
        PathConsVector new_paths;
        for(size_t i = 0; i < paths.size(); ++i) {
            for(size_t b = 0; b < 4; ++b) {
                PathCons pc = { paths[i].path + "ACGT"[b], 0.0f };
                new_paths.push_back(pc);
            }
        }

        // Score all reads
        for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
            std::vector<double> scores;
            double sum_score = -INFINITY;

            // Score all paths
            for(size_t pi = 0; pi < new_paths.size(); ++pi) {
                double curr = score_khmm_model(new_paths[pi].path, g_data.read_states[ri], AP_SEMI_KMER);
                sum_score = log(exp(sum_score) + exp(curr));
                scores.push_back(curr);
            }
            
            for(size_t pi = 0; pi < new_paths.size(); ++pi) {
                new_paths[pi].score += (scores[pi] - sum_score);
            }
        }
        
        // Cull paths
        std::sort(new_paths.begin(), new_paths.end(), sortPathConsAsc);
        std::string truth = "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTTACGCAAGCTAAAGCCCGGCAGATGATTATCTGTGCCGATATGATCAAACCGCGGTTGAATGAAAC";
        
        printf("Iteration %d\n", iteration);
        for(size_t pi = 0; pi < new_paths.size(); ++pi) {

            // Calculate the length of the matching prefix with the truth
            const std::string& s = new_paths[pi].path;

            uint32_t plen = 0;
            while(s[plen] == truth[plen] && plen < s.length() && plen < truth.length())
                plen++;

            // Match info
            char match = plen == s.length() ? '*' : ' ';
            
            printf("%zu %s %.1lf %d %c", pi, new_paths[pi].path.c_str(), new_paths[pi].score, plen, match);
            // If this is the truth path or the best path, show the scores for all reads
            if(pi == 0 || match == '*') {
                for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
                    double curr = score_khmm_model(new_paths[pi].path, g_data.read_states[ri], AP_SEMI_KMER);
                    printf("%.1lf ", curr);
                }
            }
            printf("\n");
        }

        paths.assign(new_paths.begin(), new_paths.begin() + 1024);
    }
}

int main(int argc, char** argv)
{

}
