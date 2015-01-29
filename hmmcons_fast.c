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
#include "hmmcons_khmm_parameters.h"
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
const static double P_RANDOM_SKIP = 0.05;
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
#define PRINT_TRAINING_MESSAGES

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

    // one event sequence for each strand
    CEventSequence events[2];

    // one set of parameters per strand
    KHMMParameters parameters[2];
};

double get_drift_corrected_level(const CSquiggleRead& read, uint32_t event_idx, uint32_t strand)
{
    double level = read.events[strand].level[event_idx];
    // correct level by drift
    double start = read.events[strand].time[0];
    double time = read.events[strand].time[event_idx] - start;
    return level - (time * read.pore_model[strand].drift);
}

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

struct PosteriorState
{
    uint32_t event_idx;
    uint32_t kmer_idx;
    double l_posterior;
    double l_fm;
    char state;
};

// A global vector used to store data we've received from the python code
struct HmmConsData
{
    //
    std::vector<CSquiggleRead> reads;

    //
    std::vector<HMMConsReadState> read_states;
    std::vector<std::string> candidate_consensus;
    std::string consensus_result;
};
HmmConsData g_data;
bool g_initialized = false;

// 
// Add the log-scaled values a and b using a transform to avoid
// precision errors
inline double add_logs(const double a, const double b)
{
    if(a == -INFINITY && b == -INFINITY)
        return -INFINITY;

    if(a > b) {
        double diff = b - a;
        return a + log(1.0 + exp(diff));
    } else {
        double diff = a - b;
        return b + log(1.0 + exp(diff));
    }
}

extern "C"
void initialize()
{
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
    sr.read_id = g_data.reads.size() - 1;

    for(uint32_t i = 0; i < NUM_STRANDS; ++i) {
        // Initialize pore model   
        sr.pore_model[i].scale = params.pore_model[i].scale;
        sr.pore_model[i].shift = params.pore_model[i].shift;
        sr.pore_model[i].drift = params.pore_model[i].drift;
        sr.pore_model[i].var = params.pore_model[i].var;
        
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
                                                                 sr.pore_model[i].state[0].mean, 
                                                                 sr.pore_model[i].state[0].sd);
    
        printf("First 100 events of %d\n", sr.events[i].n_events);
        for(int j = 0; j < 100; ++j)
            printf("%d: %lf\n", j, sr.events[i].level[j]);
        */
    }

    // Initialize hmm parameters for both strands of the read
    khmm_parameters_initialize(sr.parameters[0]);
    khmm_parameters_initialize(sr.parameters[1]);
}

extern "C"
void add_candidate_consensus(char* str)
{
    g_data.candidate_consensus.push_back(str);
}

// Make a unique index for the strand this read state represents
uint32_t get_strand_idx(const HMMConsReadState& rs)
{
    return rs.read->read_id + rs.strand;
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

extern "C"
void add_read_state(CReadStateInterface params)
{
    // add read state
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
// Template Matrix for POD types
//
template<typename T>
struct Matrix
{
    T* cells;
    uint32_t n_rows;
    uint32_t n_cols;
};

typedef Matrix<double> DoubleMatrix;
typedef Matrix<uint32_t> UInt32Matrix;

template<typename T>
void allocate_matrix(Matrix<T>& matrix, uint32_t n_rows, uint32_t n_cols)
{
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    
    uint32_t N = matrix.n_rows * matrix.n_cols;
    matrix.cells = (T*)malloc(N * sizeof(T));
    memset(matrix.cells, 0, N * sizeof(T));
}

template<typename T>
void free_matrix(Matrix<T>& matrix)
{
    assert(matrix.cells != NULL);
    free(matrix.cells);
    matrix.cells = NULL;
}

template<typename T>
inline uint32_t cell(const Matrix<T>& matrix, uint32_t row, uint32_t col)
{
    return row * matrix.n_cols + col;
}

template<typename T, typename U>
inline void set(Matrix<T>& matrix, uint32_t row, uint32_t col, U v)
{
    uint32_t c = cell(matrix, row, col);
    matrix.cells[c] = v;
}

template<typename T>
inline T get(const Matrix<T>& matrix, uint32_t row, uint32_t col)
{
    uint32_t c = cell(matrix, row, col);
    return matrix.cells[c];
}

//
//
//
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

// wrapper to get the rank for a kmer on the right strand
inline uint32_t get_rank(const HMMConsReadState& state, const char* s, uint32_t ki)
{
    const char* p = s + ki;
    return !state.rc ?  kmer_rank(p, K) : rc_kmer_rank(p, K);
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
    double level = get_drift_corrected_level(read, event_idx, strand);

    double m = pm.state[kmer_rank].level_mean * pm.scale + pm.shift;
    double s = pm.state[kmer_rank].level_stdv * pm.var;
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
    const KHMMParameters& parameters = state.read->parameters[state.strand];

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
                p_i_j = parameters.self_transition;
            } else {
        
                uint32_t rank_i = get_rank(state, consensus.c_str(), ki);
                uint32_t rank_j = get_rank(state, consensus.c_str(), kj);

                double level_i = (pm.state[rank_i].level_mean + pm.shift) * pm.scale;
                double level_j = (pm.state[rank_j].level_mean + pm.shift) * pm.scale;
                
                double p_skip = get_skip_probability(parameters, level_i, level_j);
                p_i_j = (1 - sum) * (1 - p_skip);
                assert(p_i_j >= 0.0f && p_i_j <= 1.0f);

#ifdef DEBUG_TRANSITION
                printf("\t\t%zu -> %zu %.2lf %.2lf p_skip: %.4lf p: %.2lf\n", ki, kj, level_i, level_j, p_skip, p_i_j);
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
        sum = add_logs(sum, t_kl + fm_k);
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
                sum = add_logs(sum, t_kl + fm_k);
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf f: %.2lf s: %.2lf\n", row, sl, sk, t_kl, fm_k, sum);
#endif
            }

            // Emission probability for event i in state sl
            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t kmer_idx = sl - 1;
            uint32_t rank = get_rank(state, sequence, kmer_idx);
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
        sum = add_logs(sum, t_kl + fm_k);
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
                uint32_t rank = get_rank(state, sequence, kmer_idx);
                double lp_e = log_probability_match(*state.read, rank, event_idx, state.strand);

                sum = add_logs(sum, lp_e + t_kl + bm_l);
#ifdef DEBUG_HMM_UPDATE
                printf("\t(%d %d %d) t: %.2lf b: %.2lf e: %.2lf s: %.2lf\n", row, sk, sl, t_kl, bm_l, lp_e, sum);
#endif
            }
            
            set(bm, row, sk, sum);

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) bm: %.2lf\n", row, sk, get(bm, row, sk));
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

std::vector<PosteriorState> posterior_decode_khmm(const std::string& sequence, const HMMConsReadState& state)
{
    uint32_t n_kmers = sequence.size() - K + 1;
    uint32_t n_states = n_kmers + 2; // one start and one end state

    DoubleMatrix tm;
    allocate_matrix(tm, n_states, n_states);

    fill_khmm_transitions(tm, sequence, state);
    
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
    double lf = fill_forward_khmm(fm, tm, sequence.c_str(), state, e_start);

    // Allocate and compute backward matrix
    DoubleMatrix bm;
    allocate_matrix(bm, n_rows, n_states);

    initialize_backward_khmm(bm, tm);
    fill_backward_khmm(bm, tm, sequence.c_str(), state, e_start);

    // posterior decode
    std::vector<PosteriorState> output;
    
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
    
        uint32_t event_idx = e_start + (row - 1) * state.stride;
        uint32_t kmer_idx = max_s - 1;
 
        double lpfm = get(fm, row, max_s);

        PosteriorState ps = { event_idx, kmer_idx, max_posterior, lpfm, 'N' };
        output.push_back(ps);
        row -= 1;
    }

    std::reverse(output.begin(), output.end());

    // First state is always a match
    output[0].state = 'M';
    uint32_t prev_ei = output[0].event_idx;
    uint32_t prev_ki = output[0].kmer_idx;

    for(uint32_t pi = 1; pi < output.size(); ++pi) {
        uint32_t ei = output[pi].event_idx;
        uint32_t ki = output[pi].kmer_idx;

        assert(abs(ei - prev_ei) == 1);

        if(ki == prev_ki) {
            output[pi].state = 'E';
        } else if(ki - prev_ki == 1) {
            output[pi].state = 'M';
        } else {
            assert(ki - prev_ki > 1);
            output[pi].state = 'K';
        }

        prev_ei = ei;
        prev_ki = ki;
    }

    free_matrix(tm);
    free_matrix(fm);
    free_matrix(bm);
    return output;
}

void update_training_khmm(const std::string& consensus, 
                          const HMMConsReadState& state)
{
    std::vector<PosteriorState> pstates = posterior_decode_khmm(consensus, state);

    const CPoreModel& pm = state.read->pore_model[state.strand];
    TrainingData& training_data = state.read->parameters[state.strand].training_data;
    size_t n_kmers = consensus.size() - K + 1;
    uint32_t strand_idx = get_strand_idx(state);

    for(size_t pi = 0; pi < pstates.size(); ++pi) {

        uint32_t ei = pstates[pi].event_idx;
        uint32_t ki = pstates[pi].kmer_idx;
        char s = pstates[pi].state;
    
        // Record transition observations
        // We do not record observations for merge states as there was no kmer transitions
        // We also do not record observations for the beginning of the matches as the
        // alignment may be poor due to edge effects
        if(pi > 5 && pi < pstates.size() - 5) {
 
            // transition           
            if(s != 'E') {
                uint32_t transition_kmer_from = pstates[pi - 1].kmer_idx;
                uint32_t transition_kmer_to = pstates[pi].kmer_idx;

                // Specially handle skips
                // We only want to record the first k-mer skipped if multiple were skipped
                if(s == 'K') {
                    transition_kmer_from = pstates[pi - 1].kmer_idx;
                    transition_kmer_to = transition_kmer_from + 1;
                }
                
                assert(transition_kmer_from < n_kmers && transition_kmer_to < n_kmers);

                uint32_t rank1 = get_rank(state, consensus.c_str(), transition_kmer_from);
                uint32_t rank2 = get_rank(state, consensus.c_str(), transition_kmer_to);
            
                double ke1 = (pm.state[rank1].level_mean + pm.shift) * pm.scale;
                double ke2 = (pm.state[rank2].level_mean + pm.shift) * pm.scale;

#ifdef PRINT_TRAINING_MESSAGES
                printf("TRAIN_SKIP\t%zu\t%.3lf\t%.3lf\t%c\n", strand_idx, ke1, ke2, s);
#endif
                TransitionObservation to = { ke1, ke2, s };
                training_data.transitions.push_back(to);
            }

            // emission
            double level = get_drift_corrected_level(*state.read, ei, state.strand);
            double sd = state.read->events[state.strand].stdv[ei];
            double start_time = state.read->events[state.strand].time[ei];
            double end_time = state.read->events[state.strand].time[ei + 1];
            if(ki >= n_kmers)
                printf("%zu %zu %zu %zu %.2lf %c\n", pi, ei, ki, n_kmers, pstates[pi].l_fm, s);
            
            assert(ki < n_kmers);
            uint32_t rank = get_rank(state, consensus.c_str(), ki);
        
            double model_m = (pm.state[rank].level_mean + pm.shift) * pm.scale;
            double model_s = pm.state[rank].level_stdv * pm.scale;
            double norm_level = (level - model_m) / model_s;

            if(s == 'M')
                training_data.emissions_for_matches.push_back(norm_level);

#ifdef PRINT_TRAINING_MESSAGES
            printf("TRAIN_EMISSION\t%zu\t%zu\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%c\n", strand_idx, ei, level, sd, model_m, model_s, norm_level, end_time - start_time, s);
#endif
        }

        // summary
        training_data.n_matches += (s == 'M');
        training_data.n_merges += (s == 'E');
        training_data.n_skips += (s == 'K');
    }
}

void debug_khmm_model(const std::string& name,
                      uint32_t id,
                      const std::string& consensus, 
                      const HMMConsReadState& state)
{
    std::vector<PosteriorState> pstates = posterior_decode_khmm(consensus, state);

    for(size_t pi = 0; pi < pstates.size(); ++pi) {

        uint32_t ei = pstates[pi].event_idx;
        uint32_t ki = pstates[pi].kmer_idx;
        char s = pstates[pi].state;
    
        double level = get_drift_corrected_level(*state.read, ei, state.strand);
        double sd = state.read->events[state.strand].stdv[ei];
        double duration = state.read->events[state.strand].time[ei];
        uint32_t rank = get_rank(state, consensus.c_str(), ki);
        
        const CPoreModel& pm = state.read->pore_model[state.strand];
        double model_m = (pm.state[rank].level_mean + pm.shift) * pm.scale;
        double model_s = pm.state[rank].level_stdv * pm.scale;
        double norm_level = (level - model_m) / model_s;
        
        double model_sd_mean = pm.state[rank].sd_mean;
        double model_sd_stdv = pm.state[rank].sd_mean;

        std::string kmer = consensus.substr(ki, K);
 
        printf("DEBUG\t%s\t%d\t", name.c_str(), id);
        printf("%c\t%zu\t%zu\t", s, ei, ki);
        printf("%s\t%.3lf\t", kmer.c_str(), duration);
        printf("(%.1lf, %.1lf, %.1lf)\t", level, model_m, norm_level);
        printf("(%.1lf, %.1lf, %.1lf)\t", sd, model_sd_mean, (sd - model_sd_mean) / model_sd_stdv);
        printf("%.2lf\t%.2lf\n", exp(pstates[pi].l_posterior), pstates[pi].l_fm);
    }
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
            uint32_t rank = get_rank(state, sequence.c_str(), kmer_idx);
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

// The indices of a k-mer match in a pair of sequences
struct kLCSPair
{
    uint32_t i;
    uint32_t j;
};
typedef std::vector<kLCSPair> kLCSResult;

// Helper to backtrack through the kLCS matrix
void _kLCSBacktrack(const UInt32Matrix& m,
                    const std::string& a, 
                    const std::string& b,
                    uint32_t row,
                    uint32_t col,
                    kLCSResult& result)
{
    if(row == 0 || col == 0)
        return;

    const char* ka = a.c_str() + row - 1;
    const char* kb = b.c_str() + col - 1;

    if(strncmp(ka, kb, K) == 0) {
        kLCSPair p = { row - 1, col - 1 };
        result.push_back(p);
        return _kLCSBacktrack(m, a, b, row - 1, col - 1, result);
    } else {

        if(get(m, row - 1, col) > get(m, row, col - 1)) {
            return _kLCSBacktrack(m, a, b, row - 1, col, result);
        } else {
            return _kLCSBacktrack(m, a, b, row, col - 1, result);
        }
    }
}

// Return the longest common subseuqence of k-mers between the two strings
kLCSResult kLCS(const std::string& a, const std::string& b)
{
    uint32_t n_kmers_a = a.size() - K + 1;
    uint32_t n_kmers_b = a.size() - K + 1;

    uint32_t n_rows = n_kmers_a + 1;
    uint32_t n_cols = n_kmers_b + 1;


    UInt32Matrix m;
    allocate_matrix(m, n_rows, n_cols);

    // Initialize first row/col to zero
    for(uint32_t row = 0; row < m.n_rows; ++row)
        set(m, row, 0, 0);
    for(uint32_t col = 0; col < m.n_cols; ++col)
        set(m, 0, col, 0);
    
    // Fill matrix
    for(uint32_t row = 1; row < m.n_rows; ++row) {
        for(uint32_t col = 1; col < m.n_cols; ++col) {
    
            const char* ka = a.c_str() + row - 1;
            const char* kb = b.c_str() + col - 1;

            uint32_t score = 0;
            if(strncmp(ka, kb, K) == 0) {
                uint32_t diag = get(m, row - 1, col - 1);
                score = diag + 1;
            } else {
                uint32_t left = get(m, row, col - 1);
                uint32_t up = get(m, row - 1, col);
                score = std::max(left, up);
            }
            set(m, row, col, score);
        }
    }

    kLCSResult result;
    _kLCSBacktrack(m, a, b, n_rows - 1, n_cols -  1, result);

    // Backtrack appends from the end to the start, reverse the vector of matches
    std::reverse(result.begin(), result.end());
    free_matrix(m);
    return result;
}

// Handy wrappers for scoring/debugging functions
// The consensus algorithms call into these so we can switch
// scoring functinos without writing a bunch of code
double score_sequence(const std::string& sequence, const HMMConsReadState& state)
{
    return score_khmm_model(sequence, state, AP_GLOBAL);
    //return score_emission_dp(sequence, state);
}

void debug_sequence(const std::string& name, uint32_t id, const std::string& sequence, const HMMConsReadState& state)
{
    return debug_khmm_model(name, id, sequence, state);
}

std::vector<PosteriorState> posterior_decode(const std::string& sequence, const HMMConsReadState& state)
{
    return posterior_decode_khmm(sequence, state);
}

extern "C"
void run_debug()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    for(uint32_t i = 0; i < g_data.read_states.size(); ++i) {
        debug_sequence("input",  i, "AACAGTCCACTATTGGATGGTAAAGCCAACAGAAATTTTTACGCAAG", g_data.read_states[i]);
        debug_sequence("oldbad", i, "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAATTTACGCAAG", g_data.read_states[i]);
        debug_sequence("final", i, "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTACGCAAG", g_data.read_states[i]);
        debug_sequence("truth", i, "AACAGTCCACTATTGGATGGTAAAGCGCTAACAGAAATTTTTACGCAAG", g_data.read_states[i]);
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

// This scores each path using the HMM and 
// sorts the paths into ascending order by score
void score_paths(PathConsVector& paths)
{
    std::string first = paths[0].path;
    
    double MIN_FIT = INFINITY;

    // Score all reads
    for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
        const HMMConsReadState& read_state = g_data.read_states[ri];
        const KHMMParameters& parameters = read_state.read->parameters[read_state.strand];
 
        if( fabs(parameters.fit_quality) > MIN_FIT)
            continue;

        std::vector<double> scores;
        double sum_score = -INFINITY;

        // Score all paths
        for(size_t pi = 0; pi < paths.size(); ++pi) {
            double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
            sum_score = add_logs(sum_score, curr);
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

        char initial = s == first ? 'I' : ' ';

        printf("%zu\t%s\t%.1lf %c %s", pi, paths[pi].path.c_str(), paths[pi].score, initial, paths[pi].mutdesc.c_str());
        // If this is the truth path or the best path, show the scores for all reads
        if(pi == 0 || initial == 'I') {
            for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
                const HMMConsReadState& read_state = g_data.read_states[ri];
                const KHMMParameters& parameters = read_state.read->parameters[read_state.strand];
                if( fabs(parameters.fit_quality) > MIN_FIT)
                    continue;

                double curr = score_sequence(paths[pi].path, g_data.read_states[ri]);
                printf("%.1lf,%.2lf ", parameters.fit_quality, curr);
            }
        }
        printf("\n");
    }
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

        score_paths(paths);

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
        
        score_paths(paths);

        sequence = paths[0].path;
    }
    
    g_data.consensus_result = sequence;
}

extern "C"
void run_splice()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }
    
    // initialize 
    std::string base = g_data.candidate_consensus[0];

    uint32_t num_rounds = 6;
    uint32_t round = 0;
    while(round++ < num_rounds) {

        PathConsVector paths;
        PathCons base_path = { base, 0.0f };
        paths.push_back(base_path);

        for(uint32_t ci = 1; ci < g_data.candidate_consensus.size(); ++ci) {
            const std::string& read = g_data.candidate_consensus[ci];
            kLCSResult result = kLCS(base, read);

            uint32_t match_idx = 0;
            while(match_idx < result.size()) {
                uint32_t last_idx = result.size() - 1;

                // advance the match to the next point of divergence
                while(match_idx != last_idx && result[match_idx].i == result[match_idx + 1].i - 1)
                    match_idx++;

                // no more divergences to process
                if(match_idx == last_idx)
                    break;

                uint32_t bl = result[match_idx + 1].i - result[match_idx].i;
                uint32_t rl = result[match_idx + 1].j - result[match_idx].j;

                std::string base_subseq = base.substr(result[match_idx].i, bl);
                std::string read_subseq = read.substr(result[match_idx].j, rl);

                // Perform the splice
                PathCons new_path = { base, 0.0f };
                new_path.path.replace(result[match_idx].i, bl, read_subseq);
                paths.push_back(new_path);
                
                //printf("REPLACE: %s %s %s\n", base_subseq.c_str(), read_subseq.c_str(), new_path.path.c_str());

                match_idx += 1;
            }
        }

        score_paths(paths);
        
        for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {
            debug_sequence("best", ri, paths[0].path, g_data.read_states[ri]); 
            debug_sequence("truth", ri, base, g_data.read_states[ri]); 
        }

        if(paths[0].path == base)
            break;
        base = paths[0].path;
    
    }

    g_data.consensus_result = base;
}

// update the training data on the current segment
extern "C"
void learn_segment()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }
    
    // initialize 
    std::string segment = g_data.candidate_consensus[0];
     
    for(uint32_t ri = 0; ri < g_data.read_states.size(); ++ri) {

        std::vector<PosteriorState> decodes = posterior_decode(segment, g_data.read_states[ri]);
        //debug_sequence("training segment", ri, segment, g_data.read_states[ri]);
        update_training_khmm(segment, g_data.read_states[ri]);
    }
}

extern "C"
void train()
{
    for(uint32_t ri = 0; ri < g_data.reads.size(); ++ri) {
        khmm_parameters_train(g_data.reads[ri].parameters[0]);
        khmm_parameters_train(g_data.reads[ri].parameters[1]);
    }
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
            debug_sequence(paths[pi].path, ri, paths[pi].path, g_data.read_states[ri]);
            sum_score = add_logs(sum_score, curr);
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
                sum_score = add_logs(sum_score, curr);
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
                sum_score = add_logs(sum_score, curr);
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
