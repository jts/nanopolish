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
#include "hmmcons_poremodel.h"
#include "hmmcons_interface.h"

// Constants

// strands
const uint8_t T_IDX = 0;
const uint8_t C_IDX = 1;
const uint8_t NUM_STRANDS = 2;

// 
const uint8_t K = 5;

const static double LOG_KMER_INSERTION = log(0.1);

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

//#define DEBUG_HMM_UPDATE 1
//#define DEBUG_HMM_EMISSION 1

struct CEventSequence
{
    uint32_t n_events;
    const double* level;
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
    uint32_t event_idx;
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
    // transitions from the match state
    params.t[0][0] = 90.f; // M
    params.t[0][1] = 5.f;  // E
    params.t[0][2] = 5.f;  // K
    params.t[0][3] = 1.f;  // terminal
    
    // transitions from the event insertion state
    params.t[1][0] = 85.f; // M
    params.t[1][1] = 10.f;  // E
    params.t[1][2] = 5.f;  // K
    params.t[1][3] = 1.f;  // terminal

    // transitions from the k-mer insertion state
    params.t[2][0] = 85.f; // M
    params.t[2][1] = 5.f;  // E
    params.t[2][2] = 10.f;  // K
    params.t[2][3] = 1.f;  // terminal

    // transitions from the terminal state
    params.t[3][0] = 0.f; // M
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
    std::vector<HMMConsReadState> read_states;

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
void add_read(CSquiggleReadInterface params)
{
    g_data.reads.push_back(CSquiggleRead());
    CSquiggleRead& sr = g_data.reads.back();
 
    for(uint32_t i = 0; i < NUM_STRANDS; ++i) {
        // Initialize pore model   
        sr.pore_model[i].scale = params.pore_model[i].scale;
        sr.pore_model[i].shift = params.pore_model[i].shift;
        
        assert(params.pore_model[i].n_states == 1024);
        for(uint32_t j = 0; j < params.pore_model[i].n_states; ++j) {
            sr.pore_model[i].state[j].mean = params.pore_model[i].mean[j];
            sr.pore_model[i].state[j].sd = params.pore_model[i].sd[j];
         }
    
        // Initialize events
        sr.events[i].n_events = params.events[i].n_events;
        sr.events[i].level = params.events[i].level;
        
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
    rs.event_idx = params.event_idx;
    rs.kmer_idx = 0;
    rs.strand = params.strand;
    rs.stride = params.stride;
    rs.rc = params.rc;
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

void free_matrix(HMMMatrix matrix)
{
    free(matrix.cells);
    matrix.cells = NULL;
}

inline uint32_t cell(const HMMMatrix& matrix, uint32_t row, uint32_t col)
{
    return row * matrix.n_cols + col;
}

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
    } while(carry > 0);
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

inline double log_probability_match(const CSquiggleRead& read,
                                    uint32_t kmer_rank,
                                    uint32_t event_idx, 
                                    uint8_t strand)
{
    const CPoreModel& pm = read.pore_model[strand];

    // Extract event
    double level = read.events[strand].level[event_idx];
    
    double m = (pm.state[kmer_rank].mean + pm.shift) * pm.scale;
    double s = pm.state[kmer_rank].sd * pm.scale;
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

void initialize_forward(HMMMatrix& matrix)
{
    //
    uint32_t c = cell(matrix, 0, 0);
    matrix.cells[c].M = log(1.0);
    matrix.cells[c].E = -INFINITY;
    matrix.cells[c].K = -INFINITY;

    // Initialize first row/column to prevent initial gaps
    for(uint32_t i = 1; i < matrix.n_rows; i++) {
        uint32_t c = cell(matrix, i, 0);
        matrix.cells[c].M = -INFINITY;
        matrix.cells[c].E = -INFINITY;
        matrix.cells[c].K = -INFINITY;
    }

    for(uint32_t j = 1; j < matrix.n_cols; j++) {
        uint32_t c = cell(matrix, 0, j);
        matrix.cells[c].M = -INFINITY;
        matrix.cells[c].E = -INFINITY;
        matrix.cells[c].K = -INFINITY;
    }
}

void initialize_backward(HMMMatrix& matrix, const HMMParameters& hmm_params)
{
    //
    uint32_t c = cell(matrix, matrix.n_rows - 1, matrix.n_cols - 1);

    // the bottom right corner of the matrix is initialized to
    // the probability of transitioning to the terminal state
    matrix.cells[c].M = hmm_params.t[0][3];
    matrix.cells[c].E = hmm_params.t[1][3];
    matrix.cells[c].K = hmm_params.t[2][3];
}

double fill_forward(HMMMatrix& matrix, 
                    const HMMParameters& hmm_params, 
                    const char* sequence,
                    const HMMConsReadState& state,
                    uint32_t e_start, 
                    uint32_t k_start)
{
    // Fill in matrix
    for(uint32_t row = 1; row < matrix.n_rows; row++) {
        for(uint32_t col = 1; col < matrix.n_cols; col++) {

            // cell indices
            uint32_t c = cell(matrix, row, col);
            uint32_t diag = cell(matrix, row - 1, col - 1);
            uint32_t up =   cell(matrix, row - 1, col);
            uint32_t left = cell(matrix, row, col - 1);

            uint32_t event_idx = e_start + (row - 1) * state.stride;
            uint32_t kmer_idx = k_start + col - 1;

            // Emission probability for a match
            uint32_t rank = !state.rc ? 
                kmer_rank(sequence + kmer_idx, K) : 
                rc_kmer_rank(sequence + kmer_idx, K);
            double l_p_m = log_probability_match(*state.read, rank, event_idx, state.strand);

            // Emission probility for an event insertion
            double l_p_e = log_probability_event_insert(*state.read, rank, event_idx, state.strand);
            
            // Emission probability for a kmer insertion
            uint32_t next_rank = !state.rc ? 
                kmer_rank(sequence + kmer_idx + 1, K) : 
                rc_kmer_rank(sequence + kmer_idx + 1, K);
            double l_p_k = log_probability_kmer_insert(*state.read, rank, event_idx, state.strand);

            // Calculate M[i, j]
            double d_m = hmm_params.t[0][0] + matrix.cells[diag].M;
            double d_e = hmm_params.t[1][0] + matrix.cells[diag].E;
            double d_k = hmm_params.t[2][0] + matrix.cells[diag].K;
            matrix.cells[c].M = l_p_m + log(exp(d_m) + exp(d_e) + exp(d_k));

            // Calculate E[i, j]
            double u_m = hmm_params.t[0][1] + matrix.cells[up].M;
            double u_e = hmm_params.t[1][1] + matrix.cells[up].E;
            double u_k = hmm_params.t[2][1] + matrix.cells[up].K;
            matrix.cells[c].E = l_p_e + log(exp(u_m) + exp(u_e) + exp(u_k));

            // Calculate K[i, j]
            double l_m = hmm_params.t[0][2] + matrix.cells[left].M;
            double l_e = hmm_params.t[1][2] + matrix.cells[left].E;
            double l_k = hmm_params.t[2][2] + matrix.cells[left].K;
            matrix.cells[c].K = l_p_k + log(exp(l_m) + exp(l_e) + exp(l_k));

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) R -- [%.2lf %.2lf %.2lf]\n", row, col, matrix.cells[c].M, matrix.cells[c].E, matrix.cells[c].K);
            printf("(%d %d) D -- e: %.2lf [%.2lf %.2lf %.2lf]\n", row, col, l_p_m, d_m, d_e, d_k);
            printf("(%d %d) U -- e: %.2lf [%.2lf %.2lf %.2lf]\n", row, col, l_p_e, u_m, u_e, u_k);
            printf("(%d %d) L -- e: %.2lf [%.2lf %.2lf %.2lf]\n", row, col, l_p_k, l_m, l_e, l_k);
#endif
        }
    }

    uint32_t c = cell(matrix, matrix.n_rows - 1, matrix.n_cols - 1);
    return log( exp(matrix.cells[c].M + hmm_params.t[0][3]) +
                exp(matrix.cells[c].E + hmm_params.t[1][3]) + 
                exp(matrix.cells[c].K + hmm_params.t[2][3]) );
}

void fill_backward(HMMMatrix& matrix, 
                   const HMMParameters& hmm_params, 
                   const char* sequence,
                   const HMMConsReadState& state,
                   uint32_t e_start, 
                   uint32_t k_start)
{
    uint32_t nr = matrix.n_rows;
    uint32_t nc = matrix.n_cols;

    // Fill in matrix
    for(uint32_t row = nr - 1; row > 0; row--) {
        for(uint32_t col = nc - 1; col > 0; col--) {

            // skip bottom right corner
            if(row == matrix.n_rows - 1 && col == matrix.n_cols - 1)
                continue;

            // cell indices
            uint32_t c = cell(matrix, row, col);
            uint32_t diag = cell(matrix, row + 1, col + 1);
            uint32_t down = cell(matrix, row + 1, col);
            uint32_t right = cell(matrix, row, col + 1);
            
            double v_m = -INFINITY;
            if(row < nr - 1 && col < nc - 1) { 
                // Emission probability for matching e_(i+1) to k_(j+1)
                // this is for row + 1 and col + 1, respectively
                uint32_t event_idx = e_start + row * state.stride;
                uint32_t kmer_idx = k_start + col;

                // Emission probability for a match
                uint32_t rank = !state.rc ? 
                    kmer_rank(sequence + kmer_idx, K) : 
                    rc_kmer_rank(sequence + kmer_idx, K);
                double l_p_m = log_probability_match(*state.read, rank, event_idx, state.strand);
                v_m = l_p_m + matrix.cells[diag].M;
            }

            double v_e = -INFINITY;
            if(row < nr - 1) {
                // Emission probability for skipping event e_(i+1)
                // this is for row + 1 and col, respectively
                uint32_t event_idx = e_start + row * state.stride;
                uint32_t kmer_idx = k_start + col - 1;

                // Emission probability for a match
                uint32_t rank = !state.rc ? 
                    kmer_rank(sequence + kmer_idx, K) : 
                    rc_kmer_rank(sequence + kmer_idx, K);
                double l_p_e = log_probability_event_insert(*state.read, rank, event_idx, state.strand);
                v_e = l_p_e + matrix.cells[down].E;
            }

            // Emission probability for skipping kmer k_(j+1)
            double v_k = -INFINITY;
            if(col < nc - 1) {
                uint32_t event_idx = e_start + (row - 1) * state.stride;
                uint32_t kmer_idx = k_start + col;
                uint32_t rank = !state.rc ? 
                    kmer_rank(sequence + kmer_idx, K) : 
                    rc_kmer_rank(sequence + kmer_idx, K);
                
                double l_p_k = log_probability_kmer_insert(*state.read, rank, event_idx, state.strand);
                v_k = l_p_k + matrix.cells[right].K;
            }

            // Calculate M[i, j], E[i, j], K[i, j]
            matrix.cells[c].M = log( exp(hmm_params.t[0][0] + v_m) +
                                     exp(hmm_params.t[0][1] + v_e) +
                                     exp(hmm_params.t[0][2] + v_k) );

            matrix.cells[c].E = log( exp(hmm_params.t[1][0] + v_m) +
                                     exp(hmm_params.t[1][1] + v_e) +
                                     exp(hmm_params.t[1][2] + v_k) );

            matrix.cells[c].K = log( exp(hmm_params.t[2][0] + v_m) +
                                     exp(hmm_params.t[2][1] + v_e) +
                                     exp(hmm_params.t[2][2] + v_k) );
        }
    }

}

ExtensionResult run_extension_hmm(const std::string& consensus, const HMMConsReadState& state)
{
    double time_start = clock();

    std::string root(consensus.c_str() + state.kmer_idx);
    printf("ROOT: %s e: %d k: %d\n", root.c_str(), state.event_idx, state.kmer_idx);
    std::string extension = root + "AAAAA";

    // Get the start/end event indices
    uint32_t e_start = state.event_idx;
    uint32_t e_end = e_start + extension.size() + 10;
    uint32_t n_events = e_end - e_start + 1;
    uint32_t k_start = 0; // this is in reference to the extension sequence
    uint32_t n_kmers = extension.size() - K + 1;
 
    // Set up HMM matrix
    HMMMatrix matrix;
    allocate_matrix(matrix, n_events + 1, n_kmers + 1);

    std::string best_path_str;
    double best_path_score = -INFINITY;
       
    ExtensionResult result;
    for(uint8_t i = 0; i < 4; ++i)
        result.b[i] = -INFINITY;

    uint32_t extension_rank = 0;
    while(extension.substr(0, extension.size() - K) == root) {
        
        initialize_forward(matrix);

        // Fill in the HMM matrix using the forward algorithm
        double l_f = fill_forward(matrix, g_data.hmm_params, extension.c_str(), state, e_start, k_start);

        // Determine the best scoring row in the last column
        uint32_t col = matrix.n_cols - 1;
        uint32_t max_row = 0;
        double max_value = -INFINITY;
        
        for(uint32_t row = 3; row < matrix.n_rows; ++row) {
            uint32_t c = cell(matrix, row, col);
            double sum = log(exp(matrix.cells[c].M) + exp(matrix.cells[c].E) + exp(matrix.cells[c].K));
            if(sum > max_value) {
                max_value = sum;
                max_row = row;
            }
        }

        printf("extensions: %s %d %.2lff\n", 
                extension.substr(extension.size() - K).c_str(), max_row, max_value);
        
        // path sum
        uint8_t br = base_rank[extension[extension.size() - K]];
        double kmer_sum = log(exp(result.b[br]) + exp(max_value));
        result.b[br] = kmer_sum;

        if(max_value > best_path_score) {
            best_path_score = max_value;
            best_path_str = extension.substr(extension.size() - K);
        }

        // Set the extension to the next string
        lexicographic_next(extension);
    }

    double time_stop = clock();
    //printf("Time: %.2lfs\n", (time_stop - time_start) / CLOCKS_PER_SEC);
    free_matrix(matrix);
    result.best_path = best_path_str;
    result.best_path_score = best_path_score;
    return result;
}

void debug_align(const std::string& consensus, const HMMConsReadState& state)
{
    printf("\nAligning read to %s\n", consensus.c_str());
    double time_start = clock();

    // Get the start/end event indices
    uint32_t e_start = state.event_idx;
    uint32_t e_end = e_start + consensus.size() + 10;
    uint32_t n_events = e_end - e_start + 1;
    uint32_t k_start = 0; // this is in reference to the extension sequence
    uint32_t n_kmers = consensus.size() - K + 1;
 
    // Set up HMM matrix
    HMMMatrix matrix;
    allocate_matrix(matrix, n_events + 1, n_kmers + 1);
    
    initialize_forward(matrix);

    // Fill in the HMM matrix using the forward algorithm
    double l_f = fill_forward(matrix, g_data.hmm_params, consensus.c_str(), state, e_start, k_start);

    // Determine the best scoring row in the last column
    uint32_t col = matrix.n_cols - 1;
    uint32_t max_row = 0;
    double max_value = -INFINITY;

    for(uint32_t row = 3; row < matrix.n_rows; ++row) {
        uint32_t c = cell(matrix, row, col);
        double sum = log(exp(matrix.cells[c].M) + exp(matrix.cells[c].E) + exp(matrix.cells[c].K));
        if(sum > max_value) {
            max_value = sum;
            max_row = row;
        }
    }

    HMMMatrix b_matrix;
    allocate_matrix(b_matrix, n_events + 1, n_kmers + 1);
    initialize_backward(b_matrix, g_data.hmm_params);
    fill_backward(b_matrix, g_data.hmm_params, consensus.c_str(), state, e_start, k_start);
    
    print_matrix(matrix);
    print_matrix(b_matrix);

    printf("l_f: %.2f posterior:\n", l_f);
    for(uint32_t col = 1; col < matrix.n_cols; ++col)
        printf("%s\t", consensus.substr(k_start + col - 1, K).c_str());
    printf("\n");
    for(uint32_t row = 1; row < matrix.n_rows; ++row) {
        for(uint32_t col = 1; col < matrix.n_cols; ++col) {
            uint32_t c = cell(matrix, row, col);

            double p_m = exp(matrix.cells[c].M + b_matrix.cells[c].M - l_f);
            double p_e = exp(matrix.cells[c].E + b_matrix.cells[c].E - l_f);
            double p_k = exp(matrix.cells[c].K + b_matrix.cells[c].K - l_f);
            printf("%.2lf,%.2lf,%.2lf\t", p_m, p_e, p_k);
        }
        printf("\n");
    }

    // posterior decode
    std::string states;
    std::vector<double> posteriors;
    uint32_t row = max_row;
    col = matrix.n_cols - 1;
    while(row > 0 && col > 0) {
        uint32_t c = cell(matrix, row, col);
        double p_m = exp(matrix.cells[c].M + b_matrix.cells[c].M - l_f);
        double p_e = exp(matrix.cells[c].E + b_matrix.cells[c].E - l_f);
        double p_k = exp(matrix.cells[c].K + b_matrix.cells[c].K - l_f);
        
        if(p_m > p_e && p_m > p_k) {
            states.append(1, 'M');
            posteriors.push_back(p_m);
            col -= 1;
            row -= 1;
        } else if (p_e > p_k) {
            states.append(1, 'E');
            posteriors.push_back(p_e);
            row -= 1;
        } else {
            states.append(1, 'K');
            posteriors.push_back(p_k);
            col -= 1;
        }
    }

    std::reverse(states.begin(), states.end());
    std::reverse(posteriors.begin(), posteriors.end());
    
    printf("Align max value: %.2lf\n", max_value);
    uint32_t ei = e_start;
    uint32_t ki = k_start;
    for(size_t i = 0; i < states.size(); ++i) {
        char s = states[i];
        
        double level = state.read->events[state.strand].level[ei];
        uint32_t rank = !state.rc ? 
            kmer_rank(consensus.c_str() + ki, K) : 
            rc_kmer_rank(consensus.c_str() + ki, K);
        
        const CPoreModel& pm = state.read->pore_model[state.strand];
        std::string kmer(consensus.c_str(), K);
        //printf("%s %d %lf %lf %lf\n", kmer.c_str(), rank, pm.state[rank].mean, pm.shift, pm.scale);
        
        double model_m = (pm.state[rank].mean + pm.shift) * pm.scale;
        double model_s = pm.state[rank].sd * pm.scale;
        printf("%c %d %d %.2lf 0.00 %s %.2lf %.2lf %.2lf\n", s, s != 'K' ? ei : -1, ki, level, consensus.substr(ki, K).c_str(), model_m, model_s, posteriors[i]);
    
        if(states[i + 1] == 'M') {
            ei += state.stride;
            ki += 1;
        } else if(states[i + 1] == 'E') {
            ei += state.stride;
        } else {
            ki += 1;
        }
    }

    free(b_matrix.cells);
    free_matrix(matrix);
}

void update_read_state(const std::string& consensus, HMMConsReadState& state)
{
    // Calculate the forward and backward matrices for this read and consensus sequence
    
    // Get the start/end event indices
    uint32_t e_start = state.event_idx;
    uint32_t e_end = e_start + 10;
    uint32_t n_events = e_end - e_start + 1;

    uint32_t k_start = state.kmer_idx;
    uint32_t n_kmers = (consensus.size() - K + 1) - k_start;
 
    // Forward algorithm
    HMMMatrix f_matrix;
    allocate_matrix(f_matrix, n_events + 1, n_kmers + 1);
    initialize_forward(f_matrix);
    double l_f = fill_forward(f_matrix, g_data.hmm_params, consensus.c_str(), state, e_start, k_start);

    // Backward algorithm
    HMMMatrix b_matrix;
    allocate_matrix(b_matrix, n_events + 1, n_kmers + 1);
    initialize_backward(b_matrix, g_data.hmm_params);
    fill_backward(b_matrix, g_data.hmm_params, consensus.c_str(), state, e_start, k_start);

    /*
    printf("forward (%.2lf):\n", l_f);
    print_matrix(f_matrix);
    printf("backward:\n");
    print_matrix(b_matrix);
    
    printf("posterior:\n");
    for(uint32_t row = 1; row < f_matrix.n_rows; ++row) {
        for(uint32_t col = 1; col < f_matrix.n_cols; ++col) {
            uint32_t c = cell(f_matrix, row, col);

            double p_m = exp(f_matrix.cells[c].M + b_matrix.cells[c].M - l_f);
            double p_e = exp(f_matrix.cells[c].E + b_matrix.cells[c].E - l_f);
            double p_k = exp(f_matrix.cells[c].K + b_matrix.cells[c].K - l_f);
            printf("%.2lf, %.2lf, %.2lf\t", p_m, p_e, p_k);
        }
        printf("\n");
    }
    */

    // Update event and kmer index if there is an event with high posterior of being matched to the new kmer
    uint32_t col = 2;
    double max_posterior = 0.0f;
    uint32_t max_row = 0;
    for(uint32_t row = 2; row < f_matrix.n_rows; ++row) {
        uint32_t c = cell(f_matrix, row, col);
        double p_m = exp(f_matrix.cells[c].M + b_matrix.cells[c].M - l_f);
        if(p_m > max_posterior) {
            max_posterior = p_m;
            max_row = row;
        }
    }

    printf("Max posterior: %.2lf in row %d\n", max_posterior, max_row);

    if(max_posterior > 0.25f) {
        state.kmer_idx += 1;
        state.event_idx += max_row - 1;
    }

    free(f_matrix.cells);
    free(b_matrix.cells);
}

extern "C"
void run_debug()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    std::string consensus = "AACAGTCCACTATTGGATG";
    for(uint32_t i = 0; i < g_data.read_states.size(); ++i) {
        debug_align(consensus, g_data.read_states[i]);
        debug_align("AACAGTCCACTATTAAGT", g_data.read_states[i]);
    }
}

extern "C"
void run_consensus()
{
    if(!g_initialized) {
        printf("ERROR: initialize() not called\n");
        exit(EXIT_FAILURE);
    }

    std::string consensus = "AACAGTCCACTATT";
    
    int iterations = 2;
    while(--iterations > 0) {

        ExtensionResult all = { 0, 0, 0, 0 };
        for(uint32_t i = 0; i < g_data.read_states.size(); ++i) {
            ExtensionResult r = run_extension_hmm(consensus, g_data.read_states[i]);

            // Normalize by the sum over all bases for this sequence
            double sequence_sum = -INFINITY;
            for(uint32_t j = 0; j < 4; ++j)
                sequence_sum = log(exp(sequence_sum) + exp(r.b[j]));
            for(uint32_t j = 0; j < 4; ++j) {
                r.b[j] -= sequence_sum;
                all.b[j] += r.b[j];
            }
            printf("seq[%d]\tLP(A): %.2lf LP(C): %.2lf LP(G): %.2lf LP(T): %.2lf path: %s %.2lf\n", i, r.b[0], r.b[1], r.b[2], r.b[3], r.best_path.c_str(), r.best_path_score);
        
            printf("Debug alignment to best\n");
            debug_align(consensus + r.best_path, g_data.read_states[i]);
            
            printf("Debug alignment to truth\n");
            debug_align("AACAGTCCACTATTGGATG", g_data.read_states[i]);
        }
        printf("seq[a]\tLP(A): %.2lf LP(C): %.2lf LP(G): %.2lf LP(T): %.2lf\n", all.b[0], all.b[1], all.b[2], all.b[3]);

        double best_lp = -INFINITY;
        char best_base;
        for(uint8_t i = 0; i < 4; ++i) {
            if(all.b[i] > best_lp) {
                best_lp = all.b[i];
                best_base = "ACGT"[i];
            }
        }
        printf("Extending to %c\n", best_base);
        consensus.append(1, best_base);
    }
    printf("Consensus: %s\n", consensus.c_str());
    //for(uint32_t i = 0; i < g_data.read_states.size(); ++i) {
    //    update_read_state(consensus, g_data.read_states[i]);
    //}
}

int main(int argc, char** argv)
{

}
