#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
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
    uint32_t read_idx;
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
};

// A global vector used to store data we've received from the python code
struct HmmConsData
{
    std::vector<CSquiggleRead> reads;
    std::vector<HMMConsReadState> read_states;
};
HmmConsData g_data;

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
    rs.read_idx = params.read_idx;
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

void print_matrix(const HMMMatrix& matrix, uint32_t n_rows, uint32_t n_cols)
{
    for(uint32_t i = 0; i < n_rows; ++i) {
        for(uint32_t j = 0; j < n_cols; ++j) {
            uint32_t c = cell(matrix, i, j);
            printf("%.2lf\t", matrix.cells[c].M);
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

void fill_forward(HMMMatrix& matrix, 
                  double t[][3], 
                  const char* sequence,
                  const CSquiggleRead& read, 
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
            uint32_t kmer_idx = col - 1;

            // Emission probability for a match
            uint32_t rank = !state.rc ? 
                kmer_rank(sequence + kmer_idx, K) : 
                rc_kmer_rank(sequence + kmer_idx, K);
            double l_p_m = log_probability_match(read, rank, event_idx, state.strand);

            // Emission probility for an event insertion
            // This is calculated using the emission probability for a match to the same kmer as the previous row
            double l_p_e = l_p_m;

            // Emission probability for a kmer insertion
            double l_p_k = LOG_KMER_INSERTION;

            // Calculate M[i, j]
            double d_m = t[0][0] + matrix.cells[diag].M;
            double d_e = t[1][0] + matrix.cells[diag].E;
            double d_k = t[2][0] + matrix.cells[diag].K;
            matrix.cells[c].M = l_p_m + log(exp(d_m) + exp(d_e) + exp(d_k));

            // Calculate E[i, j]
            double u_m = t[0][1] + matrix.cells[up].M;
            double u_e = t[1][1] + matrix.cells[up].E;
            double u_k = t[2][1] + matrix.cells[up].K;
            matrix.cells[c].E = l_p_e + log(exp(u_m) + exp(u_e) + exp(u_k));

            // Calculate K[i, j]
            double l_m = t[0][2] + matrix.cells[left].M;
            double l_e = t[1][2] + matrix.cells[left].E;
            double l_k = t[2][2] + matrix.cells[left].K;
            matrix.cells[c].K = l_p_k + log(exp(l_m) + exp(l_e) + exp(l_k));

#ifdef DEBUG_HMM_UPDATE
            printf("(%d %d) R -- [%.2lf %.2lf %.2lf]\n", row, col, matrix.cells[c].M, matrix.cells[c].E, matrix.cells[c].K);
            printf("(%d %d) D -- e: %.2lf t: [%.2lf %.2lf %.2lf] [%.2lf %.2lf %.2lf]\n", row, col, l_p_m, t[0][0], t[1][0], t[2][0], d_m, d_e, d_k);
            printf("(%d %d) U -- e: %.2lf t: [%.2lf %.2lf %.2lf] [%.2lf %.2lf %.2lf]\n", row, col, l_p_e, t[0][1], t[1][1], t[2][1], u_m, u_e, u_k);
            printf("(%d %d) L -- e: %.2lf t: [%.2lf %.2lf %.2lf] [%.2lf %.2lf %.2lf]\n", row, col, l_p_k, t[0][2], t[1][2], t[2][2], l_m, l_e, l_k);
#endif
        }
    }

}


ExtensionResult run_extension_hmm(const std::string& consensus, const HMMConsReadState& state)
{
    double time_start = clock();
    const CSquiggleRead& read = g_data.reads[state.read_idx];

    // TODO: not constant
    double LOG_KMER_INSERTION = log(0.1);

    // The root of the extension sequences to test is the last k-mer
    std::string root_kmer = consensus.substr(consensus.size() - K);
    std::string extension = root_kmer + "AAAAA";
    //std::string extension = root_kmer + "TCCAC";

    // Get the start/end event indices
    uint32_t e_start = state.event_idx;
    uint32_t e_end = e_start + 10;

    uint32_t k_start = 0;
    uint32_t n_kmers = extension.size() - K + 1;
 
    // Setup transition matrix
    static const uint32_t n_states = 3;
    double t[n_states][n_states] = { { 0.90f, 0.05f, 0.05f },   // MM, ME, MK
                                     { 0.85f, 0.10f, 0.05f },   // EM, EE, EK
                                     { 0.85f, 0.05f, 0.10f } }; // EM, EE, EK

    // Log scale the transition matrix
    for(uint32_t i = 0; i < n_states; ++i) {
        for(uint32_t j = 0; j < n_states; ++j) {
            t[i][j] = log(t[i][j]);
        }
    }
    
    // Set up HMM matrix
    HMMMatrix matrix;
    allocate_matrix(matrix, e_end - e_start + 2, n_kmers + 1);
    
    ExtensionResult result;
    for(uint8_t i = 0; i < 4; ++i)
        result.b[i] = -INFINITY;

    uint32_t extension_rank = 0;
    while(extension.substr(0, K) == root_kmer) {
        
        initialize_forward(matrix);

        // Fill in the HMM matrix using the forward algorithm
        fill_forward(matrix, t, extension.c_str(), read, state, e_start, k_start);

        // Determine the best scoring row in the last column
        uint32_t col = matrix.n_cols - 1;
        uint32_t max_row = 0;
        double max_value = -INFINITY;
        
        for(uint32_t row = 2; row < matrix.n_rows; ++row) {
            uint32_t c = cell(matrix, row, col);
            double sum = log(exp(matrix.cells[c].M) + exp(matrix.cells[c].E) + exp(matrix.cells[c].K));
            if(sum > max_value) {
                max_value = sum;
                max_row = row;
            }
        }

        /*
        printf("extensions: %s %d %.2lff\n", 
                extension.substr(extension.size() - K).c_str(), max_row, max_value);
        */
        // path sum
        uint8_t br = base_rank[extension[extension.size() - K]];
        double kmer_sum = log(exp(result.b[br]) + exp(max_value));
        result.b[br] = kmer_sum;
        
        // Set the extension to the next string
        lexicographic_next(extension);
    }

    double time_stop = clock();
    //printf("Time: %.2lfs\n", (time_stop - time_start) / CLOCKS_PER_SEC);
    free_matrix(matrix);

    return result;
}

extern "C"
void run_consensus()
{
    std::string consensus = "AACAG";

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
        printf("seq[%d]\tLP(A): %.2lf LP(C): %.2lf LP(G): %.2lf LP(T): %.2lf\n", i, r.b[0], r.b[1], r.b[2], r.b[3]);
    }
    printf("seq[a]\tLP(A): %.2lf LP(C): %.2lf LP(G): %.2lf LP(T): %.2lf\n", all.b[0], all.b[1], all.b[2], all.b[3]);
}

int main(int argc, char** argv)
{

}
