//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_klcs -- function to compute the longest
// common subsequence of k-mers for two strings
//
#include <algorithm>
#include "nanopolish_klcs.h"

// Helper to backtrack through the kLCS matrix
void _kLCSBacktrack(const UInt32Matrix& m,
                    const std::string& a, 
                    const std::string& b,
                    const int k,
                    uint32_t row,
                    uint32_t col,
                    kLCSResult& result)
{
    if(row == 0 || col == 0)
        return;

    const char* ka = a.c_str() + row - 1;
    const char* kb = b.c_str() + col - 1;

    if(strncmp(ka, kb, k) == 0) {
        kLCSPair p = { row - 1, col - 1 };
        result.push_back(p);
        return _kLCSBacktrack(m, a, b, k, row - 1, col - 1, result);
    } else {

        if(get(m, row - 1, col) > get(m, row, col - 1)) {
            return _kLCSBacktrack(m, a, b, k, row - 1, col, result);
        } else {
            return _kLCSBacktrack(m, a, b, k, row, col - 1, result);
        }
    }
}

// Return the longest common subseuqence of k-mers between the two strings
kLCSResult kLCS(const std::string& a, const std::string& b, const int k)
{
    uint32_t n_kmers_a = a.size() - k + 1;
    uint32_t n_kmers_b = b.size() - k + 1;

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
            if(strncmp(ka, kb, k) == 0) {
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
    _kLCSBacktrack(m, a, b, k, n_rows - 1, n_cols -  1, result);

    // Backtrack appends from the end to the start, reverse the vector of matches
    std::reverse(result.begin(), result.end());
    free_matrix(m);
    return result;
}
