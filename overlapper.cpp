//-------------------------------------------------------------------------------
// 
// overlapper - String-string overlap algorithm 
//
// Copyright (C) 2011 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "overlapper.h"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <limits>
#include <stdio.h>
#include <inttypes.h>

// 
OverlapperParams default_params = { 2, -6, -3, -2, ALT_OVERLAP };
OverlapperParams ungapped_params = { 2, -10000, -3, -2, ALT_OVERLAP };
OverlapperParams affine_default_params = { 2, -5, -3, -2, ALT_OVERLAP, true, true, true, true, true };

//
#define max3(x,y,z) std::max(std::max(x,y), z)
//#define DEBUG_OVERLAPPER 1
//#define DEBUG_EXTEND 1

// 
SequenceInterval::SequenceInterval() : start(0), end(-1)
{

}

SequenceOverlap::SequenceOverlap()
{
    length[0] = length[1] = 0;
    score = -1;
    edit_distance = -1;
    total_columns = -1;
}

//
bool SequenceOverlap::isValid() const
{
    return !cigar.empty() && match[0].isValid() && match[1].isValid();
}

//
double SequenceOverlap::getPercentIdentity() const
{
    return (double)(total_columns - edit_distance) * 100.0f / total_columns;
}

//
std::ostream& operator<<(std::ostream& out, const SequenceOverlap& overlap)
{
    out << "[" << overlap.match[0].start << " " << overlap.match[0].end << "] ";
    out << "[" << overlap.match[1].start << " " << overlap.match[1].end << "] ";
    out << "C:" << overlap.cigar;
    return out;
}

void SequenceOverlap::makePaddedMatches(const std::string& s1, const std::string& s2,
                                        std::string* p1, std::string* p2) const
{
    assert(isValid() && p1 != NULL && p2 != NULL);

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            p1->append(s1.substr(current_1, length));
            p2->append(s2.substr(current_2, length));
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            p1->append(s1.substr(current_1, length));
            p2->append(length, '-');
            current_1 += length;
        }
        else if(code == 'I') {
            p1->append(length, '-');
            p2->append(s2.substr(current_2, length));
            current_2 += length;
        }
        length = -1;
    }
}

//
int SequenceOverlap::calculateEditDistance(const std::string& s1, const std::string& s2) const
{
    // Recalculate the edit distance between the pair of strings, given this alignment
    int new_edit_distance = 0;

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M') {
            for(int i = 0; i < length; ++i) {
                if(s1[current_1 + i] != s2[current_2 + i])
                    new_edit_distance++;
            }
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            new_edit_distance += length;
            current_1 += length;
        }
        else if(code == 'I') {
            new_edit_distance += length;
            current_2 += length;
        }
        length = -1;
    }

    return new_edit_distance;
}

//
int SequenceOverlap::calculateTotalColumns() const
{
    // Recalculate the edit distance between the pair of strings, given this alignment
    int total_columns = 0;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        total_columns += length;
    }

    return total_columns;
}

//
void SequenceOverlap::printAlignment(const std::string& s1, const std::string& s2) const
{
    if(!isValid()) {
        std::cerr << "Alignment is not valid\n";
        return;
    }

    std::string out_1;
    std::string out_2;

    // Print out the initial part of the strings, which do not match. 
    // Typically this is the overhanging portion of one of the strings.
    std::string leader_1 = s1.substr(0, match[0].start);
    std::string leader_2 = s2.substr(0, match[1].start);

    // Pad the beginning of the output strings with spaces to align
    if(leader_1.size() < leader_2.size())
        out_1.append(leader_2.size() - leader_1.size(), ' ');

    if(leader_2.size() < leader_1.size())
        out_2.append(leader_1.size() - leader_2.size(), ' ');
    
    out_1.append(leader_1);
    out_2.append(leader_2);

    // Process the matching region using the cigar operations
    size_t current_1 = match[0].start;
    size_t current_2 = match[1].start;

    std::stringstream cigar_parser(cigar);
    int length = -1;
    char code;
    while(cigar_parser >> length >> code) {
        assert(length > 0);
        if(code == 'M' || code == 'X' || code == '=') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(s2.substr(current_2, length));
            current_1 += length;
            current_2 += length;
        }
        else if(code == 'D') {
            out_1.append(s1.substr(current_1, length));
            out_2.append(length, '-');
            current_1 += length;
        }
        else if(code == 'I') {
            out_1.append(length, '-');
            out_2.append(s2.substr(current_2, length));
            current_2 += length;
        }
        length = -1;
    }

    // Append the remainder of each string
    out_1.append(s1.substr(current_1));
    out_2.append(s2.substr(current_2));

    // Print the output strings and split long lines
    int MAX_COLUMNS = 120;
    size_t total_columns = std::max(out_1.size(), out_2.size());
    for(size_t i = 0; i < total_columns; i += MAX_COLUMNS) {
        std::string sub_1;
        std::string sub_2;
        if(i < out_1.size())
            sub_1 = out_1.substr(i, MAX_COLUMNS);
        if(i < out_2.size())
            sub_2 = out_2.substr(i, MAX_COLUMNS);
        
        std::cout << "S1\t" << sub_1 << "\n";
        std::cout << "S2\t" << sub_2 << "\n";
        std::cout << "\n";
    }
    std::cout << "Cigar: " << cigar << "\n";
    std::cout << "Score: " << score << "\n";

    printf("Identity: %2.2lf\n", getPercentIdentity());
}

typedef std::vector<int> DPCells;
typedef std::vector<DPCells> DPMatrix;

//
SequenceOverlap Overlapper::computeOverlap(const std::string& s1, const std::string& s2, const OverlapperParams params)
{
    // Exit with invalid intervals if either string is zero length
    SequenceOverlap output;
    if(s1.empty() || s2.empty()) {
        std::cerr << "Overlapper::computeOverlap error: empty input sequence\n";
        exit(EXIT_FAILURE);
    }

    // Initialize the scoring matrix
    size_t num_columns = s1.size() + 1;
    size_t num_rows = s2.size() + 1;

    DPMatrix score_matrix;
    score_matrix.resize(num_columns);
    for(size_t i = 0; i < score_matrix.size(); ++i)
        score_matrix[i].resize(num_rows);

    // Calculate scores
    for(size_t i = 1; i < num_columns; ++i) {
        for(size_t j = 1; j < num_rows; ++j) {
            // Calculate the score for entry (i,j)
            int idx_1 = i - 1;
            int idx_2 = j - 1;
            int diagonal = score_matrix[i-1][j-1] + (s1[idx_1] == s2[idx_2] ? params.match_score : params.mismatch_penalty);
            int up = score_matrix[i][j-1] + params.gap_penalty;
            int left = score_matrix[i-1][j] + params.gap_penalty;

            score_matrix[i][j] = max3(diagonal, up, left);
        }
    }
 
    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(size_t i = 1; i < num_columns; ++i) {
        int v = score_matrix[i][num_rows - 1];
        if(score_matrix[i][num_rows - 1] > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(size_t j = 1; j < num_rows; ++j) {
        int v = score_matrix[num_columns - 1][j];
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }

    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_OVERLAPPER
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = score_matrix[i - 1][j - 1] + (is_match ? params.match_score : params.mismatch_penalty);
        int up = score_matrix[i][j-1] + params.gap_penalty;
        int left = score_matrix[i-1][j] + params.gap_penalty;

        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(score_matrix[i][j] == up) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(score_matrix[i][j] == left) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
            assert(score_matrix[i][j] == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}

// Returns the index into a cell vector for for the ith column and jth row
// of a dynamic programming matrix. The band_origin gives the row in first
// column of the matrix that the bands start at. This is used to calculate
// the starting band row for each column.
inline int _getBandedCellIndex(int i, int j, int band_width, int band_origin_row)
{
    int band_start = band_origin_row + i;
    int band_row_index = j - band_start;
    return (band_row_index >= 0 && band_row_index < band_width) ? i * band_width + band_row_index : -1;
}

// Returns the score for (i,j) in the 
inline int _getBandedCellScore(const DPCells& cells, int i, int j, int band_width, int band_origin_row, int invalid_score)
{
    int band_start = band_origin_row + i;
    int band_row_index = j - band_start;
    return (band_row_index >= 0 && band_row_index < band_width) ? cells[i * band_width + band_row_index] : invalid_score;
}

SequenceOverlap Overlapper::extendMatch(const std::string& s1, const std::string& s2, 
                                        int start_1, int start_2, int band_width)
{
    SequenceOverlap output;
    int num_columns = s1.size() + 1;
    int num_rows = s2.size() + 1;

    const int MATCH_SCORE = 2;
    const int GAP_PENALTY = -5;
    const int MISMATCH_PENALTY = -3;
    
    // Calculate the number of cells off the diagonal to compute
    int half_width = band_width / 2;
    band_width = half_width * 2 + 1; // the total number of cells per band

    // Calculate the number of columns that we need to extend to for s1
    size_t num_cells_required = num_columns * band_width;

    // Allocate bands with uninitialized scores
    int INVALID_SCORE = std::numeric_limits<int>::min();
    DPCells cells(num_cells_required, 0);

    // Calculate the band center coordinates in the first
    // column of the multiple alignment. These are calculated by
    // projecting the match diagonal onto the first column. It is possible
    // that these are negative.
    int band_center = start_2 - start_1 + 1;
    int band_origin = band_center - (half_width + 1);
#ifdef DEBUG_EXTEND
    printf("Match start: [%d %d]\n", start_1, start_2);
    printf("Band center, origin: [%d %d]\n", band_center, band_origin);
    printf("Num cells: %zu\n", cells.size());
#endif

    // Fill in the bands column by column
    for(int i = 1; i < num_columns; ++i) {
        int j = band_origin + i; // start row of this band
        int end_row = j + band_width;

        // Trim band coordinates to only compute valid positions
        if(j < 1)
            j = 1;
        if(end_row > num_rows)
            end_row = num_rows;

        if(end_row <= 0 || j >= num_rows || j >= end_row)
            continue; // nothing to do for this column

#ifdef DEBUG_EXTEND
        printf("Filling column %d rows [%d %d]\n", i, j, end_row);
#endif

        // Fill in this band. To avoid the need to perform many tests whether a particular cell
        // is stored in a band, we do some of the calculations outside of the main loop below. 
        // We first calculate the score for the first cell in the band. This calculation cannot
        // involve the cell above the first row so we ignore it below. We then fill in the main
        // part of the band, which can perform valid reads from all its neighboring cells. Finally
        // we calculate the last row, which does not use the cell to its left.

        // Set up initial indices and scores
        int curr_idx = _getBandedCellIndex(i, j, band_width, band_origin);
        int left_idx = _getBandedCellIndex(i - 1, j, band_width, band_origin);
        int diagonal_idx = _getBandedCellIndex(i - 1, j - 1, band_width, band_origin);
        int diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
        int left_score = left_idx != -1 ? cells[left_idx] + GAP_PENALTY : INVALID_SCORE;
        int up_score = 0;

        // Set the first row score
        cells[curr_idx] = std::max(left_score, diagonal_score);

#ifdef DEBUG_EXTEND
        printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
        assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
        assert(diagonal_idx != -1);
#endif

        // Update indices
        curr_idx += 1;
        left_idx += 1;
        diagonal_idx += 1;
        j += 1;

        // Fill in the main part of the band, stopping before the last row
        while(j < end_row - 1) {

#ifdef DEBUG_EXTEND
            assert(diagonal_idx == _getBandedCellIndex(i - 1, j - 1, band_width, band_origin));
            assert(left_idx == _getBandedCellIndex(i - 1, j, band_width, band_origin));
            assert(curr_idx - 1 == _getBandedCellIndex(i, j - 1, band_width, band_origin));
#endif

            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            left_score = cells[left_idx] + GAP_PENALTY;
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = max3(diagonal_score, left_score, up_score);

#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
#endif
            // Update indices
            curr_idx += 1;
            left_idx += 1;
            diagonal_idx += 1;
            j += 1;
        }

        // Fill in last row, here we ignore the left cell which is now out of band
        if(j != end_row) {
            diagonal_score = cells[diagonal_idx] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            up_score = cells[curr_idx - 1] + GAP_PENALTY;
            cells[curr_idx] = std::max(diagonal_score, up_score);
#ifdef DEBUG_EXTEND
            printf("Filled [%d %d] = %d\n", i , j, cells[curr_idx]);
            assert(_getBandedCellIndex(i,j, band_width, band_origin) != -1);
#endif        
        }
    }

    // The location of the highest scoring match in the
    // last row or last column is the maximum scoring overlap
    // for the pair of strings. We start the backtracking from
    // that cell
    int max_row_value = std::numeric_limits<int>::min();
    int max_column_value = std::numeric_limits<int>::min();
    size_t max_row_index = 0;
    size_t max_column_index = 0;

    // Check every column of the last row
    // The first column is skipped to avoid empty alignments
    for(int i = 1; i < num_columns; ++i) {
        int v = _getBandedCellScore(cells, i, num_rows - 1, band_width, band_origin, INVALID_SCORE); 
        if(v > max_row_value) {
            max_row_value = v;
            max_row_index = i;
        }
    }

    // Check every row of the last column
    for(int j = 1; j < num_rows; ++j) {
        int v = _getBandedCellScore(cells, num_columns - 1, j, band_width, band_origin, INVALID_SCORE); 
        if(v > max_column_value) {
            max_column_value = v;
            max_column_index = j;
        }
    }

    // Compute the location at which to start the backtrack
    size_t i;
    size_t j;

    if(max_column_value > max_row_value) {
        i = num_columns - 1;
        j = max_column_index;
        output.score = max_column_value;
    }
    else {
        i = max_row_index;
        j = num_rows - 1;
        output.score = max_row_value;
    }    

#ifdef DEBUG_EXTEND
    printf("BEST: %zu %zu\n", i, j);
#endif

    // Backtrack to fill in the cigar string and alignment start position
    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();
#ifdef DEBUG_EXTEND
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif

    output.edit_distance = 0;
    output.total_columns = 0;

    std::string cigar;
    while(i > 0 && j > 0) {
        // Compute the possible previous locations of the path
        int idx_1 = i - 1;
        int idx_2 = j - 1;

        bool is_match = s1[idx_1] == s2[idx_2];
        int diagonal = _getBandedCellScore(cells, i - 1, j - 1, band_width, band_origin, INVALID_SCORE) + (is_match ? MATCH_SCORE : MISMATCH_PENALTY);
        int up = _getBandedCellScore(cells, i, j - 1, band_width, band_origin, INVALID_SCORE) + GAP_PENALTY;
        int left =  _getBandedCellScore(cells, i -1 , j, band_width, band_origin, INVALID_SCORE) + GAP_PENALTY;
        int curr = _getBandedCellScore(cells, i, j, band_width, band_origin, INVALID_SCORE);

        // If there are multiple possible paths to this cell
        // we break ties in order of insertion,deletion,match
        // this helps left-justify matches for homopolymer runs
        // of unequal lengths
        if(curr == up) {
            cigar.push_back('I');
            j -= 1;
            output.edit_distance += 1;
        } else if(curr == left) {
            cigar.push_back('D');
            i -= 1;
            output.edit_distance += 1;
        } else {
            assert(curr == diagonal);
            if(!is_match)
                output.edit_distance += 1;
            cigar.push_back('M');
            i -= 1;
            j -= 1;
        }

        output.total_columns += 1;
    }

    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;

    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    assert(!cigar.empty());
    output.cigar = compactCigar(cigar);
    return output;
}

static const uint8_t FROM_DIAG = 0;
static const uint8_t FROM_LEFT = 1;
static const uint8_t FROM_UP = 2;
static const uint8_t FROM_INVALID = 3;

// The score for this cell coming from a match, deletion and insertion and its direction
struct AffineCell
{
    AffineCell() : G(0), I(-std::numeric_limits<int>::max()), D(-std::numeric_limits<int>::max()) {}

    //
    int G;
    int I;
    int D;
    
    uint8_t Gt;
    uint8_t It;
    uint8_t Dt;
};

typedef std::vector<AffineCell> AffineCells;
typedef std::vector<AffineCells> AffineMatrix;

SequenceOverlap Overlapper::computeAlignmentAffine(const std::string& s1, const std::string& s2, const OverlapperParams params)
{
    // Exit with invalid intervals if either string is zero length
    SequenceOverlap output;
    if(s1.empty() || s2.empty()) {
        std::cerr << "Overlapper::computeAlignmentAffine error: empty input sequence\n";
        exit(EXIT_FAILURE);
    }

    bool gap_s1_start;
    bool gap_s1_end;
    bool gap_s2_start;
    bool gap_s2_end;
    
    // Set the bools for the explicit alignment types
    if (params.type == ALT_GLOBAL)
    {
        gap_s1_start = false;
        gap_s1_end = false;
        gap_s2_start = false;
        gap_s2_end = false;
    }
    else if (params.type == ALT_OVERLAP)
    {
        gap_s1_start = true;
        gap_s1_end = true;
        gap_s2_start = true;
        gap_s2_end = true;
    }
    else if (params.type == ALT_CONTAINMENT)
    {
        gap_s1_start = true;
        gap_s1_end = true;
        gap_s2_start = false;
        gap_s2_end = false;
    }
    else if (params.type == ALT_CUSTOM)
    {
        gap_s1_start = params.gap_s1_start;
        gap_s1_end = params.gap_s1_end;
        gap_s2_start = params.gap_s2_start;
        gap_s2_end = params.gap_s2_end;
    }
    else
    {
        // Unknown alignment type
        abort();
    }

    // Initialize the scoring matrix
    size_t num_columns = s1.size() + 1;
    size_t num_rows = s2.size() + 1;
    int gap_open = -params.gap_penalty;
    int gap_ext = -params.gap_ext_penalty;
    
    AffineMatrix score_matrix;
    score_matrix.resize(num_columns);
    
    for(size_t i = 0; i < score_matrix.size(); ++i)
        score_matrix[i].resize(num_rows);
    
    // Initialze first row and column
    // Penalties in first row iff gap_s1_start==false
    int c = (gap_s1_start == false ? 1 : 0);
    for(size_t i = 1; i < num_columns; ++i) {
        int v = -(gap_open + i * gap_ext) * c;
        score_matrix[i][0].D = v;
        score_matrix[i][0].Dt = (i == 1? FROM_DIAG : FROM_LEFT);
        score_matrix[i][0].G = v;
        score_matrix[i][0].Gt = FROM_LEFT;
    }

    // Penalties in first column iff gap_s2_start==false
    c = (gap_s2_start == false ? 1 : 0);
    for(size_t j = 1; j < num_rows; ++j) {
        int v = -(gap_open + j * gap_ext) * c;
        score_matrix[0][j].I = v;
        score_matrix[0][j].It = (j == 1? FROM_DIAG : FROM_UP);
        score_matrix[0][j].G = v;
        score_matrix[0][j].Gt = FROM_UP;
    }

    // Calculate scores
    for(size_t i = 1; i < num_columns; ++i) {
        for(size_t j = 1; j < num_rows; ++j) {
            
            // Calculate the score for entry (i,j)
            int idx_1 = i - 1;
            int idx_2 = j - 1;
            int diagonal = score_matrix[i-1][j-1].G + (s1[idx_1] == s2[idx_2] ? params.match_score : params.mismatch_penalty);
            
            AffineCell& curr = score_matrix[i][j];
            AffineCell& up = score_matrix[i][j-1];
            AffineCell& left = score_matrix[i-1][j];
            
            // When computing the score starting from the left/right cells, we have to determine
            // whether to extend an existing gap or start a new one.
            // In the last column, insertion costs are controlled by gap_s2_end
            int ins_open = (i < num_columns - 1 or not gap_s2_end? gap_open : 0);
            int ins_ext = (i < num_columns - 1 or not gap_s2_end? gap_ext : 0);
            
            // In the last row, deletion costs are controlled by gap_s1_end
            int del_open = (j < num_rows - 1 or not gap_s1_end? gap_open : 0);
            int del_ext = (j < num_rows - 1 or not gap_s1_end? gap_ext : 0);
            
            if(up.I > up.G - ins_open) {
                curr.I = up.I - ins_ext;
                curr.It = FROM_UP;
            } else {
                curr.I = up.G - (ins_open + ins_ext);
                curr.It = FROM_DIAG;
            }
            if(left.D > left.G - del_open) {
                curr.D = left.D - del_ext;
                curr.Dt = FROM_LEFT;
            } else {
                curr.D = left.G - (del_open + del_ext);
                curr.Dt = FROM_DIAG;
            }
            
            curr.G = max3(curr.D, curr.I, diagonal);
            if(curr.G == curr.I)
                curr.Gt = FROM_UP;
            else if(curr.G == curr.D)
                curr.Gt = FROM_LEFT;
            else
                curr.Gt = FROM_DIAG;
        }
    }
    
    // With the new scores, the max score is always in the bottom right cell
    size_t i = num_columns - 1;
    size_t j = num_rows - 1;
    
    output.score = score_matrix[i][j].G;
    uint8_t direction = score_matrix[i][j].Gt;
    
    // However, the alignment might contain free end gaps which we now remove
    if (gap_s2_end)
    {
        while (j >= 1 and direction == FROM_UP)
        {
            direction = score_matrix[i][j].It;
            --j;
            if (direction == FROM_DIAG)
                direction = score_matrix[i][j].Gt;
        }
    }

    if (gap_s1_end and j == num_rows - 1)
    {
        while (i >= 1 and direction == FROM_LEFT)
        {
            direction = score_matrix[i][j].Dt;
            --i;
            if (direction == FROM_DIAG)
                direction = score_matrix[i][j].Gt;
        }
    }

    // Set the alignment endpoints to be the index of the last aligned base
    output.match[0].end = i - 1;
    output.match[1].end = j - 1;
    output.length[0] = s1.length();
    output.length[1] = s2.length();

#ifdef DEBUG_OVERLAPPER
    printf("Endpoints selected: (%d %d) with score %d\n", output.match[0].end, output.match[1].end, output.score);
#endif
    
    output.edit_distance = 0;
    output.total_columns = 0;
    std::string cigar;
    
    // We stop when we hit an edge along which gaps are free
    while (not (i == 0 and j == 0) // absolute stop, regardless of free gaps
            and not (i == 0 and gap_s2_start) // stop at left edge if s2 start gaps are free
            and not (j == 0 and gap_s1_start)) // stop at top edge if s1 start gaps are free
    {
        if (direction == FROM_UP)
        {
            cigar.push_back('I');
            ++output.edit_distance;
            direction = score_matrix[i][j].It;
            --j;

            if (direction == FROM_DIAG)
                direction = score_matrix[i][j].Gt;
        }
        else if (direction == FROM_LEFT)
        {
            cigar.push_back('D');
            ++output.edit_distance;
            direction = score_matrix[i][j].Dt;
            --i;
         
            if (direction == FROM_DIAG)
                direction = score_matrix[i][j].Gt;
        }
        else
        {
            assert(direction == FROM_DIAG);
            assert(i > 0);
            assert(j > 0);
            
            bool is_match = s1[i - 1] == s2[j - 1];
            
            if (is_match)
            {
                cigar.push_back(params.use_m_ops? 'M' : '=');
            }
            else
            {
                cigar.push_back(params.use_m_ops? 'M' : 'X');
                ++output.edit_distance;
            }
            --i;
            --j;
            direction = score_matrix[i][j].Gt;
        }
        ++output.total_columns;
    }
    
    // Set the alignment startpoints
    output.match[0].start = i;
    output.match[1].start = j;
    // Compact the expanded cigar string into the canonical run length encoding
    // The backtracking produces a cigar string in reversed order, flip it
    std::reverse(cigar.begin(), cigar.end());
    output.cigar = compactCigar(cigar);
    return output;
}

// Compact an expanded CIGAR string into a regular cigar string
std::string Overlapper::compactCigar(const std::string& ecigar)
{
    if(ecigar.empty())
        return "";

    std::stringstream compact_cigar;
    char curr_symbol = ecigar[0];
    int curr_run = 1;
    for(size_t i = 1; i < ecigar.size(); ++i) {
        if(ecigar[i] == curr_symbol) {
            curr_run += 1;
        } else {
            compact_cigar << curr_run << curr_symbol;
            curr_symbol = ecigar[i];
            curr_run = 1;
        }
    }

    // Add last symbol/run
    compact_cigar << curr_run << curr_symbol;
    return compact_cigar.str();
}
