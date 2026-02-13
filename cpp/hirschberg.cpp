/**
 * Hirschberg Algorithm for Space-Efficient Global Sequence Alignment
 * 
 * This is a space-efficient divide-and-conquer algorithm for global sequence alignment,
 * developed by Dan Hirschberg in 1975. It improves upon the Needleman-Wunsch algorithm
 * by reducing space complexity from O(m * n) to O(min(m, n)) while maintaining the
 * same O(m * n) time complexity.
 * 
 * The algorithm works by:
 * 1. Using only two rows of the DP matrix at a time (space optimization)
 * 2. Dividing the problem recursively at the midpoint
 * 3. Finding the optimal split point using NW score function
 * 4. Recursively aligning left and right halves
 * 5. Concatenating the results
 * 
 * Time Complexity: O(m * n) where m and n are the lengths of the sequences
 * Space Complexity: O(min(m, n)) - key improvement over standard Needleman-Wunsch
 * 
 * Reference:
 * Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common
 * subsequences. Communications of the ACM, 18(6), 341-343.
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

using namespace std;

struct Alignment {
    string aligned_seq1;
    string aligned_seq2;
    int score;
};

/**
 * Compute the last row of Needleman-Wunsch scores using only O(n) space.
 * 
 * This function calculates the alignment scores but only keeps the last row
 * of the DP matrix, which is sufficient for finding the optimal split point
 * in Hirschberg's algorithm.
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters
 * @param mismatch_penalty Penalty for mismatching characters
 * @param gap_penalty Penalty for gaps
 * @return Vector containing the last row of DP scores (length n+1)
 */
vector<int> nw_score(const string& seq1, const string& seq2,
                     int match_score, int mismatch_penalty, int gap_penalty) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Only maintain two rows: previous and current
    vector<int> prev_row(n + 1);
    vector<int> curr_row(n + 1);
    
    // Initialize first row
    for (int j = 0; j <= n; j++) {
        prev_row[j] = j * gap_penalty;
    }
    
    // Fill rows one at a time
    for (int i = 1; i <= m; i++) {
        curr_row[0] = i * gap_penalty;
        
        for (int j = 1; j <= n; j++) {
            int match = prev_row[j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_penalty);
            int del = prev_row[j] + gap_penalty;
            int insert = curr_row[j - 1] + gap_penalty;
            
            curr_row[j] = max({match, del, insert});
        }
        
        // Swap rows for next iteration
        swap(prev_row, curr_row);
    }
    
    return prev_row;
}

/**
 * Perform space-efficient global sequence alignment using Hirschberg algorithm.
 * 
 * This is a divide-and-conquer algorithm that produces the same optimal alignment
 * as Needleman-Wunsch but uses only O(min(m, n)) space instead of O(m * n).
 * 
 * The algorithm recursively:
 * 1. Finds the midpoint of the first sequence
 * 2. Computes NW scores from both ends to find optimal split in second sequence
 * 3. Recursively aligns the left halves and right halves
 * 4. Concatenates the results
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters (default: 1)
 * @param mismatch_penalty Penalty for mismatching characters (default: -1)
 * @param gap_penalty Penalty for gaps (default: -1)
 * @return Alignment structure containing aligned sequences and score
 */
Alignment hirschberg(const string& seq1, const string& seq2,
                    int match_score = 1,
                    int mismatch_penalty = -1,
                    int gap_penalty = -1) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Base cases
    if (m == 0) {
        // All gaps in seq1
        return {string(n, '-'), seq2, n * gap_penalty};
    }
    
    if (n == 0) {
        // All gaps in seq2
        return {seq1, string(m, '-'), m * gap_penalty};
    }
    
    if (m == 1) {
        // Single character in seq1 - use simple alignment
        // Try to match seq1[0] with each position in seq2
        int best_score = numeric_limits<int>::min();
        int best_j = 0;
        
        for (int j = 0; j <= n; j++) {
            // Calculate score for aligning seq1[0] at position j in seq2
            int score = j * gap_penalty;  // gaps before
            if (j < n) {
                if (seq1[0] == seq2[j]) {
                    score += match_score;
                } else {
                    score += mismatch_penalty;
                }
                score += (n - j - 1) * gap_penalty;  // gaps after
            } else {
                score += gap_penalty;  // seq1[0] aligned to gap
            }
            
            if (score > best_score) {
                best_score = score;
                best_j = j;
            }
        }
        
        // Build alignment based on best position
        string aligned1, aligned2;
        if (best_j == n) {
            // Align seq1[0] to gap at end
            aligned1 = string(n, '-') + seq1[0];
            aligned2 = seq2 + "-";
        } else {
            // Align seq1[0] to seq2[best_j]
            aligned1 = string(best_j, '-') + seq1[0] + string(n - best_j - 1, '-');
            aligned2 = seq2;
        }
        
        return {aligned1, aligned2, best_score};
    }
    
    if (n == 1) {
        // Single character in seq2 - use simple alignment
        int best_score = numeric_limits<int>::min();
        int best_i = 0;
        
        for (int i = 0; i <= m; i++) {
            int score = i * gap_penalty;
            if (i < m) {
                if (seq1[i] == seq2[0]) {
                    score += match_score;
                } else {
                    score += mismatch_penalty;
                }
                score += (m - i - 1) * gap_penalty;
            } else {
                score += gap_penalty;
            }
            
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        
        string aligned1, aligned2;
        if (best_i == m) {
            aligned1 = seq1 + "-";
            aligned2 = string(m, '-') + seq2[0];
        } else {
            aligned1 = seq1;
            aligned2 = string(best_i, '-') + seq2[0] + string(m - best_i - 1, '-');
        }
        
        return {aligned1, aligned2, best_score};
    }
    
    // Divide and conquer
    // Find the midpoint of seq1
    int mid = m / 2;
    
    // Compute NW scores from left (seq1[:mid] vs seq2)
    string seq1_left = seq1.substr(0, mid);
    vector<int> score_left = nw_score(seq1_left, seq2, match_score, mismatch_penalty, gap_penalty);
    
    // Compute NW scores from right (seq1[mid:][::-1] vs seq2[::-1])
    // We reverse both sequences to compute scores from the end
    string seq1_right = seq1.substr(mid);
    reverse(seq1_right.begin(), seq1_right.end());
    string seq2_rev = seq2;
    reverse(seq2_rev.begin(), seq2_rev.end());
    vector<int> score_right = nw_score(seq1_right, seq2_rev, match_score, mismatch_penalty, gap_penalty);
    
    // Reverse score_right to align with seq2
    reverse(score_right.begin(), score_right.end());
    
    // Find the split point in seq2 that maximizes total score
    int max_score = numeric_limits<int>::min();
    int split = 0;
    
    for (int j = 0; j <= n; j++) {
        int total_score = score_left[j] + score_right[j];
        if (total_score > max_score) {
            max_score = total_score;
            split = j;
        }
    }
    
    // Recursively align left and right parts
    Alignment left = hirschberg(seq1.substr(0, mid), seq2.substr(0, split),
                               match_score, mismatch_penalty, gap_penalty);
    
    Alignment right = hirschberg(seq1.substr(mid), seq2.substr(split),
                                match_score, mismatch_penalty, gap_penalty);
    
    // Combine the results
    string aligned1 = left.aligned_seq1 + right.aligned_seq1;
    string aligned2 = left.aligned_seq2 + right.aligned_seq2;
    int total_score = left.score + right.score;
    
    return {aligned1, aligned2, total_score};
}

int main() {
    // Example usage demonstrating the Hirschberg algorithm
    cout << "Hirschberg Space-Efficient Global Alignment" << endl;
    cout << string(60, '=') << endl;
    
    // Example 1: DNA sequence alignment
    cout << "\nExample 1: DNA Sequence Alignment" << endl;
    cout << string(60, '-') << endl;
    string seq1 = "GATTACA";
    string seq2 = "GCATGCU";
    
    Alignment result = hirschberg(seq1, seq2);
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "\nAligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl;
    
    // Example 2: Custom scoring parameters
    cout << "\n\nExample 2: Custom Scoring Parameters" << endl;
    cout << string(60, '-') << endl;
    seq1 = "HEAGAWGHEE";
    seq2 = "PAWHEAE";
    
    result = hirschberg(seq1, seq2, 2, -1, -2);
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "\nAligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl;
    cout << "Parameters: match=2, mismatch=-1, gap=-2" << endl;
    
    // Example 3: Identical sequences
    cout << "\n\nExample 3: Identical Sequences" << endl;
    cout << string(60, '-') << endl;
    seq1 = "ACGTACGT";
    seq2 = "ACGTACGT";
    
    result = hirschberg(seq1, seq2);
    
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "\nAligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl;
    
    // Example 4: Space complexity comparison
    cout << "\n\nSpace Complexity Comparison" << endl;
    cout << string(60, '-') << endl;
    int m = 1000, n = 1000;
    cout << "For sequences of length " << m << " and " << n << ":" << endl;
    cout << "  Needleman-Wunsch space: O(" << m << " x " << n << ") = ~" << (m * n) << " integers" << endl;
    cout << "  Hirschberg space: O(min(" << m << ", " << n << ")) = ~" << min(m, n) << " integers" << endl;
    cout << "  Space reduction: ~" << ((m * n) / min(m, n)) << "x improvement" << endl;
    
    return 0;
}
