/**
 * Smith-Waterman Algorithm for Local Sequence Alignment
 * 
 * This is a dynamic programming algorithm for performing local sequence alignment.
 * It finds the optimal local alignment between two sequences by maximizing the alignment score.
 * 
 * Time Complexity: O(m * n) where m and n are the lengths of the sequences
 * Space Complexity: O(m * n)
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Alignment {
    string aligned_seq1;
    string aligned_seq2;
    int score;
};

/**
 * Perform local sequence alignment using Smith-Waterman algorithm.
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters (default: 2)
 * @param mismatch_penalty Penalty for mismatching characters (default: -1)
 * @param gap_penalty Penalty for gaps (default: -1)
 * @return Alignment structure containing aligned sequences and score
 */
Alignment smith_waterman(const string& seq1, const string& seq2,
                        int match_score = 2,
                        int mismatch_penalty = -1,
                        int gap_penalty = -1) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Initialize DP matrix with zeros (key difference from Needleman-Wunsch)
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    // Fill DP matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match = dp[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
            int del = dp[i-1][j] + gap_penalty;
            int insert = dp[i][j-1] + gap_penalty;
            
            // Key difference: allow zero (no negative scores)
            dp[i][j] = max({0, match, del, insert});
            
            // Track maximum score position
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback from maximum score position
    string aligned_seq1 = "";
    string aligned_seq2 = "";
    int i = max_i, j = max_j;
    
    while (i > 0 && j > 0 && dp[i][j] > 0) {
        int current_score = dp[i][j];
        int diagonal_score = dp[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
        
        if (current_score == diagonal_score) {
            aligned_seq1 = seq1[i-1] + aligned_seq1;
            aligned_seq2 = seq2[j-1] + aligned_seq2;
            i--;
            j--;
        } else if (current_score == dp[i-1][j] + gap_penalty) {
            aligned_seq1 = seq1[i-1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
        } else {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j-1] + aligned_seq2;
            j--;
        }
    }
    
    return {aligned_seq1, aligned_seq2, max_score};
}

int main() {
    // Example usage
    string seq1 = "GGTTGACTA";
    string seq2 = "TGTTACGG";
    
    Alignment result = smith_waterman(seq1, seq2);
    
    cout << "Smith-Waterman Local Alignment" << endl;
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "\nAligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl;
    
    return 0;
}
