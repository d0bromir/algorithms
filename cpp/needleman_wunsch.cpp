/**
 * Needleman-Wunsch Algorithm for Global Sequence Alignment
 * 
 * This is a dynamic programming algorithm for aligning two sequences globally.
 * It finds the optimal alignment between two sequences by maximizing the alignment score.
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
 * Perform global sequence alignment using Needleman-Wunsch algorithm.
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters (default: 1)
 * @param mismatch_penalty Penalty for mismatching characters (default: -1)
 * @param gap_penalty Penalty for gaps (default: -1)
 * @return Alignment structure containing aligned sequences and score
 */
Alignment needleman_wunsch(const string& seq1, const string& seq2, 
                          int match_score = 1, 
                          int mismatch_penalty = -1, 
                          int gap_penalty = -1) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Initialize DP matrix
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    
    // Initialize first row and column with gap penalties
    for (int i = 0; i <= m; i++) {
        dp[i][0] = i * gap_penalty;
    }
    for (int j = 0; j <= n; j++) {
        dp[0][j] = j * gap_penalty;
    }
    
    // Fill DP matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match = dp[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
            int del = dp[i-1][j] + gap_penalty;
            int insert = dp[i][j-1] + gap_penalty;
            
            dp[i][j] = max({match, del, insert});
        }
    }
    
    // Traceback to find alignment
    string aligned_seq1 = "";
    string aligned_seq2 = "";
    int i = m, j = n;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int current_score = dp[i][j];
            int diagonal_score = dp[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
            
            if (current_score == diagonal_score) {
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                i--;
                j--;
            } else if (i > 0 && current_score == dp[i-1][j] + gap_penalty) {
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = "-" + aligned_seq2;
                i--;
            } else {
                aligned_seq1 = "-" + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                j--;
            }
        } else if (i > 0) {
            aligned_seq1 = seq1[i-1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
        } else {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j-1] + aligned_seq2;
            j--;
        }
    }
    
    return {aligned_seq1, aligned_seq2, dp[m][n]};
}

int main() {
    // Example usage
    string seq1 = "GATTACA";
    string seq2 = "GCATGCU";
    
    Alignment result = needleman_wunsch(seq1, seq2);
    
    cout << "Needleman-Wunsch Global Alignment" << endl;
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "\nAligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl;
    
    return 0;
}
