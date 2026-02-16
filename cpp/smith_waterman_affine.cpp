/**
 * Smith-Waterman Algorithm with Affine Gap Penalties for Local Sequence Alignment
 * 
 * This is an advanced version of the Smith-Waterman algorithm that uses affine gap penalties.
 * Affine gap penalty = gap_open + k * gap_extend (where k is the gap length)
 * This is more biologically realistic as it penalizes gap opening more than gap extension.
 * 
 * Uses three matrices:
 * - M[i][j]: optimal score ending with match/mismatch
 * - I[i][j]: optimal score ending with insertion (gap in seq1)
 * - D[i][j]: optimal score ending with deletion (gap in seq2)
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
 * Perform local sequence alignment using Smith-Waterman algorithm with affine gap penalties.
 * 
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters (default: 2)
 * @param mismatch_penalty Penalty for mismatching characters (default: -1)
 * @param gap_open Penalty for opening a gap (default: -3)
 * @param gap_extend Penalty for extending a gap (default: -1)
 * @return Alignment structure containing aligned sequences and score
 */
Alignment smith_waterman_affine(const string& seq1, const string& seq2,
                                int match_score = 2,
                                int mismatch_penalty = -1,
                                int gap_open = -3,
                                int gap_extend = -1) {
    int m = seq1.length();
    int n = seq2.length();
    
    const int NEG_INF = -1000000;
    
    // Initialize three matrices: M (match), I (insert), D (delete)
    vector<vector<int>> M(m + 1, vector<int>(n + 1, 0));
    vector<vector<int>> I(m + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<int>> D(m + 1, vector<int>(n + 1, NEG_INF));
    
    int max_score = 0;
    int max_i = 0, max_j = 0;
    char max_matrix = 'M';
    
    // Fill DP matrices
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            // Calculate I[i][j] - gap in seq2 (insertion)
            I[i][j] = max({0,
                          M[i][j-1] + gap_open + gap_extend,
                          I[i][j-1] + gap_extend,
                          D[i][j-1] + gap_open + gap_extend});
            
            // Calculate D[i][j] - gap in seq1 (deletion)
            D[i][j] = max({0,
                          M[i-1][j] + gap_open + gap_extend,
                          I[i-1][j] + gap_open + gap_extend,
                          D[i-1][j] + gap_extend});
            
            // Calculate M[i][j] - match/mismatch
            int match_mismatch_score = (seq1[i-1] == seq2[j-1]) ? match_score : mismatch_penalty;
            M[i][j] = max({0,
                          M[i-1][j-1] + match_mismatch_score,
                          I[i-1][j-1] + match_mismatch_score,
                          D[i-1][j-1] + match_mismatch_score});
            
            // Track maximum score across all three matrices
            if (M[i][j] >= max_score) {
                max_score = M[i][j];
                max_i = i;
                max_j = j;
                max_matrix = 'M';
            }
            if (I[i][j] >= max_score) {
                max_score = I[i][j];
                max_i = i;
                max_j = j;
                max_matrix = 'I';
            }
            if (D[i][j] >= max_score) {
                max_score = D[i][j];
                max_i = i;
                max_j = j;
                max_matrix = 'D';
            }
        }
    }
    
    // Traceback from maximum score position
    string aligned_seq1 = "";
    string aligned_seq2 = "";
    int i = max_i, j = max_j;
    char current_matrix = max_matrix;
    
    while (i > 0 && j > 0) {
        if (current_matrix == 'M') {
            if (M[i][j] <= 0) break;
            
            int match_mismatch_score = (seq1[i-1] == seq2[j-1]) ? match_score : mismatch_penalty;
            
            if (M[i][j] == M[i-1][j-1] + match_mismatch_score) {
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                i--;
                j--;
                current_matrix = 'M';
            } else if (M[i][j] == I[i-1][j-1] + match_mismatch_score) {
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                i--;
                j--;
                current_matrix = 'I';
            } else if (M[i][j] == D[i-1][j-1] + match_mismatch_score) {
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                i--;
                j--;
                current_matrix = 'D';
            } else {
                break;
            }
        } else if (current_matrix == 'I') {
            if (I[i][j] <= 0) break;
            
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j-1] + aligned_seq2;
            j--;
            
            if (I[i][j] == M[i][j-1] + gap_open + gap_extend) {
                current_matrix = 'M';
            } else if (I[i][j] == I[i][j-1] + gap_extend) {
                current_matrix = 'I';
            } else if (I[i][j] == D[i][j-1] + gap_open + gap_extend) {
                current_matrix = 'D';
            } else {
                break;
            }
        } else { // current_matrix == 'D'
            if (D[i][j] <= 0) break;
            
            aligned_seq1 = seq1[i-1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
            
            if (D[i][j] == M[i-1][j] + gap_open + gap_extend) {
                current_matrix = 'M';
            } else if (D[i][j] == D[i-1][j] + gap_extend) {
                current_matrix = 'D';
            } else if (D[i][j] == I[i-1][j] + gap_open + gap_extend) {
                current_matrix = 'I';
            } else {
                break;
            }
        }
    }
    
    return {aligned_seq1, aligned_seq2, max_score};
}

int main() {
    // Example usage
    string seq1 = "GGTTGACTA";
    string seq2 = "TGTTACGG";
    
    cout << "Smith-Waterman with Affine Gap Penalties - Local Alignment" << endl;
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << endl;
    
    // Compare linear vs affine gap penalties
    cout << "With affine gaps (gap_open=-3, gap_extend=-1):" << endl;
    Alignment result_affine = smith_waterman_affine(seq1, seq2, 2, -1, -3, -1);
    cout << "Aligned Sequence 1: " << result_affine.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result_affine.aligned_seq2 << endl;
    cout << "Alignment Score: " << result_affine.score << endl;
    cout << endl;
    
    // Example with longer sequences showing affine gap advantage
    string seq3 = "ACGTACGTACGT";
    string seq4 = "ACGTAAACGT";
    
    cout << "Example 2 - Longer sequences:" << endl;
    cout << "Sequence 1: " << seq3 << endl;
    cout << "Sequence 2: " << seq4 << endl;
    
    Alignment result2 = smith_waterman_affine(seq3, seq4, 2, -1, -5, -1);
    cout << "\nAligned Sequence 1: " << result2.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result2.aligned_seq2 << endl;
    cout << "Alignment Score: " << result2.score << endl;
    
    return 0;
}
