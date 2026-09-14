/**
 * Needleman-Wunsch with Affine Gap Penalties (Gotoh's Algorithm)
 *
 * Global sequence alignment using affine gap penalties: gap_open + k * gap_extend
 * for a gap of length k. Gotoh's (1982) O(mn)-time reformulation of the naive
 * O(mn * max(m,n))-time affine-gap extension of Needleman-Wunsch, achieved by
 * tracking three DP matrices instead of one.
 *
 * - M[i][j]: best global-alignment score ending in a match/mismatch
 * - I[i][j]: best score ending with a gap in seq1 (insertion, consumes seq2)
 * - D[i][j]: best score ending with a gap in seq2 (deletion, consumes seq1)
 *
 * Unlike the local (Smith-Waterman) affine variant, there is no flooring at 0
 * and traceback always runs to (0, 0).
 *
 * Time Complexity: O(m * n)
 * Space Complexity: O(m * n)
 *
 * Reference:
 * Gotoh, O. (1982). An improved algorithm for matching biological sequences.
 * Journal of Molecular Biology, 162(3), 705-708.
 * https://doi.org/10.1016/0022-2836(82)90398-9
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
    double score;
};

/**
 * Perform global sequence alignment with affine gap penalties.
 *
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param match_score Score for matching characters (default: 1)
 * @param mismatch_penalty Penalty for mismatching characters (default: -1)
 * @param gap_open Additional penalty charged once when a gap is opened (default: -3)
 * @param gap_extend Penalty charged per gap character (default: -1)
 * @return Alignment structure containing aligned sequences and score
 */
Alignment needleman_wunsch_affine(const string& seq1, const string& seq2,
                                   double match_score = 1,
                                   double mismatch_penalty = -1,
                                   double gap_open = -3,
                                   double gap_extend = -1) {
    int m = seq1.length();
    int n = seq2.length();
    const double NEG_INF = -numeric_limits<double>::infinity();

    vector<vector<double>> M(m + 1, vector<double>(n + 1, NEG_INF));
    vector<vector<double>> I(m + 1, vector<double>(n + 1, NEG_INF));
    vector<vector<double>> D(m + 1, vector<double>(n + 1, NEG_INF));

    // M[0][j] (j>0) and M[i][0] (i>0) stay -inf: a match/mismatch cannot end
    // an alignment that consumed zero characters from the other sequence.
    M[0][0] = 0;
    for (int j = 1; j <= n; j++) I[0][j] = gap_open + j * gap_extend;
    for (int i = 1; i <= m; i++) D[i][0] = gap_open + i * gap_extend;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            I[i][j] = max({M[i][j - 1] + gap_open + gap_extend,
                           I[i][j - 1] + gap_extend,
                           D[i][j - 1] + gap_open + gap_extend});
            D[i][j] = max({M[i - 1][j] + gap_open + gap_extend,
                           D[i - 1][j] + gap_extend,
                           I[i - 1][j] + gap_open + gap_extend});
            double s = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;
            M[i][j] = max({M[i - 1][j - 1] + s, I[i - 1][j - 1] + s, D[i - 1][j - 1] + s});
        }
    }

    double final_score = max({M[m][n], I[m][n], D[m][n]});
    char current_matrix = (final_score == M[m][n]) ? 'M' : (final_score == I[m][n]) ? 'I' : 'D';

    string aligned_seq1, aligned_seq2;
    int i = m, j = n;

    while (i > 0 || j > 0) {
        if (current_matrix == 'M') {
            double s = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            if (M[i][j] == M[i - 1][j - 1] + s) current_matrix = 'M';
            else if (M[i][j] == I[i - 1][j - 1] + s) current_matrix = 'I';
            else current_matrix = 'D';
            i--; j--;
        } else if (current_matrix == 'I') {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            if (I[i][j] == I[i][j - 1] + gap_extend) current_matrix = 'I';
            else if (I[i][j] == M[i][j - 1] + gap_open + gap_extend) current_matrix = 'M';
            else current_matrix = 'D';
            j--;
        } else { // current_matrix == 'D'
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            if (D[i][j] == D[i - 1][j] + gap_extend) current_matrix = 'D';
            else if (D[i][j] == M[i - 1][j] + gap_open + gap_extend) current_matrix = 'M';
            else current_matrix = 'I';
            i--;
        }
    }

    return {aligned_seq1, aligned_seq2, final_score};
}

int main() {
    string seq1 = "GGTTGACTA";
    string seq2 = "TGTTACGG";

    cout << "Needleman-Wunsch with Affine Gap Penalties (Gotoh 1982) - Global Alignment" << endl;
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl << endl;

    Alignment result = needleman_wunsch_affine(seq1, seq2, 2, -1, -3, -1);
    cout << "With affine gaps (gap_open=-3, gap_extend=-1):" << endl;
    cout << "Aligned Sequence 1: " << result.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result.aligned_seq2 << endl;
    cout << "Alignment Score: " << result.score << endl << endl;

    string seq3 = "ACGTACGTACGT";
    string seq4 = "ACGTAAACGT";
    cout << "Example 2 - clustered indel:" << endl;
    cout << "Sequence 1: " << seq3 << endl;
    cout << "Sequence 2: " << seq4 << endl;
    Alignment result2 = needleman_wunsch_affine(seq3, seq4, 2, -1, -5, -1);
    cout << "Aligned Sequence 1: " << result2.aligned_seq1 << endl;
    cout << "Aligned Sequence 2: " << result2.aligned_seq2 << endl;
    cout << "Alignment Score: " << result2.score << endl;

    return 0;
}
