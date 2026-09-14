/**
 * Wavefront Alignment Algorithm (WFA) - Gap-Affine Exact Alignment
 *
 * WFA is an exact, provably-optimal gap-affine global-alignment algorithm
 * that runs in O(n * s + s^2) time and O(s^2) memory, where n = max(len(seq1),
 * len(seq2)) and s is the optimal alignment score. For similar sequences, s
 * is small, so WFA is near-linear in practice - unlike the O(mn) of
 * Needleman-Wunsch/Gotoh, which is oblivious to sequence similarity.
 *
 * WFA is indexed by score rather than by sequence position: for each
 * candidate score s = 0, 1, 2, ..., it tracks, for every diagonal
 * k = j - i, the furthest-reaching offset i reachable with exactly that
 * score (a "wavefront"). Matching characters along a diagonal are free
 * ("greedy extension"). The algorithm stops the moment a wavefront reaches
 * the bottom-right corner, so it never explores alignments worse than optimal.
 *
 * Uses a MINIMIZATION convention (mismatch/gap penalties are positive costs),
 * matching the original paper and its WFA2-lib reference implementation.
 *
 * Time Complexity: O(n * s + s^2) -- near-linear for similar sequences (small s)
 * Space Complexity: O(s^2) for this straightforward implementation
 *
 * Reference:
 * Marco-Sola, S., Moure, J. C., Moreto, M., & Espinosa, A. (2021). Fast
 * gap-affine pairwise alignment using the wavefront algorithm. Bioinformatics,
 * 37(4), 456-463. https://doi.org/10.1093/bioinformatics/btaa777
 * Code: https://github.com/smarco/WFA2-lib
 */

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>

using namespace std;

using Wavefront = unordered_map<int, int>; // diagonal -> furthest seq1-offset i

static int extend_wavefront(int diag, const string& seq1, const string& seq2, int offset) {
    int m = seq1.size(), n = seq2.size();
    int i = offset, j = offset + diag;
    while (i < m && j < n && seq1[i] == seq2[j]) {
        i++; j++;
    }
    return i;
}

/**
 * Compute the optimal gap-affine alignment score between seq1 and seq2 using
 * the Wavefront Alignment algorithm. Costs are non-negative penalties
 * (0 = match); this function MINIMIZES total cost.
 */
int wavefront_alignment(const string& seq1, const string& seq2,
                         int mismatch = 4, int gap_open = 6, int gap_extend = 2,
                         int max_score = -1) {
    int m = seq1.size(), n = seq2.size();
    if (max_score < 0) {
        max_score = (m + n) * max(mismatch, gap_open + gap_extend) + 1;
    }

    int target_diag = n - m;
    int min_diag = -m, max_diag = n;

    auto in_bounds = [&](int k, int i) {
        return k >= min_diag && k <= max_diag && i >= 0 && i <= m && (i + k) >= 0 && (i + k) <= n;
    };

    vector<Wavefront> M(max_score + 1), I(max_score + 1), D(max_score + 1);
    vector<bool> hasM(max_score + 1, false), hasI(max_score + 1, false), hasD(max_score + 1, false);

    M[0][0] = extend_wavefront(0, seq1, seq2, 0);
    hasM[0] = true;

    if (m == n && M[0][0] == m) return 0;

    for (int s = 1; s <= max_score; s++) {
        Wavefront m_s, i_s, d_s;

        // Insertion: consumes a seq2 char (j += 1) => diag += 1
        int cost_open = gap_open + gap_extend;
        if (s - cost_open >= 0 && hasM[s - cost_open]) {
            for (auto& [k, i] : M[s - cost_open]) {
                int nk = k + 1;
                if (in_bounds(nk, i) && (!i_s.count(nk) || i > i_s[nk])) i_s[nk] = i;
            }
        }
        if (s - gap_extend >= 0 && hasI[s - gap_extend]) {
            for (auto& [k, i] : I[s - gap_extend]) {
                int nk = k + 1;
                if (in_bounds(nk, i) && (!i_s.count(nk) || i > i_s[nk])) i_s[nk] = i;
            }
        }

        // Deletion: consumes a seq1 char (i += 1) => diag -= 1
        if (s - cost_open >= 0 && hasM[s - cost_open]) {
            for (auto& [k, i] : M[s - cost_open]) {
                int nk = k - 1, cand = i + 1;
                if (in_bounds(nk, cand) && (!d_s.count(nk) || cand > d_s[nk])) d_s[nk] = cand;
            }
        }
        if (s - gap_extend >= 0 && hasD[s - gap_extend]) {
            for (auto& [k, i] : D[s - gap_extend]) {
                int nk = k - 1, cand = i + 1;
                if (in_bounds(nk, cand) && (!d_s.count(nk) || cand > d_s[nk])) d_s[nk] = cand;
            }
        }

        // Match/mismatch wavefront
        Wavefront empty_wf;
        const Wavefront& prev_mismatch = (s - mismatch >= 0 && hasM[s - mismatch]) ? M[s - mismatch] : empty_wf;

        vector<int> diagonals;
        for (auto& [k, _] : prev_mismatch) diagonals.push_back(k);
        for (auto& [k, _] : i_s) diagonals.push_back(k);
        for (auto& [k, _] : d_s) diagonals.push_back(k);
        sort(diagonals.begin(), diagonals.end());
        diagonals.erase(unique(diagonals.begin(), diagonals.end()), diagonals.end());

        for (int k : diagonals) {
            int best = -1;
            if (prev_mismatch.count(k) && in_bounds(k, prev_mismatch.at(k) + 1)) {
                best = max(best, prev_mismatch.at(k) + 1);
            }
            if (i_s.count(k)) best = max(best, i_s[k]);
            if (d_s.count(k)) best = max(best, d_s[k]);
            if (best >= 0) {
                m_s[k] = extend_wavefront(k, seq1, seq2, best);
            }
        }

        if (!i_s.empty()) { I[s] = i_s; hasI[s] = true; }
        if (!d_s.empty()) { D[s] = d_s; hasD[s] = true; }
        if (!m_s.empty()) { M[s] = m_s; hasM[s] = true; }

        if (m_s.count(target_diag) && m_s[target_diag] == m) {
            int j = m_s[target_diag] + target_diag;
            if (j == n) return s;
        }
    }

    throw runtime_error("max_score exceeded; sequences may be too dissimilar");
}

int main() {
    cout << "Wavefront Alignment Algorithm (WFA) - Marco-Sola et al. 2021" << endl;
    cout << string(60, '=') << endl;

    {
        string seq1 = "GATTACA", seq2 = "GATCACA";
        int score = wavefront_alignment(seq1, seq2);
        cout << "\nSequence 1: " << seq1 << "\nSequence 2: " << seq2 << endl;
        cout << "Optimal gap-affine edit score: " << score << " (single substitution)" << endl;
    }

    {
        string seq1 = "ACGTACGTACGT", seq2 = "ACGTAAACGT";
        int score = wavefront_alignment(seq1, seq2);
        cout << "\nSequence 1: " << seq1 << "\nSequence 2: " << seq2 << endl;
        cout << "Optimal gap-affine edit score: " << score << endl;
    }

    cout << "\nWhy WFA matters: for near-identical long reads (large n, small s)," << endl;
    cout << "its O(n*s + s^2) cost is close to linear, vs. O(n^2) for classic DP." << endl;

    return 0;
}
