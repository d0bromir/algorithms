/**
 * Myers' Bit-Vector Algorithm for Approximate String Matching / Edit Distance
 *
 * Computes the edit (Levenshtein) distance between a pattern and every prefix
 * of a text by packing a full column of the O(mn) edit-distance DP matrix
 * into a machine word and updating it with O(1) bitwise operations per
 * column.
 *
 * This version packs the vectors into a single 64-bit word, so it directly
 * handles patterns up to 64 characters (the common case for k-mer / seed
 * length in short-read aligners); for longer patterns, real implementations
 * (e.g. edlib) tile this across ceil(m / 64) words.
 *
 * Time Complexity: O(n) for patterns of length <= 64 (O(n * ceil(m/64)) in general)
 * Space Complexity: O(sigma) for the character bitmask table (Peq), O(1) state
 *
 * Reference:
 * Myers, G. (1999). A fast bit-vector algorithm for approximate string matching
 * based on dynamic programming. Journal of the ACM, 46(3), 395-415.
 * https://doi.org/10.1145/316542.316550
 *
 * See also:
 * Sosic, M., & Sikic, M. (2017). Edlib: a C/C++ library for fast, exact
 * sequence alignment using edit distance. Bioinformatics, 33(9), 1394-1395.
 * https://doi.org/10.1093/bioinformatics/btw753
 */

#include <iostream>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

struct MatchResult {
    int edit_distance;
    int end_position; // text[:end_position] achieves edit_distance
};

/**
 * Compute the minimal edit distance of `pattern` against any prefix of
 * `text`, using Myers' bit-vector algorithm with a single 64-bit word.
 *
 * @param pattern Pattern string, length m <= 64
 * @param text Text string to scan, length n
 */
MatchResult myers_bit_vector_edit_distance(const string& pattern, const string& text) {
    const int m = static_cast<int>(pattern.size());
    const int n = static_cast<int>(text.size());
    if (m == 0) return {n, n};
    if (m > 64) {
        throw invalid_argument("This single-word version supports patterns up to 64 characters");
    }

    unordered_map<char, uint64_t> peq;
    for (int i = 0; i < m; i++) {
        peq[pattern[i]] |= (uint64_t(1) << i);
    }

    const uint64_t all_ones = (m == 64) ? ~uint64_t(0) : ((uint64_t(1) << m) - 1);
    const uint64_t top_bit = uint64_t(1) << (m - 1);

    uint64_t pv = all_ones;
    uint64_t mv = 0;
    int score = m;

    int best_score = m;
    int best_pos = 0;

    for (int j = 0; j < n; j++) {
        uint64_t eq = peq.count(text[j]) ? peq[text[j]] : uint64_t(0);

        uint64_t xv = eq | mv;
        uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;

        uint64_t ph = mv | ~(xh | pv);
        uint64_t mh = pv & xh;

        ph &= all_ones;
        mh &= all_ones;

        if (ph & top_bit) score += 1;
        else if (mh & top_bit) score -= 1;

        uint64_t ph_shifted = ((ph << 1) | uint64_t(1)) & all_ones;
        uint64_t mh_shifted = (mh << 1) & all_ones;

        pv = (mh_shifted | ~(xv | ph_shifted)) & all_ones;
        mv = ph_shifted & xv;

        if (score <= best_score) {
            best_score = score;
            best_pos = j + 1;
        }
    }

    return {best_score, best_pos};
}

/** Reference O(len(a) * len(b)) Levenshtein distance, for verification. */
int edit_distance_bruteforce(const string& a, const string& b) {
    int la = a.size(), lb = b.size();
    vector<vector<int>> dp(la + 1, vector<int>(lb + 1, 0));
    for (int i = 0; i <= la; i++) dp[i][0] = i;
    for (int j = 0; j <= lb; j++) dp[0][j] = j;
    for (int i = 1; i <= la; i++) {
        for (int j = 1; j <= lb; j++) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            dp[i][j] = min({dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost});
        }
    }
    return dp[la][lb];
}

int main() {
    cout << "Myers' Bit-Vector Algorithm - Edit Distance" << endl;
    cout << string(60, '=') << endl;

    {
        string pattern = "GATTACA", text = "GATTACA";
        auto result = myers_bit_vector_edit_distance(pattern, text);
        cout << "\nPattern: " << pattern << "\nText:    " << text << endl;
        cout << "Edit distance: " << result.edit_distance
             << " (ends at text[:" << result.end_position << "])" << endl;
    }

    {
        string pattern = "GATTACA", text = "GACTATA";
        auto result = myers_bit_vector_edit_distance(pattern, text);
        int brute = edit_distance_bruteforce(pattern, text);
        cout << "\nPattern: " << pattern << "\nText:    " << text << endl;
        cout << "Edit distance: " << result.edit_distance
             << " (bit-vector) vs " << brute << " (brute force)" << endl;
    }

    {
        string pattern = "ACGT", text = "TTTTACGGTTT";
        auto result = myers_bit_vector_edit_distance(pattern, text);
        cout << "\nSemi-global search: find `pattern` inside a longer `text`" << endl;
        cout << "Pattern: " << pattern << "\nText:    " << text << endl;
        cout << "Best match ends at position " << result.end_position
             << " with edit distance " << result.edit_distance << endl;
    }

    return 0;
}
