/**
 * Strobemer Seeding + MinHash Sequence-Identity Estimation
 *
 * Two related fast-mapping ideas from the newest generation of aligners:
 *
 * 1. STROBEMERS (Sahlin, 2021; used in strobealign): a strobemer links
 *    together short "strobes" chosen from successive windows by a
 *    content-dependent hash-minimizing rule, so the same strobemer tends to
 *    reappear even across indels, unlike fixed-offset k-mers. This version
 *    implements order-2 "randstrobes".
 *
 * 2. MinHash / Jaccard identity estimation (Jain et al., 2018, MashMap;
 *    building on Broder's MinHash and Ondov et al.'s Mash): estimate the
 *    Jaccard similarity of two k-mer sets from a small, fixed-size sample
 *    (sketch) of each set's minimum hash values, turning an O(nm) alignment
 *    into an O(sketch_size) set comparison.
 *
 * Time Complexity:
 * - Strobemer sketch of a sequence of length L: O(L * (w_max - w_min))
 * - MinHash sketch of a sequence of length L: O(L log(sketch_size))
 * - Jaccard/identity estimate from two sketches of size s: O(s log s)
 *
 * Reference:
 * Sahlin, K. (2022). Effective sequence similarity detection with strobemers.
 * Genome Research, 31(11), 2080-2094. https://doi.org/10.1101/gr.275648.121
 * Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and
 * accurate read alignment. Genome Biology, 23, 260.
 * https://doi.org/10.1186/s13059-022-02831-7
 * Code: https://github.com/ksahlin/strobealign
 *
 * Jain, C., Dilthey, A., Koren, S., Aluru, S., & Phillippy, A. M. (2018). A
 * fast approximate algorithm for mapping long reads to large reference
 * databases. Journal of Computational Biology, 25(7), 766-779.
 * https://doi.org/10.1089/cmb.2018.0036
 * Code: https://github.com/marbl/MashMap
 */

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <map>

using namespace std;

static uint64_t str_hash(const string& s) {
    uint64_t h = 0;
    for (char c : s) h = h * 1000003ULL + static_cast<unsigned char>(c);
    return h;
}

struct Strobemer {
    uint64_t hash;
    int i, j;
};

/** Generate order-2 randstrobes linking strobe 1 at i to strobe 2 in [i+strobe_len+w_min, ...+w_max). */
vector<Strobemer> randstrobes(const string& seq, int strobe_len = 6, int w_min = 6, int w_max = 20) {
    int n = seq.size();
    vector<Strobemer> result;
    for (int i = 0; i + strobe_len <= n; i++) {
        uint64_t h1 = str_hash(seq.substr(i, strobe_len));

        int window_start = i + strobe_len + w_min;
        int window_end = min(i + strobe_len + w_max, n - strobe_len + 1);
        if (window_start >= window_end) continue;

        int best_j = -1;
        uint64_t best_combined = 0;
        bool have_best = false;
        for (int j = window_start; j < window_end; j++) {
            uint64_t h2 = str_hash(seq.substr(j, strobe_len));
            uint64_t combined = h1 ^ h2;
            if (!have_best || combined < best_combined) {
                have_best = true;
                best_combined = combined;
                best_j = j;
            }
        }
        if (have_best) result.push_back({best_combined, i, best_j});
    }
    return result;
}

unordered_map<uint64_t, vector<int>> strobemer_index(const string& reference, int strobe_len = 6,
                                                       int w_min = 6, int w_max = 20) {
    unordered_map<uint64_t, vector<int>> index;
    for (auto& sm : randstrobes(reference, strobe_len, w_min, w_max)) index[sm.hash].push_back(sm.i);
    return index;
}

vector<pair<int, int>> find_strobemer_matches(const string& query,
                                               const unordered_map<uint64_t, vector<int>>& ref_index,
                                               int strobe_len = 6, int w_min = 6, int w_max = 20) {
    vector<pair<int, int>> matches;
    for (auto& sm : randstrobes(query, strobe_len, w_min, w_max)) {
        auto it = ref_index.find(sm.hash);
        if (it == ref_index.end()) continue;
        for (int rpos : it->second) matches.push_back({sm.i, rpos});
    }
    sort(matches.begin(), matches.end());
    return matches;
}

/** Bottom-k MinHash sketch: the `sketch_size` smallest distinct k-mer hashes. */
vector<uint64_t> minhash_sketch(const string& seq, int k = 12, size_t sketch_size = 64) {
    int n = seq.size();
    set<uint64_t> seen;
    set<uint64_t, greater<uint64_t>> top; // acts as a bounded max-first set

    for (int i = 0; i + k <= n; i++) {
        uint64_t h = str_hash(seq.substr(i, k));
        if (seen.count(h)) continue;
        seen.insert(h);
        if (top.size() < sketch_size) {
            top.insert(h);
        } else if (h < *top.begin()) {
            top.erase(top.begin());
            top.insert(h);
        }
    }

    vector<uint64_t> result(top.begin(), top.end());
    sort(result.begin(), result.end());
    return result;
}

/** Estimate Jaccard similarity from two bottom-k MinHash sketches (Mash estimator). */
double estimate_jaccard(const vector<uint64_t>& sketch_a, const vector<uint64_t>& sketch_b,
                         size_t sketch_size = 0) {
    if (sketch_a.empty() || sketch_b.empty()) return 0.0;
    if (sketch_size == 0) sketch_size = max(sketch_a.size(), sketch_b.size());

    set<uint64_t> set_a(sketch_a.begin(), sketch_a.end());
    set<uint64_t> set_b(sketch_b.begin(), sketch_b.end());
    set<uint64_t> merged;
    merged.insert(set_a.begin(), set_a.end());
    merged.insert(set_b.begin(), set_b.end());

    size_t count = 0, shared = 0;
    for (uint64_t h : merged) {
        if (count >= sketch_size) break;
        count++;
        if (set_a.count(h) && set_b.count(h)) shared++;
    }
    if (count == 0) return 0.0;
    return static_cast<double>(shared) / static_cast<double>(count);
}

/** Convert a Jaccard estimate to the Mash evolutionary distance (Ondov et al., 2016). */
double mash_distance(double jaccard, int k) {
    if (jaccard <= 0) return 1.0;
    if (jaccard >= 1) return 0.0;
    return -1.0 / k * log(2 * jaccard / (1 + jaccard));
}

int main() {
    cout << "Strobemer Seeding + MinHash Identity Estimation" << endl;
    cout << string(60, '=') << endl;

    string reference = "ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                        "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG"
                        "TTCCGGAATTCCGGATCGATCGGATCCAAGCTTGGATCCACTAGTCCAGT";
    string query = reference.substr(60, 70);

    cout << "\n--- Strobemer seeding ---" << endl;
    auto idx = strobemer_index(reference, 5, 4, 12);
    auto matches = find_strobemer_matches(query, idx, 5, 4, 12);
    cout << "Reference length: " << reference.size() << ", query length: " << query.size() << endl;
    cout << "Strobemer matches found: " << matches.size() << endl;

    if (!matches.empty()) {
        map<int, int> offset_counts;
        for (auto& [q, r] : matches) offset_counts[r - q]++;
        auto best = max_element(offset_counts.begin(), offset_counts.end(),
                                 [](auto& a, auto& b) { return a.second < b.second; });
        cout << "Most common offset: " << best->first << " (true offset: 60), supported by "
             << best->second << "/" << matches.size() << " matches" << endl;
    }

    cout << "\n--- MinHash Jaccard / identity estimation ---" << endl;
    string seq_a = reference;
    string seq_b = reference;
    srand(0);
    const char bases[] = "ACGT";
    for (int m = 0; m < 15; m++) {
        int pos = rand() % seq_b.size();
        seq_b[pos] = bases[rand() % 4];
    }

    int k = 12;
    auto sketch_a = minhash_sketch(seq_a, k, 200);
    auto sketch_b = minhash_sketch(seq_b, k, 200);
    double jaccard = estimate_jaccard(sketch_a, sketch_b);
    double distance = mash_distance(jaccard, k);
    cout << "Estimated Jaccard similarity: " << jaccard << endl;
    cout << "Estimated Mash distance (per-base substitution rate): " << distance << endl;

    return 0;
}
