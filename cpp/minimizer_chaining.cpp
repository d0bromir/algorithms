/**
 * Minimizer Sketching + Co-linear Chaining (Minimap2-style Seeding)
 *
 * Modern long-read aligners (minimap2, HISAT2, GraphAligner, ...) index only
 * a sparse, deterministic subset of k-mers called MINIMIZERS, then find
 * matching (query, reference) seed pairs and CHAIN them - a sparse dynamic
 * program over the seeds, not the underlying bases - to identify collinear
 * runs of seeds forming a candidate alignment. Only the region(s) around
 * good chains are passed to base-level DP, which is what lets these tools
 * scale to whole-genome references.
 *
 * 1. Minimizer sketching: for every window of `w` consecutive k-mers, keep
 *    only the numerically smallest k-mer hash as the window's minimizer.
 * 2. Co-linear chaining: find the highest-scoring subsequence of seeds with
 *    increasing query and reference positions (O(n^2) DP here for clarity;
 *    minimap2 uses a Fenwick-tree-backed O(N log N) version of the same DP).
 *
 * Time Complexity:
 * - Minimizer sketch of a sequence of length L: O(L)
 * - Chaining N seeds: O(N^2) here, O(N log N) in a production chainer
 *
 * Reference:
 * Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
 * Bioinformatics, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191
 * Code: https://github.com/lh3/minimap2
 *
 * Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M., & Yorke, J. A. (2004).
 * Reducing storage requirements for biological sequence comparison.
 * Bioinformatics, 20(18), 3363-3369. https://doi.org/10.1093/bioinformatics/bth408
 */

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>

using namespace std;

static uint64_t kmer_hash(const string& kmer) {
    uint64_t h = 0;
    for (char c : kmer) h = h * 1000003ULL + static_cast<unsigned char>(c);
    return h;
}

struct Minimizer {
    uint64_t hash;
    int pos;
};

/** Compute the minimizer sketch of a sequence: one entry per window-minimum k-mer. */
vector<Minimizer> compute_minimizers(const string& seq, int k, int w) {
    int n = seq.size();
    vector<Minimizer> result;
    if (n < k) return result;

    int num_kmers = n - k + 1;
    vector<uint64_t> kmer_hashes(num_kmers);
    for (int i = 0; i < num_kmers; i++) kmer_hashes[i] = kmer_hash(seq.substr(i, k));

    unordered_map<int, uint64_t> minimizers; // pos -> hash (dedup)
    for (int start = 0; start + w <= num_kmers; start++) {
        int min_idx = 0;
        for (int idx = 1; idx < w; idx++) {
            if (kmer_hashes[start + idx] < kmer_hashes[start + min_idx]) min_idx = idx;
        }
        int pos = start + min_idx;
        minimizers[pos] = kmer_hashes[pos];
    }

    result.reserve(minimizers.size());
    for (auto& [pos, hash] : minimizers) result.push_back({hash, pos});
    sort(result.begin(), result.end(), [](const Minimizer& a, const Minimizer& b) {
        return a.hash != b.hash ? a.hash < b.hash : a.pos < b.pos;
    });
    return result;
}

/** Index a reference sequence by minimizer hash -> list of positions. */
unordered_map<uint64_t, vector<int>> build_minimizer_index(const string& reference, int k, int w) {
    unordered_map<uint64_t, vector<int>> index;
    for (auto& m : compute_minimizers(reference, k, w)) index[m.hash].push_back(m.pos);
    return index;
}

struct Seed {
    int query_pos;
    int ref_pos;
};

/** Find all seed matches between a query's minimizers and a reference index. */
vector<Seed> find_seed_matches(const string& query, const unordered_map<uint64_t, vector<int>>& ref_index,
                                int k, int w) {
    vector<Seed> seeds;
    for (auto& m : compute_minimizers(query, k, w)) {
        auto it = ref_index.find(m.hash);
        if (it == ref_index.end()) continue;
        for (int rpos : it->second) seeds.push_back({m.pos, rpos});
    }
    sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
        return a.query_pos != b.query_pos ? a.query_pos < b.query_pos : a.ref_pos < b.ref_pos;
    });
    seeds.erase(unique(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
        return a.query_pos == b.query_pos && a.ref_pos == b.ref_pos;
    }), seeds.end());
    return seeds;
}

struct ChainResult {
    vector<Seed> chain;
    double score;
};

/** Highest-scoring co-linear chain of seeds via O(n^2) sparse DP. */
ChainResult chain_seeds(vector<Seed> seeds, int k, int max_gap = 50, double gap_penalty_weight = 0.5) {
    if (seeds.empty()) return {{}, 0.0};
    int n = seeds.size();

    vector<double> dp(n, static_cast<double>(k));
    vector<int> parent(n, -1);

    for (int i = 0; i < n; i++) {
        int qi = seeds[i].query_pos, ri = seeds[i].ref_pos;
        int best_j = -1;
        double best_val = k;
        for (int j = 0; j < i; j++) {
            int qj = seeds[j].query_pos, rj = seeds[j].ref_pos;
            if (qj >= qi || rj >= ri) continue;
            int gap_q = qi - (qj + k);
            int gap_r = ri - (rj + k);
            if (gap_q > max_gap || gap_r > max_gap) continue;
            double gap_cost = gap_penalty_weight * abs(gap_q - gap_r);
            double candidate = dp[j] + k - gap_cost;
            if (candidate > best_val) {
                best_val = candidate;
                best_j = j;
            }
        }
        dp[i] = best_val;
        parent[i] = best_j;
    }

    int best_i = 0;
    for (int i = 1; i < n; i++) if (dp[i] > dp[best_i]) best_i = i;

    vector<Seed> chain;
    for (int i = best_i; i != -1; i = parent[i]) chain.push_back(seeds[i]);
    reverse(chain.begin(), chain.end());

    return {chain, dp[best_i]};
}

int main() {
    cout << "Minimizer Sketching + Co-linear Chaining (Minimap2-style)" << endl;
    cout << string(60, '=') << endl;

    string reference = "ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                        "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG";
    string query = "GACCATGGGACCATTGACCTGAGGGTACCGGATT"; // exact substring of reference

    int k = 8, w = 4;
    auto ref_index = build_minimizer_index(reference, k, w);
    auto seeds = find_seed_matches(query, ref_index, k, w);

    cout << "\nReference length: " << reference.size() << endl;
    cout << "Query length: " << query.size() << endl;
    cout << "Reference minimizers indexed: " << ref_index.size() << endl;
    cout << "Seed matches found: " << seeds.size() << endl;

    auto result = chain_seeds(seeds, k);
    cout << "\nBest chain has " << result.chain.size() << " seeds, score=" << result.score << endl;
    if (!result.chain.empty()) {
        int q0 = result.chain.front().query_pos, r0 = result.chain.front().ref_pos;
        int q1 = result.chain.back().query_pos, r1 = result.chain.back().ref_pos;
        cout << "Chain spans query[" << q0 << ":" << q1 + k << "] -> reference[" << r0 << ":" << r1 + k << "]" << endl;
        size_t expected = reference.find(query);
        cout << "(sanity check: true start of query in reference = " << expected << ")" << endl;
    }

    return 0;
}
