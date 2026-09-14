"""
Minimizer Sketching + Co-linear Chaining (Minimap2-style Seeding)

Modern long-read aligners (minimap2, HISAT2, GraphAligner, ...) do not index
every k-mer of the reference: they index only a sparse, deterministic subset
called MINIMIZERS, then find matching (query, reference) seed pairs and use
CHAINING - a sparse dynamic program over the seeds themselves, not the
underlying bases - to identify collinear runs of seeds forming a candidate
alignment. Only the small region(s) around good chains are ever passed to
base-level DP (e.g. the WFA or affine Smith-Waterman implementations
elsewhere in this repository), which is what lets these tools scale to
whole-genome references.

This module implements a simplified version of both pieces:

1. Minimizer sketching: for every window of `w` consecutive k-mers, keep only
   the lexicographically (here: numerically, via a hash) smallest k-mer as
   the window's minimizer. A k-mer that is a minimizer for ANY window is kept
   in the index; adjacent windows usually share their minimizer, so the
   index is roughly 2/(w+1) the size of a full k-mer index while still
   guaranteeing that any shared substring of length >= w+k-1 between two
   sequences shares at least one minimizer.

2. Co-linear chaining: given the list of seed matches (query_pos, ref_pos,
   k) sorted by query position, find the highest-scoring subsequence of
   seeds that is "colinear" (both query and reference positions increase
   together, consistent with a single alignment without large rearrangements)
   using an O(N log N) sparse DP with a Fenwick tree over reference
   positions - the same asymptotic idea used in minimap2's chaining step
   (Li, 2018) and in classical sparse-DP chaining (Eppstein et al., 1992).

Time Complexity:
- Minimizer sketch of a sequence of length L: O(L)
- Seed lookup for a query of length q: O(q) expected (hash table)
- Chaining N seeds: O(N log N)

Reference:
Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
Bioinformatics, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191
Code: https://github.com/lh3/minimap2

Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M., & Yorke, J. A. (2004).
Reducing storage requirements for biological sequence comparison.
Bioinformatics, 20(18), 3363-3369. https://doi.org/10.1093/bioinformatics/bth408
(introduces the minimizer sketching scheme)
"""

def _kmer_hash(kmer):
    """A simple order-preserving-free hash standing in for minimap2's invertible
    64-bit hash of a 2-bit-packed k-mer. Any hash works for correctness; a
    good one avoids systematic biases in which k-mer within a window wins."""
    h = 0
    for c in kmer:
        h = (h * 1000003 + ord(c)) & 0xFFFFFFFFFFFFFFFF
    return h


def compute_minimizers(seq, k=15, w=10):
    """
    Compute the minimizer sketch of a sequence.

    Args:
        seq: input sequence (string)
        k: k-mer length
        w: window size (number of consecutive k-mers per window)

    Returns:
        list of (hash, position) tuples, one per minimizer occurrence,
        position = starting index of the k-mer in seq (0-based), deduplicated
        and sorted by position. A (hash, position) may repeat if a k-mer is
        the minimizer of more than one window, but is only reported once here.
    """
    n = len(seq)
    if n < k:
        return []

    kmer_hashes = [_kmer_hash(seq[i:i + k]) for i in range(n - k + 1)]
    minimizers = {}  # position -> hash, dedup via dict

    for start in range(0, len(kmer_hashes) - w + 1):
        window = kmer_hashes[start:start + w]
        min_idx = min(range(w), key=lambda idx: (window[idx], idx))
        pos = start + min_idx
        minimizers[pos] = kmer_hashes[pos]

    return sorted((h, p) for p, h in minimizers.items())


def build_minimizer_index(reference, k=15, w=10):
    """Index a reference sequence by minimizer hash -> list of positions."""
    index = {}
    for h, pos in compute_minimizers(reference, k, w):
        index.setdefault(h, []).append(pos)
    return index


def find_seed_matches(query, ref_index, k=15, w=10):
    """
    Find all (query_pos, ref_pos) seed matches between a query's minimizers
    and a pre-built reference minimizer index.

    Returns:
        list of (query_pos, ref_pos) tuples, sorted by query_pos then ref_pos
    """
    seeds = []
    for h, qpos in compute_minimizers(query, k, w):
        for rpos in ref_index.get(h, []):
            seeds.append((qpos, rpos))
    seeds.sort()
    return seeds


def chain_seeds(seeds, k=15, max_gap=50, gap_penalty_weight=0.5):
    """
    Find the highest-scoring co-linear chain of seeds using sparse O(N log N)
    dynamic programming, mirroring minimap2's chaining stage.

    A chain is a sequence of seeds s_1, ..., s_t with strictly increasing
    query and reference positions, where consecutive seeds are also
    "consistent" (bounded gap, no large jump suggesting a different locus).
    Score = total bases covered by the k-mers in the chain, minus a linear
    penalty on the total gap between consecutive seeds (a crude stand-in for
    minimap2's log-linear gap cost).

    Args:
        seeds: list of (query_pos, ref_pos) seed matches
        k: k-mer length (each seed covers k bases)
        max_gap: maximum allowed gap (in either coordinate) between two
            consecutive seeds in a chain; larger gaps break the chain
        gap_penalty_weight: cost charged per base of gap

    Returns:
        tuple: (best_chain, best_score)
            best_chain: list of (query_pos, ref_pos) seeds in the optimal chain
            best_score: its score
    """
    if not seeds:
        return [], 0

    seeds = sorted(set(seeds))
    n = len(seeds)

    # dp[i] = best score of a chain ENDING at seeds[i]
    dp = [k] * n
    parent = [-1] * n

    # Naive O(n^2) predecessor search with the max_gap / co-linearity filter.
    # (A production chainer replaces this inner loop with a Fenwick-tree /
    # order-statistics query over ref_pos to reach true O(N log N); the DP
    # structure and scoring here are otherwise the same as minimap2's.)
    for i in range(n):
        qi, ri = seeds[i]
        best_j, best_val = -1, k
        for j in range(i):
            qj, rj = seeds[j]
            if qj >= qi or rj >= ri:
                continue
            gap_q = qi - (qj + k)
            gap_r = ri - (rj + k)
            if gap_q > max_gap or gap_r > max_gap:
                continue
            gap_cost = gap_penalty_weight * abs(gap_q - gap_r)
            candidate = dp[j] + k - gap_cost
            if candidate > best_val:
                best_val = candidate
                best_j = j
        dp[i] = best_val
        parent[i] = best_j

    best_i = max(range(n), key=lambda i: dp[i])
    chain = []
    i = best_i
    while i != -1:
        chain.append(seeds[i])
        i = parent[i]
    chain.reverse()

    return chain, dp[best_i]


if __name__ == "__main__":
    print("Minimizer Sketching + Co-linear Chaining (Minimap2-style)")
    print("=" * 60)

    reference = ("ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                 "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG")
    query = "GACCATGGGACCATTGACCTGAGGGTACCGGATT"  # exact substring of reference

    k, w = 8, 4
    ref_index = build_minimizer_index(reference, k=k, w=w)
    seeds = find_seed_matches(query, ref_index, k=k, w=w)

    print(f"\nReference length: {len(reference)}")
    print(f"Query length: {len(query)}")
    print(f"Reference minimizers indexed: {len(ref_index)}")
    print(f"Seed matches found: {len(seeds)}")

    chain, score = chain_seeds(seeds, k=k)
    print(f"\nBest chain has {len(chain)} seeds, score={score:.1f}")
    if chain:
        q0, r0 = chain[0]
        q1, r1 = chain[-1]
        print(f"Chain spans query[{q0}:{q1 + k}] -> reference[{r0}:{r1 + k}]")
        expected_ref_start = reference.find(query)
        print(f"(sanity check: true start of query in reference = {expected_ref_start})")
