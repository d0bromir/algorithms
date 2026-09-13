"""
Strobemer Seeding + MinHash Sequence-Identity Estimation

Two related fast-mapping ideas from the newest generation of aligners, shown
together because they solve complementary problems:

1. STROBEMERS (Sahlin, 2021; used in strobealign, Sahlin 2022): fixed-length
   k-mers are brittle under indels - a single insertion shifts every
   downstream k-mer, destroying exact matches. A strobemer instead links
   together several short "strobes" (mini k-mers) chosen from successive
   windows *by a hash-minimizing rule*, so the same strobemer tends to
   reappear even when the sequence between strobes is edited, giving seeds
   that survive substitutions and indels much better than plain k-mers of
   the same total span. This module implements "randstrobes": a strobemer
   is (strobe_1, strobe_2), where strobe_2 is the position in a downstream
   window whose k-mer hash, when combined with strobe_1's hash, is minimal.

2. MinHash / Jaccard identity estimation (Jain et al., 2018, MashMap; based
   on Broder's MinHash and Ondov et al.'s Mash): rather than aligning two
   sequences to know how similar they are, estimate the Jaccard similarity
   of their k-mer sets from a small, fixed-size random sample (the
   "sketch") of each set's minimum hash values. This turns an O(nm)
   alignment into an O(sketch_size) set comparison, which is what lets
   MashMap/Mash cluster or pre-filter millions of long reads against huge
   reference databases before any base-level alignment is attempted.

Time Complexity:
- Strobemer sketch of a sequence of length L: O(L)
- MinHash sketch of a sequence of length L: O(L log(sketch_size)) (heap-based)
- Jaccard/identity estimate from two sketches of size s: O(s)

Reference:
Sahlin, K. (2022). Effective sequence similarity detection with strobemers.
Genome Research, 31(11), 2080-2094. https://doi.org/10.1101/gr.275648.121
Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and
accurate read alignment. Genome Biology, 23, 260.
https://doi.org/10.1186/s13059-022-02831-7
Code: https://github.com/ksahlin/strobealign

Jain, C., Dilthey, A., Koren, S., Aluru, S., & Phillippy, A. M. (2018). A
fast approximate algorithm for mapping long reads to large reference
databases. Journal of Computational Biology, 25(7), 766-779.
https://doi.org/10.1089/cmb.2018.0036
Code: https://github.com/marbl/MashMap
"""

import heapq


def _hash(s):
    h = 0
    for c in s:
        h = (h * 1000003 + ord(c)) & 0xFFFFFFFFFFFFFFFF
    return h


def randstrobes(seq, strobe_len=6, w_min=6, w_max=20):
    """
    Generate randstrobes: order-2 strobemers linking two short strobes.

    For each position i, the first strobe is seq[i:i+strobe_len]. The second
    strobe is chosen from the window [i + strobe_len + w_min,
    i + strobe_len + w_max] as whichever candidate k-mer minimizes the
    combined hash hash(strobe1) ^ hash(candidate) - i.e., a link function
    that depends on both strobes, so the choice adapts to content rather
    than being fixed by position alone (unlike plain spaced seeds).

    Args:
        seq: input sequence
        strobe_len: length of each strobe (mini k-mer)
        w_min: minimum gap between the end of strobe 1 and the start of the
            window strobe 2 is drawn from
        w_max: maximum offset defining the window's far edge

    Returns:
        list of (combined_hash, i, j) tuples: i = start of strobe 1,
        j = start of strobe 2
    """
    n = len(seq)
    strobemers = []
    for i in range(0, n - strobe_len + 1):
        s1 = seq[i:i + strobe_len]
        h1 = _hash(s1)

        window_start = i + strobe_len + w_min
        window_end = min(i + strobe_len + w_max, n - strobe_len + 1)
        if window_start >= window_end:
            continue

        best_j, best_combined = None, None
        for j in range(window_start, window_end):
            s2 = seq[j:j + strobe_len]
            combined = h1 ^ _hash(s2)
            if best_combined is None or combined < best_combined:
                best_combined = combined
                best_j = j

        if best_j is not None:
            strobemers.append((best_combined, i, best_j))

    return strobemers


def strobemer_index(reference, strobe_len=6, w_min=6, w_max=20):
    """Index a reference by strobemer hash -> list of strobe-1 positions."""
    index = {}
    for h, i, _j in randstrobes(reference, strobe_len, w_min, w_max):
        index.setdefault(h, []).append(i)
    return index


def find_strobemer_matches(query, ref_index, strobe_len=6, w_min=6, w_max=20):
    """Find (query_pos, ref_pos) matches between query strobemers and an index."""
    matches = []
    for h, i, _j in randstrobes(query, strobe_len, w_min, w_max):
        for rpos in ref_index.get(h, []):
            matches.append((i, rpos))
    return sorted(matches)


def minhash_sketch(seq, k=12, sketch_size=64):
    """
    Compute a bottom-k MinHash sketch: the `sketch_size` smallest distinct
    k-mer hash values observed in `seq`.

    Args:
        seq: input sequence
        k: k-mer length
        sketch_size: number of minimum hashes to retain

    Returns:
        sorted list of the smallest `sketch_size` k-mer hash values (a
        deterministic, fixed-size summary of the sequence's k-mer set)
    """
    n = len(seq)
    if n < k:
        return []

    seen = set()
    heap = []  # max-heap via negation, capped at sketch_size
    for i in range(n - k + 1):
        h = _hash(seq[i:i + k])
        if h in seen:
            continue
        seen.add(h)
        if len(heap) < sketch_size:
            heapq.heappush(heap, -h)
        elif -heap[0] > h:
            heapq.heapreplace(heap, -h)

    return sorted(-x for x in heap)


def estimate_jaccard(sketch_a, sketch_b, sketch_size=None):
    """
    Estimate the Jaccard similarity of two k-mer sets from their bottom-k
    MinHash sketches, following the Mash/MashMap estimator: merge the two
    sketches, keep the smallest `sketch_size` values overall, and the
    fraction of those that appear in BOTH original sketches approximates
    the true Jaccard index J = |A ∩ B| / |A ∪ B|.

    Returns:
        float in [0, 1]: estimated Jaccard similarity
    """
    if not sketch_a or not sketch_b:
        return 0.0
    if sketch_size is None:
        sketch_size = max(len(sketch_a), len(sketch_b))

    set_a, set_b = set(sketch_a), set(sketch_b)
    merged = sorted(set_a | set_b)[:sketch_size]
    if not merged:
        return 0.0
    shared = sum(1 for h in merged if h in set_a and h in set_b)
    return shared / len(merged)


def mash_distance(jaccard, k):
    """
    Convert a Jaccard estimate to the Mash "evolutionary distance" (an
    estimate of per-base substitution rate), following Ondov et al. (2016):
    D = -1/k * ln(2J / (1+J)).
    """
    import math
    if jaccard <= 0:
        return 1.0
    if jaccard >= 1:
        return 0.0
    return -1.0 / k * math.log(2 * jaccard / (1 + jaccard))


if __name__ == "__main__":
    print("Strobemer Seeding + MinHash Identity Estimation")
    print("=" * 60)

    reference = ("ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                 "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG"
                 "TTCCGGAATTCCGGATCGATCGGATCCAAGCTTGGATCCACTAGTCCAGT")
    query = reference[60:130]  # substring, will test seeding

    print("\n--- Strobemer seeding ---")
    idx = strobemer_index(reference, strobe_len=5, w_min=4, w_max=12)
    matches = find_strobemer_matches(query, idx, strobe_len=5, w_min=4, w_max=12)
    print(f"Reference length: {len(reference)}, query length: {len(query)}")
    print(f"Strobemer matches found: {len(matches)}")
    if matches:
        diffs = [r - q for q, r in matches]
        # Most matches should agree on the same (ref_pos - query_pos) offset.
        from collections import Counter
        common_offset, count = Counter(diffs).most_common(1)[0]
        print(f"Most common offset: {common_offset} (true offset: 60), "
              f"supported by {count}/{len(matches)} matches")

    print("\n--- MinHash Jaccard / identity estimation ---")
    seq_a = reference
    seq_b = list(reference)
    import random
    random.seed(0)
    for _ in range(15):  # point mutations scattered across the sequence
        pos = random.randrange(len(seq_b))
        seq_b[pos] = random.choice('ACGT')
    seq_b = ''.join(seq_b)

    k = 12
    sketch_a = minhash_sketch(seq_a, k=k, sketch_size=200)
    sketch_b = minhash_sketch(seq_b, k=k, sketch_size=200)
    jaccard = estimate_jaccard(sketch_a, sketch_b)
    distance = mash_distance(jaccard, k)
    print(f"Sequences differ by ~{15}/{len(reference)} = {15/len(reference):.1%} point mutations")
    print(f"Estimated Jaccard similarity: {jaccard:.3f}")
    print(f"Estimated Mash distance (per-base substitution rate): {distance:.3f}")
