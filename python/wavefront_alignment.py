"""
Wavefront Alignment Algorithm (WFA) - Gap-Affine Exact Alignment

WFA is an exact, provably-optimal gap-affine global-alignment algorithm that
runs in O(n * s + s^2) time and O(s^2) memory, where n = max(len(seq1),
len(seq2)) and s is the optimal alignment score (edit-distance-like penalty
score, or "distance"). For similar sequences, s is small, so WFA is
near-linear in practice - a fundamentally different complexity regime from
the O(mn) of Needleman-Wunsch/Gotoh, which is oblivious to how similar the
sequences actually are.

Key idea: instead of filling the DP matrix cell by cell (indexed by
sequence positions i, j), WFA is indexed by score. For each candidate score
s = 0, 1, 2, ..., it tracks, for every diagonal k = i - j, the furthest-reaching
point (i, j) reachable with exactly that score - a "wavefront". Along a
diagonal, matching characters are free ("greedy extension"), so a wavefront
snaps forward through exact matches at no cost. The wavefront for score s is
computed purely from the wavefronts of smaller scores (s - mismatch_penalty,
s - gap_open - gap_extend, s - gap_extend), so the algorithm terminates the
moment a wavefront reaches the bottom-right corner - it never explores
alignments worse than the optimum.

This module implements WFA for the gap-affine model with a MINIMIZATION
convention (mismatch/gap penalties are positive costs, lower is better),
which is the convention used in the original paper and its WFA2-lib
reference implementation.

Time Complexity: O(n * s + s^2) -- near-linear for similar sequences (small s)
Space Complexity: O(s^2) for this straightforward implementation (the
    "banded"/"BiWFA" variants in WFA2-lib reduce this further to O(s) by
    divide-and-conquer, mirroring how Hirschberg reduces Needleman-Wunsch
    space in this repository)

Reference:
Marco-Sola, S., Moure, J. C., Moreto, M., & Espinosa, A. (2021). Fast
gap-affine pairwise alignment using the wavefront algorithm. Bioinformatics,
37(4), 456-463. https://doi.org/10.1093/bioinformatics/btaa777
Code: https://github.com/smarco/WFA2-lib
"""


def _extend_wavefront(wf_diag, seq1, seq2, offset):
    """Greedily extend along a diagonal through exact matches (free operations)."""
    m, n = len(seq1), len(seq2)
    i = offset
    j = offset + wf_diag
    while i < m and j < n and seq1[i] == seq2[j]:
        i += 1
        j += 1
    return i


def wavefront_alignment(seq1, seq2, mismatch=4, gap_open=6, gap_extend=2,
                         max_score=None):
    """
    Compute the optimal gap-affine alignment score (and end diagonal/offset
    trace) between seq1 and seq2 using the Wavefront Alignment algorithm.

    Costs are non-negative penalties (0 = match, this function MINIMIZES total
    cost), matching the convention of the original WFA paper and WFA2-lib.

    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        mismatch: Cost of a substitution (default: 4)
        gap_open: Cost of opening a new gap, charged once (default: 6)
        gap_extend: Cost per gap character (default: 2)
        max_score: Optional cap on score to search (defaults to a safe upper
            bound derived from sequence lengths)

    Returns:
        int: optimal alignment score (edit cost) under the gap-affine model
    """
    m, n = len(seq1), len(seq2)
    if max_score is None:
        max_score = (m + n) * max(mismatch, gap_open + gap_extend) + 1

    target_diag = n - m  # diagonal containing the bottom-right corner (m, n)
    # Only diagonals k in [-m, n] can ever lie on a valid alignment path
    # between (0, 0) and (m, n); anything outside that range has already
    # overshot one sequence and must not be propagated further.
    min_diag, max_diag = -m, n

    def in_bounds(k, i):
        return min_diag <= k <= max_diag and 0 <= i <= m and 0 <= i + k <= n

    # Each wavefront is a dict: diagonal -> furthest-reached seq1-offset `i`
    # (the seq2-offset j = i + diagonal is implied).
    M = {}  # ends in a match/mismatch
    I = {}  # ends in an insertion (gap in seq1, i.e. consumed a seq2 char)
    D = {}  # ends in a deletion (gap in seq2, i.e. consumed a seq1 char)

    M[0] = {0: _extend_wavefront(0, seq1, seq2, 0)}

    if M[0][0] == m and n == m:
        return 0

    for s in range(1, max_score + 1):
        m_s, i_s, d_s = {}, {}, {}

        # --- Insertion wavefront: came from M or I at score s - (gap cost) ---
        # An insertion consumes a seq2 character (j += 1, i unchanged); since
        # j = offset + diag, that means diag increases by 1.
        for src, cost in ((M, gap_open + gap_extend), (I, gap_extend)):
            prev = src.get(s - cost)
            if prev is None:
                continue
            for k, i in prev.items():
                nk = k + 1
                if in_bounds(nk, i) and i > i_s.get(nk, -1):
                    i_s[nk] = i

        # --- Deletion wavefront: came from M or D at score s - (gap cost) ---
        # A deletion consumes a seq1 character (i += 1, j unchanged); since
        # j = offset + diag, that means diag decreases by 1.
        for src, cost in ((M, gap_open + gap_extend), (D, gap_extend)):
            prev = src.get(s - cost)
            if prev is None:
                continue
            for k, i in prev.items():
                nk = k - 1
                cand = i + 1
                if in_bounds(nk, cand) and cand > d_s.get(nk, -1):
                    d_s[nk] = cand

        # --- Match/mismatch wavefront ---
        prev_mismatch = M.get(s - mismatch, {})
        diagonals = set(prev_mismatch) | set(i_s) | set(d_s)
        for k in diagonals:
            best = -1
            if k in prev_mismatch and in_bounds(k, prev_mismatch[k] + 1):
                best = max(best, prev_mismatch[k] + 1)  # substitution, then re-extend below
            if k in i_s:
                best = max(best, i_s[k])
            if k in d_s:
                best = max(best, d_s[k])
            if best >= 0:
                m_s[k] = _extend_wavefront(k, seq1, seq2, best)

        if i_s:
            I[s] = i_s
        if d_s:
            D[s] = d_s
        if m_s:
            M[s] = m_s

        if target_diag in m_s and m_s[target_diag] == m:
            j = m_s[target_diag] + target_diag
            if j == n:
                return s

    raise RuntimeError("max_score exceeded; sequences may be too dissimilar")


if __name__ == "__main__":
    print("Wavefront Alignment Algorithm (WFA) - Marco-Sola et al. 2021")
    print("=" * 60)

    seq1 = "GATTACA"
    seq2 = "GATCACA"
    score = wavefront_alignment(seq1, seq2)
    print(f"\nSequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Optimal gap-affine edit score: {score} (single substitution)")

    seq1 = "ACGTACGTACGT"
    seq2 = "ACGTAAACGT"
    score = wavefront_alignment(seq1, seq2)
    print(f"\nSequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Optimal gap-affine edit score: {score}")

    print("\nWhy WFA matters: for near-identical long reads (large n, small s),")
    print("its O(n*s + s^2) cost is close to linear, vs. O(n^2) for classic DP.")
