"""
Myers' Bit-Vector Algorithm for Approximate String Matching / Edit Distance

Computes the edit (Levenshtein) distance between a pattern and every prefix of
a text - or simply between two whole sequences - by packing an entire column
of the classic O(mn) edit-distance DP matrix into a handful of machine words
and updating it with O(1) bitwise operations per column.

Key idea: instead of storing the DP column of absolute scores, store the
*differences* between vertically adjacent cells as two bit-vectors:
  Pv: positions where the column value increases by 1 going down (+1 deltas)
  Mv: positions where the column value decreases by 1 going down (-1 deltas)
A whole column update - normally m scalar operations - collapses to a fixed
number of AND/OR/XOR/ADD/shift operations on ceil(m / w)-word bit-vectors,
where w is the machine word size (e.g. 64).

Time Complexity: O(n * ceil(m / w)) -- O(n) for patterns up to word size w
Space Complexity: O(m) for the character bitmask table (Peq), O(ceil(m/w))
for the running state

This is the algorithm underlying the `edlib` library (Sosic & Sikic, 2017)
and is used as the extension kernel inside BWA-MEM, GenASM and other
high-throughput short-read aligners.

Reference:
Myers, G. (1999). A fast bit-vector algorithm for approximate string matching
based on dynamic programming. Journal of the ACM, 46(3), 395-415.
https://doi.org/10.1145/316542.316550

See also:
Sosic, M., & Sikic, M. (2017). Edlib: a C/C++ library for fast, exact sequence
alignment using edit distance. Bioinformatics, 33(9), 1394-1395.
https://doi.org/10.1093/bioinformatics/btw753
https://github.com/Martinsos/edlib
"""


def myers_bit_vector_edit_distance(pattern, text):
    """
    Compute the edit distance between `pattern` and `text` using Myers'
    bit-vector algorithm, and the minimal edit distance of `pattern` against
    any prefix of `text` (useful for semi-global / infix matching).

    This reference implementation uses Python's arbitrary-precision integers
    as the bit-vector words, so it works for patterns of any length (a real
    C/C++ implementation instead packs the vectors into fixed-width 64-bit
    words and loops over ceil(m / 64) blocks - see the accompanying .cpp file
    for that word-parallel version).

    Args:
        pattern: The (shorter) sequence to search for, length m
        text: The (longer) sequence to search within, length n

    Returns:
        tuple: (edit_distance, best_end_position)
            edit_distance: edit distance between pattern and the best-matching
                prefix of text (== classic Levenshtein distance if
                best_end_position == len(text))
            best_end_position: index into text (1-based, i.e. text[:pos]) of
                the prefix achieving the minimal distance
    """
    m = len(pattern)
    n = len(text)
    if m == 0:
        return n, n

    # Peq[c] = bitmask with bit i set iff pattern[i] == c
    peq = {}
    for i, c in enumerate(pattern):
        peq[c] = peq.get(c, 0) | (1 << i)

    all_ones = (1 << m) - 1
    top_bit = 1 << (m - 1)

    pv = all_ones   # Pv[i] = 1 => vertical delta at row i is +1
    mv = 0          # Mv[i] = 1 => vertical delta at row i is -1
    score = m       # current edit distance for the column so far (top-to-bottom sum)

    best_score = m
    best_pos = 0

    for j, c in enumerate(text, start=1):
        eq = peq.get(c, 0)

        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq

        ph = mv | ~(xh | pv)
        mh = pv & xh

        # Keep values within the m-bit window (Python ints are unbounded).
        ph &= all_ones
        mh &= all_ones

        if ph & top_bit:
            score += 1
        elif mh & top_bit:
            score -= 1

        ph_shifted = ((ph << 1) | 1) & all_ones
        mh_shifted = (mh << 1) & all_ones

        pv = (mh_shifted | ~(xv | ph_shifted)) & all_ones
        mv = ph_shifted & xv

        if score <= best_score:
            best_score = score
            best_pos = j

    return best_score, best_pos


def edit_distance_bruteforce(a, b):
    """Reference O(len(a) * len(b)) Levenshtein distance, for verification."""
    la, lb = len(a), len(b)
    dp = [[0] * (lb + 1) for _ in range(la + 1)]
    for i in range(la + 1):
        dp[i][0] = i
    for j in range(lb + 1):
        dp[0][j] = j
    for i in range(1, la + 1):
        for j in range(1, lb + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)
    return dp[la][lb]


if __name__ == "__main__":
    print("Myers' Bit-Vector Algorithm - Edit Distance")
    print("=" * 60)

    pattern = "GATTACA"
    text = "GATTACA"
    dist, pos = myers_bit_vector_edit_distance(pattern, text)
    print(f"\nPattern: {pattern}")
    print(f"Text:    {text}")
    print(f"Edit distance: {dist} (ends at text[:{pos}])")

    pattern = "GATTACA"
    text = "GACTATA"
    dist, pos = myers_bit_vector_edit_distance(pattern, text)
    brute = edit_distance_bruteforce(pattern, text)
    print(f"\nPattern: {pattern}")
    print(f"Text:    {text}")
    print(f"Edit distance: {dist} (bit-vector) vs {brute} (brute force)")

    print("\nSemi-global search: find `pattern` inside a longer `text`")
    pattern = "ACGT"
    text = "TTTTACGGTTT"  # 'ACGT' with one mismatch embedded at position 4-7
    dist, pos = myers_bit_vector_edit_distance(pattern, text)
    print(f"Pattern: {pattern}")
    print(f"Text:    {text}")
    print(f"Best match ends at position {pos} with edit distance {dist}")
    print(f"Matched region: '{text[max(0, pos - len(pattern) - dist):pos]}'")
