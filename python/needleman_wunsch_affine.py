"""
Needleman-Wunsch with Affine Gap Penalties (Gotoh's Algorithm)

Global sequence alignment using affine gap penalties: gap_open + k * gap_extend
for a gap of length k. This is Gotoh's (1982) O(mn)-time reformulation of the
naive O(mn * max(m,n))-time affine-gap extension of Needleman-Wunsch, achieved
by tracking three DP matrices instead of one.

Uses three matrices:
- M[i][j]: optimal global-alignment score of seq1[:i], seq2[:j] ending in a
  match/mismatch
- I[i][j]: optimal score ending with a gap in seq1 (insertion, moves along seq2)
- D[i][j]: optimal score ending with a gap in seq2 (deletion, moves along seq1)

Unlike Smith-Waterman-affine (local alignment, matrices floored at 0 and
traceback stops at 0), this is a *global* alignment: there is no flooring, the
first row/column are initialized with gap-opening costs, and traceback always
runs to (0, 0).

Time Complexity: O(m * n)
Space Complexity: O(m * n) (reducible to O(n) for score-only via Hirschberg-style
recursion, as in the WFA and Hirschberg implementations in this repository)

Reference:
Gotoh, O. (1982). An improved algorithm for matching biological sequences.
Journal of Molecular Biology, 162(3), 705-708.
https://doi.org/10.1016/0022-2836(82)90398-9
"""


def needleman_wunsch_affine(seq1, seq2, match_score=1, mismatch_penalty=-1,
                             gap_open=-3, gap_extend=-1):
    """
    Perform global sequence alignment with affine gap penalties.

    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match_score: Score for matching characters (default: 1)
        mismatch_penalty: Penalty for mismatching characters (default: -1)
        gap_open: Additional penalty charged once when a gap is opened (default: -3)
        gap_extend: Penalty charged per gap character (default: -1)
            A gap of length k costs gap_open + k * gap_extend.

    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)
    """
    m, n = len(seq1), len(seq2)
    NEG_INF = float('-inf')

    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    I = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    D = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    # Initialization: aligning a prefix against nothing but gaps.
    # M[0][j] and M[i][0] stay -inf (j>0 / i>0): a match/mismatch cannot end
    # an alignment that consumed zero characters from the other sequence.
    M[0][0] = 0
    for j in range(1, n + 1):
        I[0][j] = gap_open + j * gap_extend
    for i in range(1, m + 1):
        D[i][0] = gap_open + i * gap_extend

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Gap in seq1 (insertion relative to seq1): move along seq2.
            I[i][j] = max(
                M[i][j - 1] + gap_open + gap_extend,
                I[i][j - 1] + gap_extend,
                D[i][j - 1] + gap_open + gap_extend,
            )
            # Gap in seq2 (deletion relative to seq1): move along seq1.
            D[i][j] = max(
                M[i - 1][j] + gap_open + gap_extend,
                D[i - 1][j] + gap_extend,
                I[i - 1][j] + gap_open + gap_extend,
            )
            s = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            M[i][j] = max(
                M[i - 1][j - 1] + s,
                I[i - 1][j - 1] + s,
                D[i - 1][j - 1] + s,
            )

    # Best score is whichever matrix ends the alignment at (m, n).
    final = max(M[m][n], I[m][n], D[m][n])
    if final == M[m][n]:
        current_matrix = 'M'
    elif final == I[m][n]:
        current_matrix = 'I'
    else:
        current_matrix = 'D'

    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n

    while i > 0 or j > 0:
        if current_matrix == 'M':
            s = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            if M[i][j] == M[i - 1][j - 1] + s:
                current_matrix = 'M'
            elif M[i][j] == I[i - 1][j - 1] + s:
                current_matrix = 'I'
            else:
                current_matrix = 'D'
            i -= 1
            j -= 1
        elif current_matrix == 'I':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            if I[i][j] == I[i][j - 1] + gap_extend:
                current_matrix = 'I'
            elif I[i][j] == M[i][j - 1] + gap_open + gap_extend:
                current_matrix = 'M'
            else:
                current_matrix = 'D'
            j -= 1
        else:  # current_matrix == 'D'
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            if D[i][j] == D[i - 1][j] + gap_extend:
                current_matrix = 'D'
            elif D[i][j] == M[i - 1][j] + gap_open + gap_extend:
                current_matrix = 'M'
            else:
                current_matrix = 'I'
            i -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return ''.join(aligned_seq1), ''.join(aligned_seq2), final


if __name__ == "__main__":
    seq1 = "GGTTGACTA"
    seq2 = "TGTTACGG"

    print("Needleman-Wunsch with Affine Gap Penalties (Gotoh 1982) - Global Alignment")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()

    aligned1, aligned2, score = needleman_wunsch_affine(seq1, seq2, 2, -1, -3, -1)
    print("With affine gaps (gap_open=-3, gap_extend=-1):")
    print(f"Aligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
    print()

    # A run of indels should cost less than the same number of scattered gaps.
    seq3 = "ACGTACGTACGT"
    seq4 = "ACGTAAACGT"
    print("Example 2 - clustered indel:")
    print(f"Sequence 1: {seq3}")
    print(f"Sequence 2: {seq4}")
    aligned3, aligned4, score2 = needleman_wunsch_affine(seq3, seq4, 2, -1, -5, -1)
    print(f"Aligned Sequence 1: {aligned3}")
    print(f"Aligned Sequence 2: {aligned4}")
    print(f"Alignment Score: {score2}")
