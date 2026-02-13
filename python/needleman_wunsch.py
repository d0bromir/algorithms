"""
Needleman-Wunsch Algorithm for Global Sequence Alignment

This is a dynamic programming algorithm for aligning two sequences globally.
It finds the optimal alignment between two sequences by maximizing the alignment score.

Time Complexity: O(m * n) where m and n are the lengths of the sequences
Space Complexity: O(m * n)
"""


def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """
    Perform global sequence alignment using Needleman-Wunsch algorithm.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match_score: Score for matching characters (default: 1)
        mismatch_penalty: Penalty for mismatching characters (default: -1)
        gap_penalty: Penalty for gaps (default: -1)
    
    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize DP matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column with gap penalties
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match = dp[i - 1][j - 1] + match_score
            else:
                match = dp[i - 1][j - 1] + mismatch_penalty
            
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            
            dp[i][j] = max(match, delete, insert)
    
    # Traceback to find alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            current_score = dp[i][j]
            diagonal_score = dp[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            
            if current_score == diagonal_score:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and current_score == dp[i - 1][j] + gap_penalty:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1
        elif i > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
    
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    return ''.join(aligned_seq1), ''.join(aligned_seq2), dp[m][n]


if __name__ == "__main__":
    # Example usage
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    
    aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
    
    print("Needleman-Wunsch Global Alignment")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
