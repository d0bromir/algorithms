"""
Smith-Waterman Algorithm for Local Sequence Alignment

This is a dynamic programming algorithm for performing local sequence alignment.
It finds the optimal local alignment between two sequences by maximizing the alignment score.

Time Complexity: O(m * n) where m and n are the lengths of the sequences
Space Complexity: O(m * n)
"""


def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    """
    Perform local sequence alignment using Smith-Waterman algorithm.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match_score: Score for matching characters (default: 2)
        mismatch_penalty: Penalty for mismatching characters (default: -1)
        gap_penalty: Penalty for gaps (default: -1)
    
    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize DP matrix with zeros (key difference from Needleman-Wunsch)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = (0, 0)
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match = dp[i - 1][j - 1] + match_score
            else:
                match = dp[i - 1][j - 1] + mismatch_penalty
            
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            
            # Key difference: allow zero (no negative scores)
            dp[i][j] = max(0, match, delete, insert)
            
            # Track maximum score position
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_pos = (i, j)
    
    # Traceback from maximum score position
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos
    
    while i > 0 and j > 0 and dp[i][j] > 0:
        current_score = dp[i][j]
        diagonal_score = dp[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
        
        if current_score == diagonal_score:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == dp[i - 1][j] + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
    
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    return ''.join(aligned_seq1), ''.join(aligned_seq2), max_score


if __name__ == "__main__":
    # Example usage
    seq1 = "GGTTGACTA"
    seq2 = "TGTTACGG"
    
    aligned1, aligned2, score = smith_waterman(seq1, seq2)
    
    print("Smith-Waterman Local Alignment")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
