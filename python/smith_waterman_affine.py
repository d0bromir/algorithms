"""
Smith-Waterman Algorithm with Affine Gap Penalties for Local Sequence Alignment

This is an advanced version of the Smith-Waterman algorithm that uses affine gap penalties.
Affine gap penalty = gap_open + k * gap_extend (where k is the gap length)
This is more biologically realistic as it penalizes gap opening more than gap extension.

Uses three matrices:
- M[i][j]: optimal score ending with match/mismatch
- I[i][j]: optimal score ending with insertion (gap in seq1)
- D[i][j]: optimal score ending with deletion (gap in seq2)

Time Complexity: O(m * n) where m and n are the lengths of the sequences
Space Complexity: O(m * n)
"""


def smith_waterman_affine(seq1, seq2, match_score=2, mismatch_penalty=-1, 
                          gap_open=-3, gap_extend=-1):
    """
    Perform local sequence alignment using Smith-Waterman algorithm with affine gap penalties.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match_score: Score for matching characters (default: 2)
        mismatch_penalty: Penalty for mismatching characters (default: -1)
        gap_open: Penalty for opening a gap (default: -3)
        gap_extend: Penalty for extending a gap (default: -1)
    
    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)
    """
    m, n = len(seq1), len(seq2)
    
    NEG_INF = float('-inf')
    
    # Initialize three matrices: M (match), I (insert), D (delete)
    M = [[0] * (n + 1) for _ in range(m + 1)]
    I = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    D = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    
    max_score = 0
    max_pos = (0, 0, 'M')  # (i, j, matrix_type)
    
    # Fill DP matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate I[i][j] - gap in seq2 (insertion)
            I[i][j] = max(
                0,
                M[i][j-1] + gap_open + gap_extend,
                I[i][j-1] + gap_extend,
                D[i][j-1] + gap_open + gap_extend
            )
            
            # Calculate D[i][j] - gap in seq1 (deletion)
            D[i][j] = max(
                0,
                M[i-1][j] + gap_open + gap_extend,
                I[i-1][j] + gap_open + gap_extend,
                D[i-1][j] + gap_extend
            )
            
            # Calculate M[i][j] - match/mismatch
            match_mismatch_score = match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty
            M[i][j] = max(
                0,
                M[i-1][j-1] + match_mismatch_score,
                I[i-1][j-1] + match_mismatch_score,
                D[i-1][j-1] + match_mismatch_score
            )
            
            # Track maximum score across all three matrices
            if M[i][j] >= max_score:
                max_score = M[i][j]
                max_pos = (i, j, 'M')
            if I[i][j] >= max_score:
                max_score = I[i][j]
                max_pos = (i, j, 'I')
            if D[i][j] >= max_score:
                max_score = D[i][j]
                max_pos = (i, j, 'D')
    
    # Traceback from maximum score position
    aligned_seq1 = []
    aligned_seq2 = []
    i, j, current_matrix = max_pos
    
    while i > 0 and j > 0:
        if current_matrix == 'M':
            if M[i][j] <= 0:
                break
            
            match_mismatch_score = match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty
            
            if M[i][j] == M[i-1][j-1] + match_mismatch_score:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'M'
            elif M[i][j] == I[i-1][j-1] + match_mismatch_score:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'I'
            elif M[i][j] == D[i-1][j-1] + match_mismatch_score:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'D'
            else:
                break
                
        elif current_matrix == 'I':
            if I[i][j] <= 0:
                break
            
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
            
            if I[i][j] == M[i][j-1] + gap_open + gap_extend:
                current_matrix = 'M'
            elif I[i][j] == I[i][j-1] + gap_extend:
                current_matrix = 'I'
            elif I[i][j] == D[i][j-1] + gap_open + gap_extend:
                current_matrix = 'D'
            else:
                break
                
        else:  # current_matrix == 'D'
            if D[i][j] <= 0:
                break
            
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
            
            if D[i][j] == M[i-1][j] + gap_open + gap_extend:
                current_matrix = 'M'
            elif D[i][j] == D[i-1][j] + gap_extend:
                current_matrix = 'D'
            elif D[i][j] == I[i-1][j] + gap_open + gap_extend:
                current_matrix = 'I'
            else:
                break
    
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    return ''.join(aligned_seq1), ''.join(aligned_seq2), max_score


if __name__ == "__main__":
    # Example usage
    seq1 = "GGTTGACTA"
    seq2 = "TGTTACGG"
    
    print("Smith-Waterman with Affine Gap Penalties - Local Alignment")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()
    
    # Compare linear vs affine gap penalties
    print("With affine gaps (gap_open=-3, gap_extend=-1):")
    aligned1, aligned2, score = smith_waterman_affine(seq1, seq2, 2, -1, -3, -1)
    print(f"Aligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
    print()
    
    # Example with longer sequences showing affine gap advantage
    seq3 = "ACGTACGTACGT"
    seq4 = "ACGTAAACGT"
    
    print("Example 2 - Longer sequences:")
    print(f"Sequence 1: {seq3}")
    print(f"Sequence 2: {seq4}")
    
    aligned3, aligned4, score2 = smith_waterman_affine(seq3, seq4, 2, -1, -5, -1)
    print(f"\nAligned Sequence 1: {aligned3}")
    print(f"Aligned Sequence 2: {aligned4}")
    print(f"Alignment Score: {score2}")
