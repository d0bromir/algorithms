"""
Hirschberg Algorithm for Space-Efficient Global Sequence Alignment

This is a space-efficient divide-and-conquer algorithm for global sequence alignment,
developed by Dan Hirschberg in 1975. It improves upon the Needleman-Wunsch algorithm
by reducing space complexity from O(m * n) to O(min(m, n)) while maintaining the
same O(m * n) time complexity.

The algorithm works by:
1. Using only two rows of the DP matrix at a time (space optimization)
2. Dividing the problem recursively at the midpoint
3. Finding the optimal split point using NW score function
4. Recursively aligning left and right halves
5. Concatenating the results

Time Complexity: O(m * n) where m and n are the lengths of the sequences
Space Complexity: O(min(m, n)) - key improvement over standard Needleman-Wunsch

Reference:
Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common
subsequences. Communications of the ACM, 18(6), 341-343.
"""


def nw_score(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """
    Compute the last row of Needleman-Wunsch scores using only O(n) space.
    
    This function calculates the alignment scores but only keeps the last row
    of the DP matrix, which is sufficient for finding the optimal split point
    in Hirschberg's algorithm.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        match_score: Score for matching characters (default: 1)
        mismatch_penalty: Penalty for mismatching characters (default: -1)
        gap_penalty: Penalty for gaps (default: -1)
    
    Returns:
        list: Last row of DP scores (length n+1)
    """
    m, n = len(seq1), len(seq2)
    
    # Only maintain two rows: previous and current
    prev_row = [j * gap_penalty for j in range(n + 1)]
    curr_row = [0] * (n + 1)
    
    for i in range(1, m + 1):
        curr_row[0] = i * gap_penalty
        
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match = prev_row[j - 1] + match_score
            else:
                match = prev_row[j - 1] + mismatch_penalty
            
            delete = prev_row[j] + gap_penalty
            insert = curr_row[j - 1] + gap_penalty
            
            curr_row[j] = max(match, delete, insert)
        
        # Swap rows for next iteration
        prev_row, curr_row = curr_row, prev_row
    
    return prev_row


def hirschberg(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """
    Perform space-efficient global sequence alignment using Hirschberg algorithm.
    
    This is a divide-and-conquer algorithm that produces the same optimal alignment
    as Needleman-Wunsch but uses only O(min(m, n)) space instead of O(m * n).
    
    The algorithm recursively:
    1. Finds the midpoint of the first sequence
    2. Computes NW scores from both ends to find optimal split in second sequence
    3. Recursively aligns the left halves and right halves
    4. Concatenates the results
    
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
    
    # Base cases
    if m == 0:
        # All gaps in seq1
        return '-' * n, seq2, n * gap_penalty
    
    if n == 0:
        # All gaps in seq2
        return seq1, '-' * m, m * gap_penalty
    
    if m == 1:
        # Single character in seq1 - use simple alignment
        # Try to match seq1[0] with each position in seq2
        best_score = float('-inf')
        best_j = 0
        
        for j in range(n + 1):
            # Calculate score for aligning seq1[0] at position j in seq2
            score = j * gap_penalty  # gaps before
            if j < n:
                if seq1[0] == seq2[j]:
                    score += match_score
                else:
                    score += mismatch_penalty
                score += (n - j - 1) * gap_penalty  # gaps after
            else:
                score += gap_penalty  # seq1[0] aligned to gap
            
            if score > best_score:
                best_score = score
                best_j = j
        
        # Build alignment based on best position
        if best_j == n:
            # Align seq1[0] to gap at end
            aligned1 = '-' * n + seq1[0]
            aligned2 = seq2 + '-'
        else:
            # Align seq1[0] to seq2[best_j]
            aligned1 = '-' * best_j + seq1[0] + '-' * (n - best_j - 1)
            aligned2 = seq2
        
        return aligned1, aligned2, best_score
    
    if n == 1:
        # Single character in seq2 - use simple alignment
        best_score = float('-inf')
        best_i = 0
        
        for i in range(m + 1):
            score = i * gap_penalty
            if i < m:
                if seq1[i] == seq2[0]:
                    score += match_score
                else:
                    score += mismatch_penalty
                score += (m - i - 1) * gap_penalty
            else:
                score += gap_penalty
            
            if score > best_score:
                best_score = score
                best_i = i
        
        if best_i == m:
            aligned1 = seq1 + '-'
            aligned2 = '-' * m + seq2[0]
        else:
            aligned1 = seq1
            aligned2 = '-' * best_i + seq2[0] + '-' * (m - best_i - 1)
        
        return aligned1, aligned2, best_score
    
    # Divide and conquer
    # Find the midpoint of seq1
    mid = m // 2
    
    # Compute NW scores from left (seq1[:mid] vs seq2)
    score_left = nw_score(seq1[:mid], seq2, match_score, mismatch_penalty, gap_penalty)
    
    # Compute NW scores from right (seq1[mid:][::-1] vs seq2[::-1])
    # We reverse both sequences to compute scores from the end
    score_right = nw_score(seq1[mid:][::-1], seq2[::-1], match_score, mismatch_penalty, gap_penalty)
    
    # Reverse score_right to align with seq2
    score_right.reverse()
    
    # Find the split point in seq2 that maximizes total score
    max_score = float('-inf')
    split = 0
    
    for j in range(n + 1):
        total_score = score_left[j] + score_right[j]
        if total_score > max_score:
            max_score = total_score
            split = j
    
    # Recursively align left and right parts
    left1, left2, score1 = hirschberg(
        seq1[:mid], seq2[:split], 
        match_score, mismatch_penalty, gap_penalty
    )
    
    right1, right2, score2 = hirschberg(
        seq1[mid:], seq2[split:],
        match_score, mismatch_penalty, gap_penalty
    )
    
    # Combine the results
    aligned1 = left1 + right1
    aligned2 = left2 + right2
    total_score = score1 + score2
    
    return aligned1, aligned2, total_score


if __name__ == "__main__":
    # Example usage demonstrating the Hirschberg algorithm
    print("Hirschberg Space-Efficient Global Alignment")
    print("=" * 60)
    
    # Example 1: DNA sequence alignment
    print("\nExample 1: DNA Sequence Alignment")
    print("-" * 60)
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    
    aligned1, aligned2, score = hirschberg(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
    
    # Example 2: Protein sequence alignment with custom scoring
    print("\n\nExample 2: Custom Scoring Parameters")
    print("-" * 60)
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    
    # Higher match score, lower penalties
    aligned1, aligned2, score = hirschberg(
        seq1, seq2, 
        match_score=2, 
        mismatch_penalty=-1, 
        gap_penalty=-2
    )
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
    print(f"Parameters: match=2, mismatch=-1, gap=-2")
    
    # Example 3: Comparing with identical sequences
    print("\n\nExample 3: Identical Sequences")
    print("-" * 60)
    seq1 = "ACGTACGT"
    seq2 = "ACGTACGT"
    
    aligned1, aligned2, score = hirschberg(seq1, seq2)
    
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"\nAligned Sequence 1: {aligned1}")
    print(f"Aligned Sequence 2: {aligned2}")
    print(f"Alignment Score: {score}")
    
    # Example 4: Demonstrate space efficiency advantage
    print("\n\nSpace Complexity Comparison")
    print("-" * 60)
    m, n = 1000, 1000
    print(f"For sequences of length {m} and {n}:")
    print(f"  Needleman-Wunsch space: O({m} × {n}) = ~{m*n:,} integers")
    print(f"  Hirschberg space: O(min({m}, {n})) = ~{min(m,n):,} integers")
    print(f"  Space reduction: ~{(m*n)/(min(m,n)):.0f}x improvement")
    
    # Verification example
    print("\n\nExample 5: Verification Against Needleman-Wunsch")
    print("-" * 60)
    
    # Import needleman_wunsch for comparison
    try:
        from needleman_wunsch import needleman_wunsch
        
        seq1 = "AGTACGCA"
        seq2 = "TATGC"
        
        # Run both algorithms
        h_aligned1, h_aligned2, h_score = hirschberg(seq1, seq2)
        nw_aligned1, nw_aligned2, nw_score = needleman_wunsch(seq1, seq2)
        
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")
        print(f"\nHirschberg alignment:")
        print(f"  Aligned 1: {h_aligned1}")
        print(f"  Aligned 2: {h_aligned2}")
        print(f"  Score: {h_score}")
        print(f"\nNeedleman-Wunsch alignment:")
        print(f"  Aligned 1: {nw_aligned1}")
        print(f"  Aligned 2: {nw_aligned2}")
        print(f"  Score: {nw_score}")
        print(f"\nScores match: {h_score == nw_score}")
        
    except ImportError:
        print("needleman_wunsch module not available for comparison")
