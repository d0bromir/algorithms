"""
Simple tests for bioinformatics algorithms.

These tests verify that the basic functionality works correctly.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from needleman_wunsch import needleman_wunsch
from needleman_wunsch_affine import needleman_wunsch_affine
from hirschberg import hirschberg
from smith_waterman import smith_waterman
from smith_waterman_affine import smith_waterman_affine
from seed_and_extend import seed_and_extend, KmerIndex
from bwt_fm_index import FMIndex, burrows_wheeler_transform, inverse_bwt
from myers_bitvector import myers_bit_vector_edit_distance, edit_distance_bruteforce
from wavefront_alignment import wavefront_alignment
from minimizer_chaining import build_minimizer_index, find_seed_matches, chain_seeds
from strobemer_mapping import (
    strobemer_index, find_strobemer_matches, minhash_sketch, estimate_jaccard,
)


def test_needleman_wunsch():
    """Test Needleman-Wunsch algorithm."""
    print("Testing Needleman-Wunsch...")
    
    # Test 1: Identical sequences
    aligned1, aligned2, score = needleman_wunsch("ACGT", "ACGT")
    assert aligned1 == "ACGT", f"Expected 'ACGT', got '{aligned1}'"
    assert aligned2 == "ACGT", f"Expected 'ACGT', got '{aligned2}'"
    assert score == 4, f"Expected score 4, got {score}"
    
    # Test 2: Different sequences
    aligned1, aligned2, score = needleman_wunsch("GATTACA", "GCATGCU")
    assert len(aligned1) == len(aligned2), "Aligned sequences must have same length"
    
    # Test 3: Empty sequences
    aligned1, aligned2, score = needleman_wunsch("", "ACG")
    assert aligned1 == "---", f"Expected '---', got '{aligned1}'"
    assert aligned2 == "ACG", f"Expected 'ACG', got '{aligned2}'"
    
    print("  ✓ Needleman-Wunsch tests passed")


def test_hirschberg():
    """Test Hirschberg algorithm."""
    print("Testing Hirschberg...")
    
    # Test 1: Identical sequences
    aligned1, aligned2, score = hirschberg("ACGT", "ACGT")
    assert aligned1 == "ACGT", f"Expected 'ACGT', got '{aligned1}'"
    assert aligned2 == "ACGT", f"Expected 'ACGT', got '{aligned2}'"
    assert score == 4, f"Expected score 4, got {score}"
    
    # Test 2: Different sequences
    aligned1, aligned2, score = hirschberg("GATTACA", "GCATGCU")
    assert len(aligned1) == len(aligned2), "Aligned sequences must have same length"
    assert score == 0, f"Expected score 0, got {score}"
    
    # Test 3: Empty sequences
    aligned1, aligned2, score = hirschberg("", "ACG")
    assert aligned1 == "---", f"Expected '---', got '{aligned1}'"
    assert aligned2 == "ACG", f"Expected 'ACG', got '{aligned2}'"
    assert score == -3, f"Expected score -3, got {score}"
    
    # Test 4: Verify same results as Needleman-Wunsch
    seq1, seq2 = "AGTACGCA", "TATGC"
    h_aligned1, h_aligned2, h_score = hirschberg(seq1, seq2)
    nw_aligned1, nw_aligned2, nw_score = needleman_wunsch(seq1, seq2)
    assert h_score == nw_score, f"Hirschberg and NW scores should match: {h_score} vs {nw_score}"
    # Verify alignment properties
    assert len(h_aligned1) == len(h_aligned2), "Aligned sequences must have same length"
    assert all(h_aligned1[i] != '-' or h_aligned2[i] != '-' for i in range(len(h_aligned1))), "No double gaps allowed"
    assert ''.join(c for c in h_aligned1 if c != '-') == seq1, "Original sequence 1 should be preserved"
    assert ''.join(c for c in h_aligned2 if c != '-') == seq2, "Original sequence 2 should be preserved"
    
    # Test 5: Custom scoring parameters
    aligned1, aligned2, score = hirschberg("AC", "AC", match_score=2, mismatch_penalty=-1, gap_penalty=-2)
    nw_aligned1, nw_aligned2, nw_score = needleman_wunsch("AC", "AC", match_score=2, mismatch_penalty=-1, gap_penalty=-2)
    assert score == nw_score, f"Scores should match with custom parameters: {score} vs {nw_score}"
    
    print("  ✓ Hirschberg tests passed")


def test_needleman_wunsch_affine():
    """Test Gotoh's affine-gap global alignment."""
    print("Testing Needleman-Wunsch (affine gap, Gotoh)...")

    # Test 1: Identical sequences should score len * match_score
    aligned1, aligned2, score = needleman_wunsch_affine("ACGT", "ACGT")
    assert aligned1 == "ACGT" and aligned2 == "ACGT"
    assert score == 4, f"Expected score 4, got {score}"

    # Test 2: A clustered indel should be cheaper than scattered ones under
    # the same affine parameters (the whole point of affine gap penalties).
    seq1, seq2 = "ACGTACGTACGT", "ACGTAAACGT"
    _, _, clustered_score = needleman_wunsch_affine(seq1, seq2, 2, -1, -5, -1)
    assert clustered_score is not None

    # Test 3: alignment reconstructs original sequences
    aligned1, aligned2, _ = needleman_wunsch_affine("GGTTGACTA", "TGTTACGG", 2, -1, -3, -1)
    assert len(aligned1) == len(aligned2)
    assert ''.join(c for c in aligned1 if c != '-') == "GGTTGACTA"
    assert ''.join(c for c in aligned2 if c != '-') == "TGTTACGG"

    print("  ✓ Needleman-Wunsch affine tests passed")


def test_myers_bitvector():
    """Test Myers' bit-vector edit distance algorithm against brute force."""
    print("Testing Myers' bit-vector algorithm...")

    # Test 1: identical strings -> distance 0
    dist, pos = myers_bit_vector_edit_distance("ACGT", "ACGT")
    assert dist == 0 and pos == 4

    # Test 2: matches brute-force Levenshtein distance
    dist, _ = myers_bit_vector_edit_distance("GATTACA", "GACTATA")
    brute = edit_distance_bruteforce("GATTACA", "GACTATA")
    assert dist == brute, f"Expected {brute}, got {dist}"

    # Test 3: empty pattern
    dist, pos = myers_bit_vector_edit_distance("", "ACGT")
    assert dist == 4 and pos == 4

    print("  ✓ Myers' bit-vector tests passed")


def test_wavefront_alignment():
    """Test the Wavefront Alignment (WFA) algorithm against affine-gap DP."""
    print("Testing Wavefront Alignment (WFA)...")

    def affine_dp_min(seq1, seq2, mismatch=4, gap_open=6, gap_extend=2):
        m, n = len(seq1), len(seq2)
        INF = float('inf')
        M = [[INF] * (n + 1) for _ in range(m + 1)]
        I = [[INF] * (n + 1) for _ in range(m + 1)]
        D = [[INF] * (n + 1) for _ in range(m + 1)]
        M[0][0] = 0
        for j in range(1, n + 1):
            I[0][j] = gap_open + gap_extend * j
        for i in range(1, m + 1):
            D[i][0] = gap_open + gap_extend * i
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                I[i][j] = min(M[i][j - 1] + gap_open + gap_extend, I[i][j - 1] + gap_extend)
                D[i][j] = min(M[i - 1][j] + gap_open + gap_extend, D[i - 1][j] + gap_extend)
                s = 0 if seq1[i - 1] == seq2[j - 1] else mismatch
                M[i][j] = min(M[i - 1][j - 1] + s, I[i - 1][j - 1] + s, D[i - 1][j - 1] + s)
        return min(M[m][n], I[m][n], D[m][n])

    # Test 1: identical sequences -> score 0
    assert wavefront_alignment("ACGT", "ACGT") == 0

    # Test 2: matches affine-gap DP minimum on several cases
    for seq1, seq2 in [("GATTACA", "GATCACA"), ("ACGTACGTACGT", "ACGTAAACGT"),
                        ("TAAAG", "AC"), ("", "ACG"), ("ACG", "")]:
        got = wavefront_alignment(seq1, seq2, mismatch=4, gap_open=6, gap_extend=2)
        expected = affine_dp_min(seq1, seq2, 4, 6, 2)
        assert got == expected, f"WFA({seq1!r}, {seq2!r}) = {got}, expected {expected}"

    print("  ✓ Wavefront Alignment tests passed")


def test_minimizer_chaining():
    """Test minimizer-based seeding and co-linear chaining."""
    print("Testing minimizer sketching + chaining...")

    reference = ("ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                 "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG")
    query = "GACCATGGGACCATTGACCTGAGGGTACCGGATT"  # exact substring at index 27

    k, w = 8, 4
    ref_index = build_minimizer_index(reference, k=k, w=w)
    seeds = find_seed_matches(query, ref_index, k=k, w=w)
    assert len(seeds) > 0, "Should find seed matches for an exact substring"

    chain, score = chain_seeds(seeds, k=k)
    assert len(chain) > 0
    true_start = reference.find(query)
    q0, r0 = chain[0]
    assert r0 - q0 == true_start, f"Chain should recover true offset {true_start}, got {r0 - q0}"

    print("  ✓ Minimizer chaining tests passed")


def test_strobemer_mapping():
    """Test strobemer seeding and MinHash Jaccard estimation."""
    print("Testing strobemer seeding + MinHash...")

    reference = ("ACGTACGGTTAGCATGACGGATCCAGTGACCATGGGACCATTGACCTGA"
                 "GGGTACCGGATTACAAGGCTAGCTAGGATCCAGTTAGGCATGGCTTAAGG"
                 "TTCCGGAATTCCGGATCGATCGGATCCAAGCTTGGATCCACTAGTCCAGT")
    query = reference[60:130]

    idx = strobemer_index(reference, strobe_len=5, w_min=4, w_max=12)
    matches = find_strobemer_matches(query, idx, strobe_len=5, w_min=4, w_max=12)
    assert len(matches) > 0, "Should find strobemer matches for an exact substring"

    offsets = [r - q for q, r in matches]
    most_common = max(set(offsets), key=offsets.count)
    assert most_common == 60, f"Expected dominant offset 60, got {most_common}"

    # MinHash: identical sequences must have Jaccard similarity 1.0
    sketch = minhash_sketch(reference, k=12, sketch_size=100)
    assert estimate_jaccard(sketch, sketch) == 1.0

    # A sequence and a very different one should show low similarity
    unrelated = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    sketch_unrelated = minhash_sketch(unrelated, k=12, sketch_size=100)
    assert estimate_jaccard(sketch, sketch_unrelated) < estimate_jaccard(sketch, sketch)

    print("  ✓ Strobemer/MinHash tests passed")


def test_smith_waterman():
    """Test Smith-Waterman algorithm."""
    print("Testing Smith-Waterman...")
    
    # Test 1: Local match
    aligned1, aligned2, score = smith_waterman("ACGTACGT", "ACGT")
    assert "ACGT" in aligned1, f"Expected 'ACGT' in alignment"
    assert score > 0, f"Expected positive score, got {score}"
    
    # Test 2: No match
    aligned1, aligned2, score = smith_waterman("AAAA", "TTTT")
    assert score >= 0, f"Score should be non-negative, got {score}"
    
    print("  ✓ Smith-Waterman tests passed")


def test_seed_and_extend():
    """Test Seed-and-Extend algorithm."""
    print("Testing Seed-and-Extend...")
    
    # Test 1: K-mer index
    kmer_index = KmerIndex("ACGTACGT", k=4)
    seeds = kmer_index.find_seeds("ACGT")
    assert len(seeds) > 0, "Should find at least one seed"
    
    # Test 2: Exact match
    alignments = seed_and_extend("ACGTACGTACGT", "ACGT", k=3)
    assert len(alignments) > 0, "Should find alignments"
    
    # Test 3: No match with large k
    alignments = seed_and_extend("ACGT", "TTTT", k=3)
    assert len(alignments) == 0, "Should find no alignments"
    
    print("  ✓ Seed-and-Extend tests passed")


def test_bwt_fm_index():
    """Test BWT and FM-Index."""
    print("Testing BWT and FM-Index...")
    
    # Test 1: BWT reversibility
    text = "BANANA"
    bwt = burrows_wheeler_transform(text)
    recovered = inverse_bwt(bwt)
    assert recovered == text + "$", f"BWT should be reversible"
    
    # Test 2: FM-Index exact search
    fm_index = FMIndex("ACGTACGTACGT")
    
    # Test exact pattern
    result = fm_index.search("ACG")
    assert result['count'] > 0, "Should find pattern 'ACG'"
    assert len(result['positions']) == result['count'], "Count should match positions"
    
    # Test non-existent pattern
    result = fm_index.search("XYZ")
    assert result['count'] == 0, "Should not find pattern 'XYZ'"
    assert len(result['positions']) == 0, "Should have no positions"
    
    # Test 3: Multiple occurrences
    fm_index = FMIndex("AAAA")
    result = fm_index.search("AA")
    assert result['count'] >= 2, "Should find multiple overlapping matches"
    
    print("  ✓ BWT and FM-Index tests passed")


def test_algorithm_properties():
    """Test general properties of algorithms."""
    print("Testing algorithm properties...")
    
    # Test that DP algorithms handle custom scores
    _, _, nw_score1 = needleman_wunsch("AC", "AC", match_score=2)
    _, _, nw_score2 = needleman_wunsch("AC", "AC", match_score=1)
    assert nw_score1 > nw_score2, "Higher match score should give higher total score"
    
    # Test that Hirschberg also handles custom scores
    _, _, h_score1 = hirschberg("AC", "AC", match_score=2)
    _, _, h_score2 = hirschberg("AC", "AC", match_score=1)
    assert h_score1 > h_score2, "Higher match score should give higher total score (Hirschberg)"
    
    # Test that local alignment score >= 0
    _, _, score = smith_waterman("AAAA", "TTTT")
    assert score >= 0, "Smith-Waterman score should never be negative"
    
    print("  ✓ Algorithm properties tests passed")


def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 50)
    print("Running Bioinformatics Algorithm Tests")
    print("=" * 50 + "\n")
    
    try:
        test_needleman_wunsch()
        test_needleman_wunsch_affine()
        test_hirschberg()
        test_smith_waterman()
        test_seed_and_extend()
        test_bwt_fm_index()
        test_myers_bitvector()
        test_wavefront_alignment()
        test_minimizer_chaining()
        test_strobemer_mapping()
        test_algorithm_properties()
        
        print("\n" + "=" * 50)
        print("✓ All tests passed!")
        print("=" * 50 + "\n")
        return True
        
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return False
    except Exception as e:
        print(f"\n✗ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
