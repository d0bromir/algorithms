"""
Simple tests for bioinformatics algorithms.

These tests verify that the basic functionality works correctly.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman
from seed_and_extend import seed_and_extend, KmerIndex
from bwt_fm_index import FMIndex, burrows_wheeler_transform, inverse_bwt


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
    aligned1, aligned2, score1 = needleman_wunsch("AC", "AC", match_score=2)
    aligned1, aligned2, score2 = needleman_wunsch("AC", "AC", match_score=1)
    assert score1 > score2, "Higher match score should give higher total score"
    
    # Test that local alignment score >= 0
    aligned1, aligned2, score = smith_waterman("AAAA", "TTTT")
    assert score >= 0, "Smith-Waterman score should never be negative"
    
    print("  ✓ Algorithm properties tests passed")


def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 50)
    print("Running Bioinformatics Algorithm Tests")
    print("=" * 50 + "\n")
    
    try:
        test_needleman_wunsch()
        test_smith_waterman()
        test_seed_and_extend()
        test_bwt_fm_index()
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
