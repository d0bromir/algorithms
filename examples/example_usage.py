"""
Example usage of all bioinformatics algorithms.

This file demonstrates how to use each algorithm with realistic examples.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman
from seed_and_extend import seed_and_extend
from bwt_fm_index import FMIndex, burrows_wheeler_transform


def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70 + "\n")


def example_needleman_wunsch():
    """Example: Global sequence alignment."""
    print_section("Needleman-Wunsch Global Alignment")
    
    # Example 1: DNA sequences
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    
    aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
    
    print("DNA Sequence Alignment:")
    print(f"  Sequence 1: {seq1}")
    print(f"  Sequence 2: {seq2}")
    print(f"\n  Aligned 1:  {aligned1}")
    print(f"  Aligned 2:  {aligned2}")
    print(f"  Score: {score}")
    
    # Show matching positions
    matches = ''.join(['|' if aligned1[i] == aligned2[i] else ' ' 
                      for i in range(len(aligned1))])
    print(f"  Matches:    {matches}")
    
    # Example 2: Protein sequences
    print("\n" + "-" * 70)
    protein1 = "PAWHEAE"
    protein2 = "HEAGAWGHEE"
    
    aligned1, aligned2, score = needleman_wunsch(protein1, protein2)
    
    print("\nProtein Sequence Alignment:")
    print(f"  Sequence 1: {protein1}")
    print(f"  Sequence 2: {protein2}")
    print(f"\n  Aligned 1:  {aligned1}")
    print(f"  Aligned 2:  {aligned2}")
    print(f"  Score: {score}")


def example_smith_waterman():
    """Example: Local sequence alignment."""
    print_section("Smith-Waterman Local Alignment")
    
    # Example: Finding conserved region
    seq1 = "GGTTGACTA"
    seq2 = "TGTTACGG"
    
    aligned1, aligned2, score = smith_waterman(seq1, seq2)
    
    print("Finding Best Local Match:")
    print(f"  Sequence 1: {seq1}")
    print(f"  Sequence 2: {seq2}")
    print(f"\n  Local Alignment:")
    print(f"    Region 1: {aligned1}")
    print(f"    Region 2: {aligned2}")
    print(f"    Score: {score}")
    
    # Show matching positions
    matches = ''.join(['|' if aligned1[i] == aligned2[i] else ' ' 
                      for i in range(len(aligned1))])
    print(f"    Matches:  {matches}")
    
    # Example 2: Different scoring
    print("\n" + "-" * 70)
    seq1 = "ACACACTA"
    seq2 = "AGCACACA"
    
    aligned1, aligned2, score = smith_waterman(
        seq1, seq2, 
        match_score=3, 
        mismatch_penalty=-2, 
        gap_penalty=-2
    )
    
    print("\nWith Custom Scoring:")
    print(f"  Match: +3, Mismatch: -2, Gap: -2")
    print(f"  Sequence 1: {seq1}")
    print(f"  Sequence 2: {seq2}")
    print(f"\n  Local Alignment:")
    print(f"    Region 1: {aligned1}")
    print(f"    Region 2: {aligned2}")
    print(f"    Score: {score}")


def example_seed_and_extend():
    """Example: Fast alignment with k-mers."""
    print_section("Seed-and-Extend with K-mer Hashing")
    
    # Example: Short read alignment
    reference = "ACGTACGTACGTAAACCCGGGTTTACGTACGTCCGGTTAA"
    query = "ACGTAAACCCGGG"
    
    print("Aligning Short Read to Reference:")
    print(f"  Reference: {reference}")
    print(f"  Query:     {query}")
    print(f"  K-mer size: 5")
    
    alignments = seed_and_extend(reference, query, k=5)
    
    print(f"\n  Found {len(alignments)} seed(s)")
    
    # Show top alignment
    if alignments:
        best = alignments[0]
        print(f"\n  Best Alignment:")
        print(f"    Query position:  {best['query_start']}-{best['query_end']}")
        print(f"    Ref position:    {best['ref_start']}-{best['ref_end']}")
        print(f"    Score:           {best['score']}")
        print(f"    Query sequence:  {best['query_seq']}")
        print(f"    Ref sequence:    {best['ref_seq']}")
        
        # Show alignment in context
        print(f"\n  Alignment in Context:")
        context_start = max(0, best['ref_start'] - 5)
        context_end = min(len(reference), best['ref_end'] + 5)
        context = reference[context_start:context_end]
        offset = best['ref_start'] - context_start
        
        print(f"    Reference: {context}")
        print(f"               {' ' * offset}{best['ref_seq']}")
        print(f"    Query:     {' ' * offset}{best['query_seq']}")


def example_bwt_fm_index():
    """Example: Fast exact pattern matching."""
    print_section("BWT + FM-Index Pattern Matching")
    
    # Example: Genome indexing and search
    genome = "ACGTACGTACGTTAGCTAGCTAGCT"
    
    print("Indexing Genome Sequence:")
    print(f"  Genome: {genome}")
    print(f"  Length: {len(genome)} bp")
    
    # Compute BWT
    bwt = burrows_wheeler_transform(genome)
    print(f"  BWT:    {bwt}")
    
    # Build FM-index
    fm_index = FMIndex(genome)
    
    # Search for patterns
    print("\n  Pattern Searches:")
    patterns = ["ACG", "TAG", "GCT", "ACGT", "XYZ", "TAGCT"]
    
    for pattern in patterns:
        result = fm_index.search(pattern)
        print(f"\n    Pattern: '{pattern}'")
        print(f"      Count: {result['count']}")
        if result['positions']:
            print(f"      Positions: {result['positions']}")
            
            # Show first occurrence in context
            pos = result['positions'][0]
            context_start = max(0, pos - 3)
            context_end = min(len(genome), pos + len(pattern) + 3)
            before = genome[context_start:pos]
            match = genome[pos:pos + len(pattern)]
            after = genome[pos + len(pattern):context_end]
            
            print(f"      Context: ...{before}[{match}]{after}...")


def example_comparison():
    """Example: Comparing different algorithms."""
    print_section("Algorithm Comparison")
    
    # Use same sequences for all algorithms
    reference = "ACGTACGTTAGCTAGCT"
    query = "ACGTTAGC"
    
    print(f"Reference: {reference}")
    print(f"Query:     {query}")
    
    # Global alignment
    print("\n1. Needleman-Wunsch (Global):")
    aligned1, aligned2, score = needleman_wunsch(query, reference)
    print(f"   Score: {score}")
    print(f"   Aligned query: {aligned1}")
    print(f"   Aligned ref:   {aligned2}")
    
    # Local alignment
    print("\n2. Smith-Waterman (Local):")
    aligned1, aligned2, score = smith_waterman(query, reference)
    print(f"   Score: {score}")
    print(f"   Best local match:")
    print(f"     Query region: {aligned1}")
    print(f"     Ref region:   {aligned2}")
    
    # Seed-and-extend
    print("\n3. Seed-and-Extend:")
    alignments = seed_and_extend(reference, query, k=4)
    if alignments:
        best = alignments[0]
        print(f"   Score: {best['score']}")
        print(f"   Position in ref: {best['ref_start']}-{best['ref_end']}")
        print(f"   Matched: {best['ref_seq']}")
    
    # FM-index
    print("\n4. FM-Index (Exact Match):")
    fm_index = FMIndex(reference)
    result = fm_index.search(query)
    if result['count'] > 0:
        print(f"   Exact matches: {result['count']}")
        print(f"   Positions: {result['positions']}")
    else:
        # Try shorter exact match
        shorter_query = query[:6]
        result = fm_index.search(shorter_query)
        print(f"   No exact match for full query")
        print(f"   Searching for '{shorter_query}': {result['count']} match(es)")
        if result['positions']:
            print(f"   Positions: {result['positions']}")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("  BIOINFORMATICS ALGORITHMS - Example Usage")
    print("=" * 70)
    
    try:
        example_needleman_wunsch()
        example_smith_waterman()
        example_seed_and_extend()
        example_bwt_fm_index()
        example_comparison()
        
        print("\n" + "=" * 70)
        print("  All examples completed successfully!")
        print("=" * 70 + "\n")
        
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
