# Bioinformatics Algorithms Documentation

## Algorithm Details

This document provides detailed information about each algorithm implementation.

## 1. Needleman-Wunsch (Global Alignment)

### Overview
The Needleman-Wunsch algorithm is a dynamic programming algorithm used for global sequence alignment. It finds the optimal alignment between two complete sequences by maximizing the alignment score.

### Algorithm Steps
1. **Initialize** - Create a (m+1) × (n+1) matrix where m and n are sequence lengths
2. **Fill Matrix** - For each cell (i,j), compute:
   - Match/Mismatch: score[i-1,j-1] + match_score or mismatch_penalty
   - Deletion: score[i-1,j] + gap_penalty
   - Insertion: score[i,j-1] + gap_penalty
   - Take maximum of these three values
3. **Traceback** - Start from bottom-right, follow path of maximum scores to top-left

### Parameters
- `match_score`: Score added for matching characters (default: 1)
- `mismatch_penalty`: Score subtracted for mismatches (default: -1)
- `gap_penalty`: Score subtracted for gaps (default: -1)

### Time and Space Complexity
- Time: O(m × n)
- Space: O(m × n)

### Use Cases
- Comparing complete protein sequences
- Aligning two genes to find evolutionary relationships
- When you need to align entire sequences end-to-end

### Limitations
- Slow for large sequences (>10,000 bp)
- Requires quadratic memory
- Forces alignment of entire sequences even if only part matches

## 2. Smith-Waterman (Local Alignment)

### Overview
The Smith-Waterman algorithm is a dynamic programming algorithm for local sequence alignment. It finds the best matching subsequence regions between two sequences.

### Key Differences from Needleman-Wunsch
1. Matrix initialized with zeros (not gap penalties)
2. Negative scores reset to zero
3. Traceback starts from maximum score (not bottom-right)
4. Traceback stops when reaching zero

### Algorithm Steps
1. **Initialize** - Create (m+1) × (n+1) matrix, all cells = 0
2. **Fill Matrix** - Same as Needleman-Wunsch but take max(0, match, delete, insert)
3. **Find Maximum** - Track position of highest score
4. **Traceback** - Start from maximum, stop at zero

### Parameters
- `match_score`: Score for matches (default: 2, higher than Needleman-Wunsch)
- `mismatch_penalty`: Penalty for mismatches (default: -1)
- `gap_penalty`: Penalty for gaps (default: -1)

### Time and Space Complexity
- Time: O(m × n)
- Space: O(m × n)

### Use Cases
- Finding conserved domains in proteins
- Identifying similar regions in divergent sequences
- Database searches where you want local matches
- Finding motifs or patterns

### Advantages over Needleman-Wunsch
- Better for divergent sequences
- Ignores poorly matching regions
- More sensitive for partial matches

## 2.5. Smith-Waterman with Affine Gap Penalties

### Overview
An enhanced version of the Smith-Waterman algorithm that uses affine gap penalties instead of linear gap penalties. Affine gaps are more biologically realistic as they distinguish between opening a gap (more costly) and extending an existing gap (less costly).

### Affine Gap Model
- **Linear gap penalty** (basic version): gap_penalty × gap_length
- **Affine gap penalty** (this version): gap_open + (gap_extend × gap_length)

For example, a gap of length 3:
- Linear: -1 × 3 = -3
- Affine: -3 + (-1 × 3) = -6 (gap opening is more costly)

### Three-Matrix Approach
Uses three dynamic programming matrices:
1. **M[i][j]**: Best score ending with match/mismatch at position (i,j)
2. **I[i][j]**: Best score ending with insertion (gap in seq1) at position (i,j)
3. **D[i][j]**: Best score ending with deletion (gap in seq2) at position (i,j)

### Algorithm Steps
1. **Initialize** - Create three (m+1) × (n+1) matrices
   - M initialized with 0
   - I and D initialized with negative infinity (except allowing gaps from 0)
2. **Fill Matrices** - For each cell (i,j):
   - I[i][j] = max(0, M[i][j-1] + gap_open + gap_extend, I[i][j-1] + gap_extend, D[i][j-1] + gap_open + gap_extend)
   - D[i][j] = max(0, M[i-1][j] + gap_open + gap_extend, D[i-1][j] + gap_extend, I[i-1][j] + gap_open + gap_extend)
   - M[i][j] = max(0, M[i-1][j-1] + s(i,j), I[i-1][j-1] + s(i,j), D[i-1][j-1] + s(i,j))
   - Where s(i,j) is match_score or mismatch_penalty
3. **Find Maximum** - Track highest score across all three matrices
4. **Traceback** - Follow the path through appropriate matrices until score reaches 0

### Parameters
- `match_score`: Score for matching characters (default: 2)
- `mismatch_penalty`: Penalty for mismatches (default: -1)
- `gap_open`: Penalty for opening a new gap (default: -3)
- `gap_extend`: Penalty for extending existing gap (default: -1)

### Time and Space Complexity
- Time: O(m × n) - same as basic Smith-Waterman
- Space: O(3 × m × n) = O(m × n) - three matrices instead of one

### Use Cases
- Protein sequence alignment (gaps often appear in runs)
- Aligning sequences with indels
- When biological accuracy is important
- Sequences where insertions/deletions tend to cluster

### Advantages over Linear Gap Penalty
- More biologically realistic
- Penalizes scattered gaps more than clustered gaps
- Better models evolutionary insertions/deletions
- Produces cleaner alignments with fewer scattered gaps

### When to Use Affine vs Linear
- **Use Affine** when:
  - Working with biological sequences (especially proteins)
  - You expect runs of insertions/deletions
  - Alignment quality is critical
- **Use Linear** when:
  - Speed is more important than accuracy
  - Gap structure is not important
  - Quick similarity scoring is needed

### Typical Parameter Values
**For DNA:**
- match_score: 1 to 5
- mismatch_penalty: -1 to -4
- gap_open: -3 to -10
- gap_extend: -1 to -2

**For Proteins:**
- Use substitution matrix (BLOSUM62, PAM250)
- gap_open: -10 to -12
- gap_extend: -1 to -2

## 3. Seed-and-Extend (K-mer Hashing)

### Overview
A heuristic algorithm that uses exact k-mer matches as "seeds" and extends them to find longer alignments. This is the basis for BLAST, MAQ, and SOAP.

### Algorithm Steps
1. **Index Reference** - Create hash table of all k-mers in reference
   - Key: k-mer sequence
   - Value: list of positions where k-mer occurs
2. **Find Seeds** - For each k-mer in query, look up in hash table
3. **Extend Seeds** - For each seed, extend in both directions:
   - Add match_score for matches
   - Add mismatch_penalty for mismatches
   - Stop when score drops too far below maximum
4. **Rank Results** - Sort alignments by score

### Parameters
- `k`: K-mer length (default: 11 for DNA)
  - Larger k: fewer false positives, faster, less sensitive
  - Smaller k: more sensitive, slower, more false positives
- `match_score`: Score for matches during extension
- `mismatch_penalty`: Penalty for mismatches during extension

### Time and Space Complexity
- Indexing: O(n) time, O(n) space for reference of length n
- Query: O(m) time for query of length m (average case)
- Worst case: O(m × n) if many k-mer matches

### K-mer Size Guidelines
- **k=8**: Very sensitive, many false hits, slow
- **k=11**: Good balance for short reads (BLAST default for DNA)
- **k=15**: Fast, good for similar sequences
- **k=20+**: Very fast, only finds highly similar regions

### Use Cases
- Fast database searches (BLAST)
- Aligning millions of short reads
- Finding similar sequences in large databases
- Pre-filtering before exact alignment

### Advantages
- Much faster than dynamic programming
- Scales to large databases
- Tunable sensitivity/speed tradeoff

### Limitations
- Heuristic, not guaranteed optimal
- May miss alignments without exact k-mer matches
- Sensitive to k-mer choice

## 4. Burrows-Wheeler Transform (BWT) + FM-Index

### Overview
The BWT reorganizes text to make it more compressible, and the FM-index uses this for ultra-fast exact pattern matching. Used in BWA, Bowtie, and HISAT2.

### Burrows-Wheeler Transform

#### Steps
1. Add sentinel character '$' (lexicographically smallest)
2. Generate all rotations of the text
3. Sort rotations lexicographically
4. BWT is the last column of sorted rotations

#### Properties
- Reversible transformation
- Groups similar characters together
- Enables compression and fast search

### FM-Index

#### Components
1. **BWT**: The transformed text
2. **C Array**: Count of characters lexicographically smaller than each character
3. **Occurrence Array (Occ)**: Count of each character up to each position

#### Backward Search Algorithm
For pattern P = p₁p₂...pₘ:
1. Start with full range [0, n-1]
2. For each character from right to left:
   - top = C[c] + Occ[c][top]
   - bottom = C[c] + Occ[c][bottom+1] - 1
3. If top > bottom, pattern not found
4. Otherwise, [top, bottom] gives range in suffix array

### Time and Space Complexity
- Construction: O(n log n) time, O(n) space
- Search: O(m) time for pattern of length m (independent of text size!)
- Space: O(n) with compression possible

### Parameters
- Text can be preprocessed once
- Search is parameter-free (exact match only)

### Use Cases
- Aligning millions of short reads to genome (BWA, Bowtie)
- Finding all exact occurrences of pattern
- Compressed full-text search
- When you need to search same reference many times

### Advantages
- Search time independent of text length
- Very memory efficient with compression
- Supports backward search and other advanced queries
- Can find ALL occurrences in O(m) time

### Extensions
- **Seeding**: Use FM-index to find exact matches (seeds)
- **MEMs**: Maximal Exact Matches
- **SMEMs**: Super-Maximal Exact Matches
- **Inexact Matching**: Allow mismatches with backtracking

### Limitations
- Only finds exact matches (extensions needed for mismatches)
- Construction slower than simple hashing
- More complex to implement than other methods

## Choosing the Right Algorithm

### Decision Tree

```
Are sequences very similar (>95% identity)?
├─ YES → Use Seed-and-Extend or FM-Index
│         - Very fast
│         - Good for reads alignment
│
└─ NO → Are sequences short (<1000 bp)?
        ├─ YES → Use Smith-Waterman or Needleman-Wunsch
        │         - Optimal alignment
        │         - Handles divergent sequences
        │
        └─ NO → Use Seed-and-Extend
                  - Dynamic programming too slow
                  - May need to reduce k for sensitivity
```

### By Use Case

| Use Case | Algorithm | Why |
|----------|-----------|-----|
| Database search | Seed-and-Extend | Fast, scalable |
| NGS read alignment | BWT + FM-Index | Ultra-fast exact matching |
| Protein comparison | Smith-Waterman (Affine) | Finds functional domains, realistic gaps |
| Gene comparison | Needleman-Wunsch | Complete gene alignment |
| Finding motifs | Smith-Waterman | Local pattern matching |
| SNP calling | BWT + FM-Index | Align millions of reads |
| Sequences with indels | Smith-Waterman (Affine) | Better models clustered gaps |

### By Sequence Properties

| Property | Best Algorithm |
|----------|---------------|
| Very long (>100 Mbp) | BWT + FM-Index |
| Long (10K-100K bp) | Seed-and-Extend |
| Medium (1K-10K bp) | Smith-Waterman |
| Short (<1K bp) | Needleman-Wunsch or Smith-Waterman |
| Highly similar | BWT + FM-Index or Seed-and-Extend |
| Divergent | Smith-Waterman |

## Performance Characteristics

### Typical Running Times (approximate)

For aligning two sequences of length n:

| Algorithm | n=100 | n=1,000 | n=10,000 | n=100,000 |
|-----------|-------|---------|----------|-----------|
| Needleman-Wunsch | <1ms | 10ms | 1s | 100s |
| Smith-Waterman | <1ms | 10ms | 1s | 100s |
| Seed-and-Extend | <1ms | 1ms | 10ms | 100ms |
| FM-Index (search) | <1ms | <1ms | <1ms | <1ms |

*Note: FM-Index construction is O(n log n), but search is O(m) independent of reference size*

## References and Further Reading

1. Needleman, S.B. & Wunsch, C.D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology.

2. Smith, T.F. & Waterman, M.S. (1981). "Identification of common molecular subsequences". Journal of Molecular Biology.

3. Altschul, S.F. et al. (1990). "Basic local alignment search tool (BLAST)". Journal of Molecular Biology.

4. Burrows, M. & Wheeler, D.J. (1994). "A block-sorting lossless data compression algorithm". Technical Report 124, Digital Equipment Corporation.

5. Ferragina, P. & Manzini, G. (2000). "Opportunistic data structures with applications". Proceedings of FOCS.

6. Li, H. & Durbin, R. (2009). "Fast and accurate short read alignment with Burrows-Wheeler transform". Bioinformatics.

7. Langmead, B. et al. (2009). "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome". Genome Biology.
