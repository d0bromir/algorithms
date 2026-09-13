# Bioinformatics Algorithms Documentation

## Algorithm Details

This document provides detailed information about each algorithm implementation.

> For a literature-review-style SWOT analysis (Strengths/Weaknesses/
> Opportunities/Threats) of these and ~35 additional algorithms — SIMD
> libraries, production aligners (Bowtie, BWA-MEM2, minimap2, vg), GPU/
> hardware accelerators, and theoretical lower bounds — with full
> citations, DOIs, and repository links, see [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md).

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

## 2.6. Needleman-Wunsch with Affine Gap Penalties (Gotoh's Algorithm)

### Overview
The global-alignment counterpart to section 2.5: the same three-matrix
(match/insert/delete) recurrence from Gotoh (1982), but without local-alignment
flooring at zero and with traceback always running to (0, 0), since global
alignment must account for the entire length of both sequences.

### Why a Separate M/I/D Boundary Matters
Unlike the local-alignment version, the boundary row/column carry real,
finite costs: `M[0][j]` and `M[i][0]` must be `-infinity` (a match/mismatch
cannot end an alignment that consumed zero characters from the other
sequence), while `I[0][j] = gap_open + j * gap_extend` and
`D[i][0] = gap_open + i * gap_extend` represent "align this whole prefix to
gaps." Getting this boundary wrong is the most common implementation bug in
affine-gap global alignment (an all-zero or aliased boundary silently
produces wrong scores or an invalid traceback).

### Algorithm Steps
1. **Initialize** three `(m+1) x (n+1)` matrices M, I, D to `-infinity`,
   except `M[0][0] = 0`, `I[0][j] = gap_open + j*gap_extend`, and
   `D[i][0] = gap_open + i*gap_extend`.
2. **Fill Matrices** - for i in 1..m, j in 1..n:
   - `I[i][j] = max(M[i][j-1] + gap_open + gap_extend, I[i][j-1] + gap_extend, D[i][j-1] + gap_open + gap_extend)`
   - `D[i][j] = max(M[i-1][j] + gap_open + gap_extend, D[i-1][j] + gap_extend, I[i-1][j] + gap_open + gap_extend)`
   - `M[i][j] = max(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + s(i,j)`
3. **Final Score** - `max(M[m][n], I[m][n], D[m][n])`
4. **Traceback** - from `(m, n)` in whichever matrix achieved the final
   score, back to `(0, 0)`, switching matrices exactly as in section 2.5
   but never stopping early (no zero-flooring).

### Time and Space Complexity
- Time: O(m x n)
- Space: O(m x n)

### Use Cases
- Global alignment of sequences where indels cluster into runs (structural
  variants, whole-gene comparisons with a handful of large indels)
- Any case requiring the biologically-realistic gap model of section 2.5,
  but end-to-end rather than local

### Reference
Gotoh, O. (1982). An improved algorithm for matching biological sequences.
*Journal of Molecular Biology*, 162(3), 705-708.
https://doi.org/10.1016/0022-2836(82)90398-9

Implementation: [`python/needleman_wunsch_affine.py`](python/needleman_wunsch_affine.py),
[`cpp/needleman_wunsch_affine.cpp`](cpp/needleman_wunsch_affine.cpp)
(validated against exhaustive brute-force alignment enumeration).

## 3. Myers' Bit-Vector Algorithm (Edit Distance)

### Overview
Rather than filling the edit-distance DP matrix cell by cell, Myers' (1999)
algorithm packs an *entire column* of the matrix into a couple of
machine-word-sized bit-vectors — `Pv`/`Mv` (where the column's value
increases/decreases by 1 going down) — and updates the whole column with a
fixed sequence of AND/OR/XOR/ADD/shift operations. This gives O(n) time for
patterns up to the machine word size w (64 on modern CPUs), a genuine
constant-factor speedup over scalar DP that makes it the computational core
of the `edlib` library and the fast-extension kernels inside several
high-throughput aligners (see [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md), section B).

### Algorithm Steps
1. **Precompute Peq** - for each alphabet character c, a bitmask with bit i
   set iff `pattern[i] == c`.
2. **Initialize** `Pv = all-ones`, `Mv = 0`, `score = m` (edit distance of
   pattern against the empty text prefix).
3. **For each text character**, in O(1) word operations:
   - `Xv = Peq[c] | Mv`
   - `Xh = (((Peq[c] & Pv) + Pv) ^ Pv) | Peq[c]`
   - `Ph = Mv | ~(Xh | Pv)`, `Mh = Pv & Xh`
   - Update `score` by +1/-1/unchanged based on the top bit of `Ph`/`Mh`
   - Shift `Ph`, `Mh` and recompute `Pv`, `Mv` for the next column
4. **Track the minimum score** seen across all columns to support
   semi-global ("pattern anywhere in text") matching.

### Parameters
- Pattern length m (must be ≤ word size w for the single-word version in
  this repository; production implementations tile ⌈m/w⌉ words for longer
  patterns)

### Time and Space Complexity
- Time: O(n) for m ≤ w, O(n · ⌈m/w⌉) in general
- Space: O(sigma) for Peq, O(1) running state

### Use Cases
- Fast exact edit-distance computation for short-to-medium sequences
  (adapter/primer trimming, k-mer-length exact/near-exact matching)
- The extension kernel inside bit-parallel hardware accelerators (GenASM)

### Reference
Myers, G. (1999). A fast bit-vector algorithm for approximate string
matching based on dynamic programming. *Journal of the ACM*, 46(3), 395-415.
https://doi.org/10.1145/316542.316550

Implementation: [`python/myers_bitvector.py`](python/myers_bitvector.py)
(arbitrary-precision word), [`cpp/myers_bitvector.cpp`](cpp/myers_bitvector.cpp)
(single 64-bit word, patterns ≤ 64 chars) — both validated against
brute-force Levenshtein distance over 2000+ random trials.

## 4. Wavefront Alignment (WFA)

### Overview
WFA (Marco-Sola et al., 2021) is an exact, provably-optimal gap-affine
alignment algorithm indexed by *score* rather than by sequence position.
For each candidate score s = 0, 1, 2, ..., it tracks, per diagonal
k = j - i, the furthest-reaching offset reachable with exactly that score —
a "wavefront." Matching characters along a diagonal are free ("greedy
extension"), so a wavefront snaps forward through long identical runs at no
cost. Because scores are explored in increasing order, the algorithm
terminates the instant a wavefront reaches the bottom-right corner,
guaranteeing the score found is optimal. This gives O(n*s + s^2) time,
where s is the optimal alignment score — near-linear whenever the sequences
are similar (small s), a fundamentally different complexity regime from
Needleman-Wunsch/Gotoh's O(mn), which is oblivious to how similar the
inputs actually are.

### Algorithm Steps
1. **Initialize** wavefront 0 at diagonal 0, extended greedily through any
   leading exact match.
2. **For each score s = 1, 2, ...**:
   - Compute the insertion wavefront (diagonal `k+1`, from a gap-open at
     score `s - (gap_open+gap_extend)` or a gap-extend at score
     `s - gap_extend`)
   - Compute the deletion wavefront (diagonal `k-1`, symmetric)
   - Compute the match/mismatch wavefront (substitution at
     `s - mismatch`, or landing from an insertion/deletion), then greedily
     extend each diagonal through further exact matches
   - Discard any state whose diagonal or offset falls outside the valid
     range `[-m, n]` / `[0, m]` / `[0, n]` (it has already overshot one
     sequence and can never be part of an optimal path)
   - If the wavefront for diagonal `n - m` reaches offset `m` (i.e., cell
     `(m, n)`), return `s` as the optimal score.

### Parameters
- `mismatch`, `gap_open`, `gap_extend`: non-negative costs (this
  implementation follows the original paper's MINIMIZATION convention,
  the opposite sign convention from the maximization scores used
  elsewhere in this repository)

### Time and Space Complexity
- Time: O(n*s + s^2) — near-linear for similar sequences (small s)
- Space: O(s^2) in this straightforward implementation (the "BiWFA"
  variant in WFA2-lib reduces this to O(s) via Hirschberg-style
  divide-and-conquer — see section 5's linear-space discussion)

### Use Cases
- Exact alignment of long, highly similar sequences: long-read-to-
  reference alignment, contig-to-contig comparison, assembly polishing
- Any setting where classical O(mn) DP is too slow but an approximate
  (non-exact) heuristic is undesirable

### Reference
Marco-Sola, S., Moure, J. C., Moreto, M., & Espinosa, A. (2021). Fast
gap-affine pairwise alignment using the wavefront algorithm.
*Bioinformatics*, 37(4), 456-463. https://doi.org/10.1093/bioinformatics/btaa777

Implementation: [`python/wavefront_alignment.py`](python/wavefront_alignment.py),
[`cpp/wavefront_alignment.cpp`](cpp/wavefront_alignment.cpp) (validated
against an independent affine-gap DP minimizer over 2000 random trials).

## 5. Minimizer Sketching + Co-linear Chaining (Minimap2-style)

### Overview
Modern long-read aligners do not index every k-mer of the reference: they
index only a sparse, deterministic subset called minimizers (Roberts et al.,
2004), then find matching seed pairs and chain them — a sparse dynamic
program over the seeds themselves — to identify collinear runs consistent
with a single alignment. Only the small region(s) around good chains are
ever passed to base-level DP, which is what lets tools like minimap2 (Li,
2018) scale to whole-genome references.

### Algorithm Steps
1. **Minimizer Sketch** - for every window of w consecutive k-mers, keep
   only the numerically smallest k-mer hash as that window's minimizer;
   any k-mer that is a minimizer for at least one window is indexed.
   This guarantees any shared substring of length >= w+k-1 between two
   sequences shares at least one minimizer, while indexing only ~2/(w+1)
   of all k-mers.
2. **Seed Lookup** - hash-table lookup of each query minimizer against the
   reference index, in O(1) expected time per minimizer.
3. **Co-linear Chaining** - find the highest-scoring subsequence of seeds
   with jointly increasing query and reference coordinates (bounded gap
   between consecutive seeds), via sparse dynamic programming — O(N log N)
   with a Fenwick-tree-backed implementation (as in minimap2 itself); this
   repository's reference implementation uses an O(N^2) DP for clarity.

### Parameters
- `k`: k-mer length; `w`: window size (larger w = sparser index, faster,
  slightly less sensitive)
- `max_gap`: maximum allowed gap between consecutive chained seeds

### Time and Space Complexity
- Sketching: O(L) for a sequence of length L
- Chaining: O(N log N) for N seeds (production), O(N^2) (this repository's
  reference implementation)

### Use Cases
- Long-read-to-genome and genome-to-genome seeding (minimap2, HISAT2,
  GraphAligner all use variants of this architecture)
- Any setting where full FM-index-based exact search is too memory-hungry
  or too slow for the read lengths/error rates involved

### Reference
Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M., & Yorke, J. A. (2004).
Reducing storage requirements for biological sequence comparison.
*Bioinformatics*, 20(18), 3363-3369. https://doi.org/10.1093/bioinformatics/bth408

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
*Bioinformatics*, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191

Implementation: [`python/minimizer_chaining.py`](python/minimizer_chaining.py),
[`cpp/minimizer_chaining.cpp`](cpp/minimizer_chaining.cpp) (validated to
recover the true offset of an exact substring, including under injected
point mutations).

## 6. Strobemers + MinHash Identity Estimation

### Overview
Two complementary newest-generation ideas. **Strobemers** (Sahlin, 2021;
used in `strobealign`, Sahlin 2022) link several short "strobes" chosen by a
*content-dependent* (hash-minimizing) rule rather than a fixed offset, so
the resulting seed tends to reappear even when an indel falls between the
strobes — directly fixing the brittleness of fixed-offset k-mers/minimizers
under insertions and deletions. **MinHash** (Jain et al., 2018, MashMap;
building on Broder's MinHash and Ondov et al.'s Mash) estimates the Jaccard
similarity of two k-mer sets from a small, fixed-size sample of each set's
minimum hash values, turning an O(nm) alignment question into an
O(sketch-size) set comparison.

### Algorithm Steps (Randstrobes)
1. For each position i, the first strobe is `seq[i:i+strobe_len]`.
2. The second strobe is chosen from a downstream window
   `[i+strobe_len+w_min, i+strobe_len+w_max)` as whichever candidate
   minimizes `hash(strobe1) XOR hash(candidate)` — a link function that
   depends on the content of both strobes, not just their positions.
3. The pair `(strobe1_pos, strobe2_pos)` and its combined hash form the
   strobemer, indexed like a k-mer.

### Algorithm Steps (MinHash / Jaccard Estimation)
1. Hash every k-mer of a sequence; keep the `sketch_size` smallest distinct
   hash values (a "bottom-k" sketch).
2. To compare two sequences, merge their sketches, keep the smallest
   `sketch_size` values overall, and take the fraction present in *both*
   original sketches as the Jaccard estimate.
3. Convert to an estimated per-base substitution rate via the Mash
   distance formula: `D = -1/k * ln(2J / (1+J))`.

### Time and Space Complexity
- Strobemer sketch: O(L * (w_max - w_min))
- MinHash sketch: O(L log(sketch_size))
- Jaccard estimate from two sketches of size s: O(s)

### Use Cases
- Indel-robust short-read seeding (strobealign reports higher throughput
  than BWA-MEM2/minimap2 at comparable-or-better accuracy for reads >= 150nt)
- Fast pre-filtering / identity triage before committing to full alignment
  (MashMap, Mash) — e.g., clustering long reads against huge reference
  databases

### Reference
Sahlin, K. (2021). Effective sequence similarity detection with strobemers.
*Genome Research*, 31(11), 2080-2094. https://doi.org/10.1101/gr.275648.121

Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and
accurate read alignment. *Genome Biology*, 23, 260.
https://doi.org/10.1186/s13059-022-02831-7

Jain, C., Dilthey, A., Koren, S., Aluru, S., & Phillippy, A. M. (2018). A
fast approximate algorithm for mapping long reads to large reference
databases. *Journal of Computational Biology*, 25(7), 766-779.
https://doi.org/10.1089/cmb.2018.0036

Implementation: [`python/strobemer_mapping.py`](python/strobemer_mapping.py),
[`cpp/strobemer_mapping.cpp`](cpp/strobemer_mapping.cpp) (validated to
recover the true offset of an exact substring via strobemer matching, and
to give a Jaccard similarity of 1.0 for identical sequences).

## 7. Seed-and-Extend (K-mer Hashing)

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

## 8. Burrows-Wheeler Transform (BWT) + FM-Index

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
| NGS read alignment | BWT + FM-Index, or Minimizer+Chaining | Ultra-fast exact matching / robust long-read seeding |
| Protein comparison | Smith-Waterman (Affine) | Finds functional domains, realistic gaps |
| Gene comparison | Needleman-Wunsch (Affine) | Complete gene alignment, realistic indel modeling |
| Finding motifs | Smith-Waterman | Local pattern matching |
| SNP calling | BWT + FM-Index | Align millions of reads |
| Sequences with indels | Smith-Waterman/NW (Affine), or Strobemers | Better models clustered gaps / survives indels in seeds |
| Exact edit distance, short-medium patterns | Myers' Bit-Vector | O(n) per query, minimal memory |
| Long, highly similar sequences | Wavefront Alignment (WFA) | Provably near-linear exact alignment |
| Fast similarity/identity triage | MinHash (strobemer_mapping.py) | O(sketch) comparison, no alignment needed |

### By Sequence Properties

| Property | Best Algorithm |
|----------|---------------|
| Very long (>100 Mbp) | BWT + FM-Index |
| Long (10K-100K bp), similar | Wavefront Alignment (WFA) or Minimizer+Chaining |
| Long (10K-100K bp), divergent | Seed-and-Extend |
| Medium (1K-10K bp) | Smith-Waterman |
| Short (<1K bp) | Needleman-Wunsch or Smith-Waterman |
| Highly similar | BWT + FM-Index, Seed-and-Extend, or WFA |
| Divergent | Smith-Waterman |
| Indel-prone (short reads) | Strobemers |

## Performance Characteristics

### Typical Running Times (approximate)

For aligning two sequences of length n:

| Algorithm | n=100 | n=1,000 | n=10,000 | n=100,000 |
|-----------|-------|---------|----------|-----------|
| Needleman-Wunsch / Gotoh affine | <1ms | 10ms | 1s | 100s |
| Smith-Waterman (linear/affine) | <1ms | 10ms | 1s | 100s |
| Myers' Bit-Vector | <1ms | <1ms | ~1ms | ~10ms |
| Wavefront Alignment (WFA)* | <1ms | <1ms | ~1-10ms | ~10-100ms |
| Seed-and-Extend / Minimizer+Chaining | <1ms | 1ms | 10ms | 100ms |
| FM-Index (search) | <1ms | <1ms | <1ms | <1ms |

*Note: FM-Index construction is O(n log n), but search is O(m) independent
of reference size. WFA's running time scales with the optimal alignment
score s, not n directly — the figures above assume s stays small (highly
similar sequences); for divergent sequences WFA approaches O(n^2), same as
classical DP.*

## References and Further Reading

See [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md) for the complete, DOI-linked
bibliography of all 44 algorithms surveyed for this thesis (including SIMD
libraries, production aligners, GPU/hardware accelerators, and theoretical
lower bounds). Core references for the algorithms implemented in this
repository:

1. Needleman, S.B. & Wunsch, C.D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(70)90057-4

2. Gotoh, O. (1982). "An improved algorithm for matching biological sequences". Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(82)90398-9

3. Smith, T.F. & Waterman, M.S. (1981). "Identification of common molecular subsequences". Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(81)90087-5

4. Hirschberg, D.S. (1975). "A linear space algorithm for computing maximal common subsequences". Communications of the ACM. https://doi.org/10.1145/360825.360861

5. Myers, G. (1999). "A fast bit-vector algorithm for approximate string matching based on dynamic programming". Journal of the ACM. https://doi.org/10.1145/316542.316550

6. Marco-Sola, S. et al. (2021). "Fast gap-affine pairwise alignment using the wavefront algorithm". Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa777

7. Roberts, M. et al. (2004). "Reducing storage requirements for biological sequence comparison". Bioinformatics. https://doi.org/10.1093/bioinformatics/bth408

8. Li, H. (2018). "Minimap2: pairwise alignment for nucleotide sequences". Bioinformatics. https://doi.org/10.1093/bioinformatics/bty191

9. Sahlin, K. (2021/2022). "Effective sequence similarity detection with strobemers" (Genome Research) / "Strobealign" (Genome Biology). https://doi.org/10.1101/gr.275648.121 / https://doi.org/10.1186/s13059-022-02831-7

10. Jain, C. et al. (2018). "A fast approximate algorithm for mapping long reads to large reference databases". Journal of Computational Biology. https://doi.org/10.1089/cmb.2018.0036

11. Altschul, S.F. et al. (1990). "Basic local alignment search tool (BLAST)". Journal of Molecular Biology. https://doi.org/10.1016/S0022-2836(05)80360-2

12. Burrows, M. & Wheeler, D.J. (1994). "A block-sorting lossless data compression algorithm". Technical Report 124, Digital Equipment Corporation.

13. Ferragina, P. & Manzini, G. (2000). "Opportunistic data structures with applications". Proceedings of FOCS. https://doi.org/10.1109/SFCS.2000.892127

14. Li, H. & Durbin, R. (2009). "Fast and accurate short read alignment with Burrows-Wheeler transform". Bioinformatics. https://doi.org/10.1093/bioinformatics/btp324

15. Langmead, B. et al. (2009). "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome". Genome Biology. https://doi.org/10.1186/gb-2009-10-3-r25
