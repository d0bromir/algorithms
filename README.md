# Bioinformatics Algorithms

Python and C++ implementations of common bioinformatics algorithms for sequence alignment and genome analysis.

These are educational primitives, not validated clinical read mappers. The
literature review is selective, not a demonstrated complete survey. No
end-to-end benchmark or clinical superiority result is established here.

## Overview

This repository contains implementations of algorithms used for alignment of sequencing reads (FASTQ) to a reference genome (FASTA):

1. **Dynamic Programming (DP)** – Exact alignment algorithms
   - Needleman-Wunsch (global alignment)
   - Needleman-Wunsch with Affine Gap Penalties (Gotoh's algorithm)
   - Smith-Waterman (local alignment)
   - Smith-Waterman with Affine Gap Penalties (biologically realistic gaps)
   - Hirschberg's algorithm (linear-space global alignment)

2. **Bit-Parallel / Score-Parameterized Exact Alignment**
   - Myers' bit-vector algorithm (O(n·⌈m/w⌉) edit distance)
   - Wavefront Alignment (WFA) – O(n·s + s²) exact gap-affine alignment

3. **Seed-and-Extend with Hash Tables / k-mers / Sketches** – Fast approximate alignment
   - Classic BLAST-style k-mer seed-and-extend (used in BLAST, MAQ, SOAP)
   - Minimizer sketching + co-linear chaining (minimap2-style)
   - Strobemers + MinHash identity estimation (strobealign/MashMap-style)

4. **Burrows-Wheeler Transform (BWT) + FM-index** – Ultra-fast exact matching
   - Foundational for BWA, Bowtie, HISAT2
   - Exact backward search via FM-index
   - Memory-efficient pattern matching

A selective literature review and SWOT analysis of these and additional
algorithms (SIMD libraries, production aligners, GPU/hardware accelerators,
and theoretical lower bounds) is in [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md).

## Directory Structure

```
algorithms/
├── python/               # Python implementations
│   ├── needleman_wunsch.py
│   ├── needleman_wunsch_affine.py
│   ├── smith_waterman.py
│   ├── smith_waterman_affine.py
│   ├── seed_and_extend.py
│   ├── bwt_fm_index.py
│   ├── hirschberg.py
│   ├── myers_bitvector.py
│   ├── wavefront_alignment.py
│   ├── minimizer_chaining.py
│   └── strobemer_mapping.py
├── cpp/                  # C++ implementations
│   ├── needleman_wunsch.cpp
│   ├── needleman_wunsch_affine.cpp
│   ├── smith_waterman.cpp
│   ├── smith_waterman_affine.cpp
│   ├── seed_and_extend.cpp
│   ├── bwt_fm_index.cpp
│   ├── hirschberg.cpp
│   ├── myers_bitvector.cpp
│   ├── wavefront_alignment.cpp
│   ├── minimizer_chaining.cpp
│   └── strobemer_mapping.cpp
├── examples/             # Example usage
└── SWOT_ANALYSIS.md      # Literature review + SWOT analysis for the PhD thesis
```

## Algorithms

### 1. Needleman-Wunsch (Global Alignment)

Dynamic programming algorithm for finding the optimal global alignment between two sequences.

- **Time Complexity:** O(m × n)
- **Space Complexity:** O(m × n)
- **Use Case:** Aligning complete sequences

**Python Usage:**
```python
from needleman_wunsch import needleman_wunsch

seq1 = "GATTACA"
seq2 = "GCATGCU"
aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o nw cpp/needleman_wunsch.cpp
./nw
```

### 2. Smith-Waterman (Local Alignment)

Dynamic programming algorithm for finding the optimal local alignment between two sequences.

- **Time Complexity:** O(m × n)
- **Space Complexity:** O(m × n)
- **Use Case:** Finding similar regions in sequences

**Python Usage:**
```python
from smith_waterman import smith_waterman

seq1 = "GGTTGACTA"
seq2 = "TGTTACGG"
aligned1, aligned2, score = smith_waterman(seq1, seq2)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o sw cpp/smith_waterman.cpp
./sw
```

#### 2.1 Smith-Waterman with Affine Gap Penalties

Enhanced version using affine gap penalties (gap_open + gap_extend × length) for more biologically realistic alignments.

- **Time Complexity:** O(m × n)
- **Space Complexity:** O(m × n) - uses three matrices
- **Use Case:** Protein alignment, sequences with indels

**Python Usage:**
```python
from smith_waterman_affine import smith_waterman_affine

seq1 = "GGTTGACTA"
seq2 = "TGTTACGG"
aligned1, aligned2, score = smith_waterman_affine(seq1, seq2, 
                                                   gap_open=-3, 
                                                   gap_extend=-1)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o sw_affine cpp/smith_waterman_affine.cpp
./sw_affine
```

#### 2.2 Needleman-Wunsch with Affine Gap Penalties (Gotoh's Algorithm)

The global-alignment counterpart to Smith-Waterman-affine: same three-matrix
(match/insert/delete) recurrence, but without local-alignment flooring, so
gaps are penalized `gap_open + k * gap_extend` while the whole sequences are
still aligned end-to-end.

- **Time Complexity:** O(m × n)
- **Space Complexity:** O(m × n)
- **Use Case:** Global alignment where indels cluster (e.g., structural
  variants), matching [Gotoh (1982)](https://doi.org/10.1016/0022-2836(82)90398-9)

**Python Usage:**
```python
from needleman_wunsch_affine import needleman_wunsch_affine

seq1 = "GGTTGACTA"
seq2 = "TGTTACGG"
aligned1, aligned2, score = needleman_wunsch_affine(seq1, seq2,
                                                     gap_open=-3,
                                                     gap_extend=-1)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o nw_affine cpp/needleman_wunsch_affine.cpp
./nw_affine
```

### 3. Myers' Bit-Vector Algorithm (Edit Distance)

Packs an entire DP column into machine words so a whole column update
becomes a handful of bitwise operations, giving O(n) time for patterns up to
the machine word size. Underlies the `edlib` library; BWA-MEM uses affine-gap
DP, while GenASM is based on modified Bitap, not this Myers recurrence.

- **Time Complexity:** O(n · ⌈m/w⌉) — O(n) for patterns ≤ word size w
- **Space Complexity:** O(σ) for the character bitmask table, O(1) state
- **Use Case:** Fast exact edit distance / semi-global search, matching
  [Myers (1999)](https://doi.org/10.1145/316542.316550)

**Python Usage:**
```python
from myers_bitvector import myers_bit_vector_edit_distance

dist, end_pos = myers_bit_vector_edit_distance("GATTACA", "GACTATA")
```

**C++ Usage:**
```bash
g++ -std=c++17 -o myers cpp/myers_bitvector.cpp
./myers
```

### 4. Wavefront Alignment (WFA)

A score-parameterized exact algorithm for gap-affine global alignment:
O(n·s + s²) time, where `s` is the optimal alignment cost with fixed
penalties. This is linear in length for bounded `s`, but remains quadratic
when `s` grows proportionally to length, even at a low fixed error rate.
Exactness requires compatible scoring/boundaries and no heuristic pruning;
it does not guarantee the correct genomic mapping location.

- **Time Complexity:** O(n·s + s²), plus input scanning at score zero
- **Space Complexity:** O(s²) in this simplified version (O(s) in the
  BiWFA variant of the original library)
- **Use Case:** Exact alignment of long, highly similar sequences (long
  reads, contigs), matching
  [Marco-Sola et al. (2021)](https://doi.org/10.1093/bioinformatics/btaa777)

**Python Usage:**
```python
from wavefront_alignment import wavefront_alignment

score = wavefront_alignment("GATTACA", "GATCACA", mismatch=4, gap_open=6, gap_extend=2)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o wfa cpp/wavefront_alignment.cpp
./wfa
```

### 5. Minimizer Sketching + Co-linear Chaining (Minimap2-style)

Indexes only a sparse, deterministic subset of k-mers (minimizers), then
chains matching seeds with a sparse O(N log N) dynamic program to find
candidate alignment regions — the architecture behind minimap2, minigraph,
GraphAligner and vg Giraffe.

- **Time Complexity:** O(L) sketching; chaining O(N log N) exact
  (RMQ-based) or O(N·h) heuristic as in minimap2's default (O(N²) in this
  simplified reference implementation)
- **Space Complexity:** O(L / w) for the minimizer index
- **Use Case:** Long-read / whole-genome seeding, matching
  [Li (2018)](https://doi.org/10.1093/bioinformatics/bty191)

**Python Usage:**
```python
from minimizer_chaining import build_minimizer_index, find_seed_matches, chain_seeds

ref_index = build_minimizer_index(reference, k=15, w=10)
seeds = find_seed_matches(query, ref_index, k=15, w=10)
chain, score = chain_seeds(seeds, k=15)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o minchain cpp/minimizer_chaining.cpp
./minchain
```

### 6. Strobemers + MinHash Identity Estimation

Two newest-generation ideas: strobemers (content-linked seeds that survive
indels far better than fixed k-mers, used in `strobealign`) and MinHash
sketching for O(sketch-size) Jaccard/identity estimation without alignment
(used in `MashMap`/`Mash`).

- **Time Complexity:** O(L) strobemer sketch; O(L log(sketch size)) MinHash sketch
- **Space Complexity:** O(sketch size)
- **Use Case:** Indel-robust short-read seeding and fast similarity
  triage, matching [Sahlin (2021/2022)](https://doi.org/10.1186/s13059-022-02831-7)
  and [Jain et al. (2018)](https://doi.org/10.1089/cmb.2018.0036)

**Python Usage:**
```python
from strobemer_mapping import strobemer_index, find_strobemer_matches, minhash_sketch, estimate_jaccard

idx = strobemer_index(reference, strobe_len=6, w_min=6, w_max=20)
matches = find_strobemer_matches(query, idx, strobe_len=6, w_min=6, w_max=20)

sketch_a = minhash_sketch(seq_a, k=12, sketch_size=200)
sketch_b = minhash_sketch(seq_b, k=12, sketch_size=200)
similarity = estimate_jaccard(sketch_a, sketch_b)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o strobe cpp/strobemer_mapping.cpp
./strobe
```

### 7. Seed-and-Extend (K-mer Hashing)

Fast approximate alignment using exact k-mer matches followed by extension.

- **Time Complexity:** O(n) for indexing, O(m) for query
- **Space Complexity:** O(n) for index
- **Use Case:** Fast similarity search (BLAST-like)

**Python Usage:**
```python
from seed_and_extend import seed_and_extend

reference = "ACGTACGTACGTAAACCCGGGTTTACGTACGT"
query = "ACGTAAACCCGGG"
alignments = seed_and_extend(reference, query, k=5)
```

**C++ Usage:**
```bash
g++ -std=c++17 -o sae cpp/seed_and_extend.cpp
./sae
```

### 8. BWT + FM-Index (Exact Pattern Matching)

Burrows-Wheeler Transform with FM-index for memory-efficient exact pattern matching.

- **Time Complexity:** O(m) for pattern of length m
- **Space Complexity:** O(n) compressed
- **Use Case:** Fast read alignment (BWA, Bowtie)

**Python Usage:**
```python
from bwt_fm_index import FMIndex

text = "ACGTACGTACGT"
fm_index = FMIndex(text)
result = fm_index.search("ACG")
print(f"Count: {result['count']}, Positions: {result['positions']}")
```

**C++ Usage:**
```bash
g++ -std=c++17 -o bwt cpp/bwt_fm_index.cpp
./bwt
```

## Running the Examples

### Python
```bash
# Run individual algorithms
python3 python/needleman_wunsch.py
python3 python/needleman_wunsch_affine.py
python3 python/smith_waterman.py
python3 python/smith_waterman_affine.py
python3 python/hirschberg.py
python3 python/myers_bitvector.py
python3 python/wavefront_alignment.py
python3 python/minimizer_chaining.py
python3 python/strobemer_mapping.py
python3 python/seed_and_extend.py
python3 python/bwt_fm_index.py
```

### C++
```bash
# Compile and run
g++ -std=c++17 -o nw cpp/needleman_wunsch.cpp && ./nw
g++ -std=c++17 -o nw_affine cpp/needleman_wunsch_affine.cpp && ./nw_affine
g++ -std=c++17 -o sw cpp/smith_waterman.cpp && ./sw
g++ -std=c++17 -o sw_affine cpp/smith_waterman_affine.cpp && ./sw_affine
g++ -std=c++17 -o hirsch cpp/hirschberg.cpp && ./hirsch
g++ -std=c++17 -o myers cpp/myers_bitvector.cpp && ./myers
g++ -std=c++17 -o wfa cpp/wavefront_alignment.cpp && ./wfa
g++ -std=c++17 -o minchain cpp/minimizer_chaining.cpp && ./minchain
g++ -std=c++17 -o strobe cpp/strobemer_mapping.cpp && ./strobe
g++ -std=c++17 -o sae cpp/seed_and_extend.cpp && ./sae
g++ -std=c++17 -o bwt cpp/bwt_fm_index.cpp && ./bwt
```

## Algorithm Comparison

| Algorithm | Type | Speed | Memory | Accuracy | Use Case |
|-----------|------|-------|--------|----------|----------|
| Needleman-Wunsch | DP | Slow | High | Optimal | Small sequences, global alignment |
| Needleman-Wunsch (Affine, Gotoh) | DP | Slow | High | Optimal | Global alignment, clustered indels |
| Smith-Waterman | DP | Slow | High | Optimal | Small sequences, local alignment |
| Smith-Waterman (Affine) | DP | Slow | High | Optimal | Proteins, realistic gap modeling |
| Hirschberg | DP (linear-space) | Slow | Low | Optimal | Long sequences, memory-limited global alignment |
| Myers' Bit-Vector | Bit-parallel exact | Very Fast | Low | Exact | Short-to-medium patterns, edit distance |
| Wavefront Alignment (WFA) | Score-parameterized exact | Fast* | Medium | Optimal | Long, similar sequences (*fast when score s is small) |
| Minimizer + Chaining | Sketching + heuristic | Very Fast | Low | Approximate | Long-read / whole-genome seeding |
| Strobemers + MinHash | Sketching + heuristic | Very Fast | Low | Approximate | Indel-robust seeding, identity estimation |
| Seed-and-Extend | Heuristic | Fast | Medium | Approximate | Medium sequences, BLAST-like |
| BWT + FM-index | Exact | Very Fast | Low | Exact | Large genomes, read alignment |

A full literature review with Strengths/Weaknesses/Opportunities/Threats for
each of these plus ~50 additional algorithms and tools (SIMD libraries,
production aligners such as Bowtie/BWA-MEM2/minibwa/strobealign/minimap2/vg,
clinical pipelines such as DRAGEN/Parabricks, GPU aligners, and hardware
accelerators) is in [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md). Future research
directions, a proposed adaptive aligner, and a clinical benchmarking
protocol are in [RESEARCH_ROADMAP.md](RESEARCH_ROADMAP.md). A prototype of
that aligner's certified short-read fast path (C++17 for x86-64/ARM64, plus
CUDA) is in [certa/](certa/README.md). The proposal for an exact, provably
optimal aligner that builds on it is in
[EXACT_ALIGNMENT_PROPOSAL.md](EXACT_ALIGNMENT_PROPOSAL.md).

## Applications

- **Read Alignment:** Align millions of short reads from sequencing to a reference genome
- **Homology Search:** Find similar sequences in databases
- **Variant Calling:** Identify differences between sequences and reference
- **Assembly:** Overlap detection for sequence assembly

## References

See [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md) for the complete, DOI-linked
bibliography (60+ algorithms/papers, DOIs verified Sept. 2026). Core references for algorithms
implemented in this repository:

- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. https://doi.org/10.1016/0022-2836(70)90057-4
- Gotoh, O. (1982). An improved algorithm for matching biological sequences. https://doi.org/10.1016/0022-2836(82)90398-9
- Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. https://doi.org/10.1016/0022-2836(81)90087-5
- Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common subsequences. https://doi.org/10.1145/360825.360861
- Myers, G. (1999). A fast bit-vector algorithm for approximate string matching based on dynamic programming. https://doi.org/10.1145/316542.316550
- Marco-Sola, S., et al. (2021). Fast gap-affine pairwise alignment using the wavefront algorithm. https://doi.org/10.1093/bioinformatics/btaa777
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. https://doi.org/10.1093/bioinformatics/bty191
- Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and accurate read alignment. https://doi.org/10.1186/s13059-022-02831-7
- Jain, C., et al. (2018). A fast approximate algorithm for mapping long reads to large reference databases. https://doi.org/10.1089/cmb.2018.0036
- Altschul, S. F., et al. (1990). Basic local alignment search tool (BLAST). https://doi.org/10.1016/S0022-2836(05)80360-2
- Burrows, M., & Wheeler, D. J. (1994). A block-sorting lossless data compression algorithm.
- Ferragina, P., & Manzini, G. (2000). Opportunistic data structures with applications. https://doi.org/10.1109/SFCS.2000.892127

## License

See LICENSE file for details.
