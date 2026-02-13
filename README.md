# Bioinformatics Algorithms

Python and C++ implementations of common bioinformatics algorithms for sequence alignment and genome analysis.

## Overview

This repository contains implementations of algorithms used for alignment of sequencing reads (FASTQ) to a reference genome (FASTA):

1. **Dynamic Programming (DP)** – Exact alignment algorithms
   - Needleman-Wunsch (global alignment)
   - Smith-Waterman (local alignment)

2. **Seed-and-Extend with Hash Tables / k-mers** – Fast approximate alignment
   - Used in BLAST, MAQ, SOAP
   - K-mer indexing with hash tables
   - Seed extension with scoring

3. **Burrows-Wheeler Transform (BWT) + FM-index** – Ultra-fast exact matching
   - Foundational for BWA, Bowtie, HISAT2
   - Exact backward search via FM-index
   - Memory-efficient pattern matching

## Directory Structure

```
algorithms/
├── python/               # Python implementations
│   ├── needleman_wunsch.py
│   ├── smith_waterman.py
│   ├── seed_and_extend.py
│   └── bwt_fm_index.py
├── cpp/                  # C++ implementations
│   ├── needleman_wunsch.cpp
│   ├── smith_waterman.cpp
│   ├── seed_and_extend.cpp
│   └── bwt_fm_index.cpp
└── examples/             # Example usage
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

### 3. Seed-and-Extend (K-mer Hashing)

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

### 4. BWT + FM-Index (Exact Pattern Matching)

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
python3 python/smith_waterman.py
python3 python/seed_and_extend.py
python3 python/bwt_fm_index.py
```

### C++
```bash
# Compile and run
g++ -std=c++17 -o nw cpp/needleman_wunsch.cpp && ./nw
g++ -std=c++17 -o sw cpp/smith_waterman.cpp && ./sw
g++ -std=c++17 -o sae cpp/seed_and_extend.cpp && ./sae
g++ -std=c++17 -o bwt cpp/bwt_fm_index.cpp && ./bwt
```

## Algorithm Comparison

| Algorithm | Type | Speed | Memory | Accuracy | Use Case |
|-----------|------|-------|--------|----------|----------|
| Needleman-Wunsch | DP | Slow | High | Optimal | Small sequences, global alignment |
| Smith-Waterman | DP | Slow | High | Optimal | Small sequences, local alignment |
| Seed-and-Extend | Heuristic | Fast | Medium | Approximate | Medium sequences, BLAST-like |
| BWT + FM-index | Exact | Very Fast | Low | Exact | Large genomes, read alignment |

## Applications

- **Read Alignment:** Align millions of short reads from sequencing to a reference genome
- **Homology Search:** Find similar sequences in databases
- **Variant Calling:** Identify differences between sequences and reference
- **Assembly:** Overlap detection for sequence assembly

## References

- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins.
- Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences.
- Altschul, S. F., et al. (1990). Basic local alignment search tool (BLAST).
- Burrows, M., & Wheeler, D. J. (1994). A block-sorting lossless data compression algorithm.
- Ferragina, P., & Manzini, G. (2000). Opportunistic data structures with applications.

## License

See LICENSE file for details.
