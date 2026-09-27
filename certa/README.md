# CERTA: certified short-read fast path (prototype)

Prototype of tiers **S0/S1** from [RESEARCH_ROADMAP.md §3.3](../RESEARCH_ROADMAP.md#33-short-read-engine-illumina-and-other-sbs-platforms).
It maps each read with a **correctness certificate** when it can, and
passes every other read to a full aligner (minibwa, BWA-MEM2, …). One C++17
core (`include/certa/core.h`) runs unchanged on the CPU (x86-64, ARM64) and
under CUDA (A100 and other GPUs), so both back-ends give byte-identical
output.

> **Status:** research prototype for the feasibility aim (Aim 1). It is
> not a clinical tool. Single-end mapping only; paired-end logic, MAPQ
> calibration and a validated merge with the fallback aligner are next steps.

## What the certificate guarantees

For a read of length L and edit budget k, the read is split into k+1 parts.
Under unit-cost edit distance, every locus within k edits contains one part
*exactly* (pigeonhole). The index stores the q-mer at every s-th reference
position, and every part is looked up at all s shifts, so each exact
occurrence of a part is found (parts need length ≥ q+s−1). Parts with
more than `--cap` hits are not enumerated. With c such parts, the
guarantee shrinks to the **certified radius R = k − c**; it is never
silently weakened. Every candidate cluster is then verified by exact banded
edit distance.

For a certified read CERTA reports:

- `NM` (= d1): the **minimum edit distance over the whole reference, both
  strands**, together with an alignment that achieves it.
- `XB`: how many distinct loci tie at d1 (ties are broken deterministically
  from the read name, never by thread timing).
- `XD`: the exact second-best distance if it is ≤ R, otherwise −1, meaning no
  other locus within R.
- `XR` = R and `XT:Z:S0` (exact match) or `S1` (1..R edits).

`tests/test_certa.cpp` checks this against a brute-force oracle: Sellers'
DP over the entire reference, on a genome with diverged repeat copies, an
exact duplicate, a tandem repeat and N runs. For every read and every
(k, cap) setting it checks:

1. **Soundness:** the CIGAR re-scores to NM.
2. **Optimality:** NM equals the oracle's global minimum.
3. **Losslessness:** no read is rejected as "not found" while a locus
   within R exists.
4. **Uniqueness:** "no second locus" is never wrong.
5. **Determinism:** results are identical across thread counts.

Mutation testing shows the suite fails when the sampling shifts or the
diagonal band are broken.

What it does **not** guarantee: that the minimum-edit locus is the read's
*true* origin. A read from a diverged repeat copy can align better to
another copy, which is true of every aligner. The certificate is about
the optimization problem; MAPQ reflects how much margin there is.

## Build

Requirements: CMake ≥ 3.18, a C++17 compiler (GCC ≥ 9 or Clang ≥ 10), zlib
(optional, for `.gz` input), and CUDA ≥ 11.x for the GPU back-end.

```bash
cd certa
scripts/build.sh              # auto-detects CPU arch and nvcc; runs the certificate test
```

`build.sh` puts the build in `build-<host>-<arch>/`, so galaxy and a2 can
share one checkout. It builds the GPU back-end when `nvcc` is on `PATH` or
under `/usr/local/cuda*`, and it compiles for `sm_80` (A100) and `sm_86` by
default. Use `CUDA_ARCHS=80 scripts/build.sh` for A100 only, or
`CERTA_CUDA=OFF` to force a CPU-only build. `-march=native` (or
`-mcpu=native` on older aarch64 compilers) is used unless
`-DCERTA_NATIVE=OFF`.

### On galaxy and a2

```bash
ssh galaxy   # or a2
git clone https://github.com/d0bromir/algorithms.git && cd algorithms/certa
scripts/build.sh
scripts/demo.sh 50000000 2000000   # synthetic 50 Mbp genome, 2 M reads; CPU, plus GPU if present
```

On a GPU host, `demo.sh` runs the CPU and GPU back-ends on the same reads
and fails if their outputs differ.

## Run on real data

```bash
# Index once. GRCh38 no-alt analysis set: ~390 M entries; allow ~20 GB RAM
# while building, ~9 GB index file.
build-*/certa index GRCh38_no_alt_analysis_set.fa -o grch38.cidx -t 32

# Map. The CPU uses all cores by default; add --gpu [--device N] for CUDA.
build-*/certa map grch38.cidx sample_R1.fastq.gz -k 2 --gpu \
    -o certified.sam -u uncertified.fq --stats stats.json

# Fallback for the uncertified reads (example):
minibwa ... uncertified.fq > fallback.sam      # or bwa-mem2 mem ...
```

`stats.json` records the certified fraction by tier, why the other reads
were not certified (radius < 0 = repetitive, not found = more than R edits
or no locus, …), and a timing breakdown (index load, GPU upload,
read/encode, map, output). The **certified fraction on real NovaSeq X
data is the Aim 1 measurement**. It bounds the speedup the fast path can
deliver (Amdahl's law).

Parameters:

| option | default | meaning |
|---|---|---|
| `-q` (index) | 22 | q-mer length, 8–32 |
| `-s` (index) | 8 | sampling step, 1–16 (memory ∝ 1/s) |
| `-k` (map) | 2 | edit budget, 0–5; reads need length ≥ (k+1)(q+s−1) (87 bp for k=2 with the defaults; 145 bp for k=4) |
| `--cap` | 32 | max hits enumerated per part (≤ 32) |
| `--batch` | 200 000 CPU / 1 000 000 GPU | reads per batch |

## Results so far (local development machine only)

| Setting | Result |
|---|---|
| Certificate test (oracle), 10 seeds × 480 reads × 5 (k, cap) configs, x86-64 | all checks pass |
| Mutation tests (drop a sampling shift; narrow the band) | caught by the test |
| Synthetic 20 Mbp genome (repeats, SDs, STRs, N gaps), 500 k × 150 bp reads, 0.3 % sub + 0.01 % indel, k=2 | 95.7 % certified (60.8 % S0 exact, 34.9 % S1); 2.5 % repetitive (radius < 0), 1.7 % > R edits |
| Placement vs. simulated origin, certified reads | MAPQ 60: 100 %, MAPQ 40: 100 %, MAPQ 20: 99.94 %, MAPQ 0 (ties): 43 % |
| CPU, 8 threads (laptop, WSL2) | 357 k reads/s (map step) |
| GPU, RTX 3050 Laptop (sm_86) | 1.57 M reads/s (map step, including transfers); output byte-identical to CPU |
| Kernel resources (sm_80) | 48 registers, 27 KB stack/thread (≈6 GB local-memory reservation on an A100) |

These are synthetic-data numbers on a laptop. They show that the mechanics
work and that CPU and GPU give identical output. They are **not**
evidence about real genomes, where repeat content, error profiles and
the variant spectrum differ. Not yet run: galaxy/a2, ARM64 hardware (CI
covers it via GitHub's ARM runners), A100, GRCh38, real FASTQ.

## Layout

```
include/certa/core.h   per-read algorithm (host + device)
src/index.*            FASTA loading, sampled q-mer index, save/load
src/seqio.*            FASTQ/FASTA reading (plain or gzip)
src/mapper_cpu.cpp     multi-threaded CPU driver
src/mapper_gpu.cu      CUDA driver (one thread per read)
src/main.cpp           CLI, SAM/FASTQ output, statistics
tests/test_certa.cpp   brute-force certificate test
tools/simulate.py      synthetic genome/reads and placement evaluation
scripts/               build.sh, demo.sh
```

## Known limitations and next steps

- **Single-end only.** Add mate-pair consistency and rescue.
- **MAPQ:** a provisional 20 × (distance gap) rule, capped at 60. Calibrate
  it on simulated reads from HPRC assemblies (roadmap §2.4).
- **Reference storage:** 1 byte per base (3.2 GB for GRCh38). Packing to
  2 bits and memory-mapping the index would cut load time and memory.
- **GPU:** one thread per read, with local-memory workspaces. Warp-cooperative
  lookups, overlapping transfers with compute, and a device-side output
  path are the obvious optimizations.
- **Clusters wider than 48 diagonals** (long tandem repeats) are left
  uncertified instead of being split.
- **Fallback merge:** combine the certified SAM with the fallback
  aligner's output (e.g., `samtools merge` + sort) in one fused pipeline
  (roadmap §2.6).
