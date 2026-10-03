# CERTA benchmark: GIAB HG002, Illumina NovaSeq X (September 2026)

Single-end benchmark of the CERTA prototype against the fastest open-source
short-read mappers on real clinical-grade data. Raw data are in
[`results/hg002_2026-09/`](results/hg002_2026-09/). The scripts are in this
directory (`subsample.sh`, `build_indexes.sh`, `run_bench.sh`,
`accuracy.sh`, `compare_sam.py`).

## Headline

1. **Certified fraction:** CERTA certifies **83.4 % of all 444.5 M HG002
   NovaSeq X reads** (k = 2): 68.7 % exact and 14.8 % with 1–2 edits. The
   roadmap's Aim 1 threshold was 60 %.
2. **The certificate holds on real data.** 1,672,120 certified reads have
   an unclipped BWA-MEM2 alignment, and BWA-MEM2 never finds fewer edits
   except in 3 cases, all explained by a single `N` in GRCh38.
   (BWA indexes replace `N` with a random base; CERTA scores `N` as a
   mismatch.)
3. **Complete mapping** (certified reads plus a fallback for the rest):
   only the concurrent GPU + CPU design beats the fastest tool. On galaxy,
   CERTA on the A100, piping uncertified reads into minibwa on 64 CPU
   threads, maps all reads in **24.4 s vs 26.5 s for minibwa (1.08×)**,
   using 23 % less CPU time. On CPU only (a2, x86-64), CERTA + fallback
   is **~1.5× slower** than minibwa.
4. **Why the gain is small:** the **16.6 % of reads CERTA cannot certify
   account for ~70 % of minibwa's CPU time.** The reads that are cheap to
   certify are the reads minibwa already maps cheaply. Aim 1 should
   therefore be judged by the *share of competitor cost on certifiable
   reads*, not by the share of reads (see §5).

## 1. Setup

| | galaxy | a2 |
|---|---|---|
| CPU | Ampere Altra (Neoverse-N1), 128 cores, 1 NUMA node | 4 × Intel Xeon Gold 5218 (Cascade Lake, AVX-512), 64 cores, 4 NUMA nodes |
| RAM | 246 GB | 376 GB |
| GPU | 2 × NVIDIA A100 80GB PCIe (GPU 1 used; ollama servers resident but idle) | none |
| OS / compiler | Ubuntu 25.10, GCC 15.2, CUDA 13.1 | Ubuntu 24.04, GCC 13.3 |

**Data**
- **Reads:** GIAB HG002, NovaSeq X 25B, run SRR37356339 (study
  PRJNA1427896, Weill Cornell), R1 only. 444,524,794 × 150 bp;
  md5-verified from ENA.
- **Samples:** a seeded (seed 42) uniform random sample of 20,001,073
  reads for timing, and its first 2,000,000 reads for accuracy.
- **Reference:** Ensembl GRCh38 primary assembly (194 sequences, no ALT
  contigs). It was already installed on galaxy.

**Tools** (Bioconda, versions from package metadata): minibwa 0.7,
strobealign 0.18.0, BWA-MEM2 2.3 (its binary reports 2.2.1), minimap2
2.31 (`-x sr`), bowtie2 2.5.5, samtools 1.24.
- All tools run single-end with default settings, 64 threads, reading the
  same uncompressed FASTQ.
- SAM goes to `/dev/null` for timing: formatting is timed, disk writes are
  not.
- Each tool gets a warm-up run (index in page cache), then 3 timed
  repeats. Tables show the median, with the range in brackets.

**CERTA modes**
- `certa`: certified reads only; the rest are written to FASTQ. This is
  **not** a complete mapping.
- `+fallback`: CERTA, then minibwa on the uncertified reads, sequentially.
- `-pipe`: CERTA streams uncertified reads through a pipe into minibwa,
  which runs concurrently. On the CPU the 64 threads are split 32 + 32; with
  the GPU, CERTA runs on the A100 and minibwa gets 64 threads.
- **Versions:** v0.1 is the initial prototype. v0.2 parses and encodes
  input on several threads. v0.3 also slices records with one copy. Output
  is byte-identical across versions and across CPU and GPU.

## 2. Timing (20,001,073 reads, 64 threads)

### galaxy (ARM64 + A100)

| Configuration | Complete? | Wall (s) | CPU (s) | Peak RSS (GiB) | vs minibwa |
|---|---|---|---|---|---|
| **CERTA v0.3 GPU → pipe → minibwa** | yes | **24.40** [24.27–24.41] | 1058 | 12.1 | **1.08×** |
| minibwa | yes | 26.46 [26.39–26.52] | 1368 | 9.1 | 1.00× |
| CERTA v0.3 GPU + fallback (sequential) | yes | 33.93 [33.91–34.36] | 1043 | 12.2 | 0.78× |
| CERTA v0.3 CPU + fallback (sequential) | yes | 35.16 [34.32–35.54] | 1687 | 12.2 | 0.75× |
| strobealign | yes | 38.59 [38.59–39.83] | 1090 | 13.4 | 0.69× |
| CERTA v0.3 CPU → pipe → minibwa (32 + 32 threads) | yes | 40.28 [40.23–40.30] | 1664 | 11.8 | 0.66× |
| minimap2 `-x sr` | yes | 68.80 [68.73–68.91] | 3117 | 11.9 | 0.38× |
| BWA-MEM2 | yes | 99.69 [99.62–100.10] | 5224 | 31.1 | 0.27× |
| bowtie2 | yes | 187.80 [187.76–187.82] | 12026 | 3.9 | 0.14× |
| *CERTA v0.3 GPU, certified only* | *83.4 % of reads* | *15.81* | *68* | *10.2* | — |
| *CERTA v0.3 CPU, certified only* | *83.4 % of reads* | *17.28* | *716* | *8.6* | — |

Earlier CERTA versions on galaxy (complete mapping, sequential fallback):
v0.1 46.2 s (CPU) / 46.5 s (GPU), v0.2 39.0 s / 38.5 s. The certified-only
step went from 27.3 s (v0.1) to 20.6 s (v0.2) to 17.3 s (v0.3) on the CPU,
all from input parsing. Mapping itself takes ~10.6 s (1.9 M reads/s).

### a2 (x86-64, CPU only)

| Configuration | Complete? | Wall (s) | CPU (s) | Peak RSS (GiB) | vs minibwa |
|---|---|---|---|---|---|
| minibwa | yes | 37.99 [37.81–38.25] | 1551 | 9.3 | 1.00× |
| strobealign | yes | 57.55 [56.52–58.05] | 997 | 13.4 | 0.66× |
| CERTA v0.3 CPU → pipe → minibwa (32 + 32) | yes | 57.61 [55.81–58.96] | 2158 | 11.7 | 0.66× |
| CERTA v0.3 CPU + fallback (sequential) | yes | 58.74 [52.09–64.63] | 2280 | 12.2 | 0.65× |
| minimap2 `-x sr` | yes | 91.96 [91.19–93.88] | 3737 | 12.9 | 0.41× |
| BWA-MEM2 | yes | 117.30 [115.21–118.25] | 5786 | 32.0 | 0.32× |
| bowtie2 | yes | 156.58 [156.48–157.77] | 9788 | 3.9 | 0.24× |
| *CERTA v0.3 CPU, certified only* | *83.4 % of reads* | *35.87* | *1425* | *8.6* | — |

a2's CERTA timings vary more (4 NUMA nodes; no thread pinning was used).

## 3. Full run: all 444.5 M reads (galaxy, A100, v0.1)

| | Reads | Share |
|---|---|---|
| Certified | 370,927,327 | 83.44 % |
| S0 (exact) | 305,256,410 | 68.67 % |
| S1 (1–2 edits) | 65,670,917 | 14.77 % |
| Ties at best distance (within certified) | 10,794,424 | 2.43 % |
| Uncertified: repetitive (radius < 0) | 34,589,072 | 7.78 % |
| Uncertified: > R edits or no locus | 38,814,803 | 8.73 % |
| Uncertified: too short / cluster too wide / crosses contig | 193,592 | 0.04 % |

- **GPU time:** 42 s of GPU mapping (10.6 M reads/s).
- **End-to-end time:** 14 min 8 s, almost all of it single-threaded
  gunzip of the 29.4 GB input (v0.1 parser; overlapped with the index
  builds). The 20 M-read sample gives the same fractions (83.44 %).

## 4. Accuracy (2 M reads; reference = BWA-MEM2)

There is no per-read truth for real reads. The table below measures
**agreement with BWA-MEM2** (same contig, unclipped start within 10 bp).
It is not correctness, and GIAB variant-level accuracy is future work.

| Tool | Mapped | Agrees with BWA-MEM2 | Agrees where BWA-MEM2 MAPQ ≥ 30 |
|---|---|---|---|
| minibwa | 99.74 % | 94.48 % | 99.96 % |
| strobealign | 99.92 % | 93.97 % | 99.86 % |
| minimap2 | 99.52 % | 93.95 % | 99.81 % |
| bowtie2 | 98.58 % | 93.60 % | 99.66 % |

CERTA's certified reads (83.66 % of the 2 M), by tier and MAPQ:

| Tier | MAPQ | Reads | Agrees with BWA-MEM2 |
|---|---|---|---|
| S0 | 60 | 1,114,511 | 100.000 % |
| S0 | 40 | 134,280 | 100.000 % |
| S0 | 20 | 108,064 | 99.999 % |
| S1 | 40 | 192,891 | 99.998 % |
| S1 | 20 | 75,176 | 99.944 % |
| S0 / S1 | 0 (ties) | 48,221 | ~35 % (ties broken differently) |

The accuracy reports are **byte-identical on galaxy (ARM64) and a2
(x86-64)**.

## 5. Interpretation for the thesis

- **The certificate is sound on real data.** No unexplained violation
  among 1.67 M checked reads.
- **Most reads are easy:** 83 % are certifiable at k = 2, and 69 % are
  exact matches.
- **But easy reads are already cheap for good mappers.** minibwa spends
  about 70 % of its CPU time on the hard 17 %. An infinitely fast
  certified path could therefore save at most about 30 % of minibwa's CPU
  time. So "fraction of reads certified" was the wrong go/no-go metric;
  use the **share of the best competitor's cost on certifiable reads**.
- **Where CERTA does help:**
  1. Heterogeneous hardware: the GPU certifies easy reads while all CPU
     cores map hard ones (1.08× wall and 23 % less CPU than minibwa here).
  2. Workloads dominated by exact or near-exact reads.
  3. As a provable-correctness layer (certified NM, exact second-best
     distance) rather than a speed layer.
- **Next steps, ordered by expected gain:**
  1. A faster fallback for the *hard* reads. That is where the time goes.
  2. A fused in-process fallback, without the pipe and FASTQ round-trip.
  3. Paired-end support. Mate information certifies more reads and is the
     clinical norm.
  4. GIAB variant-level accuracy (DeepVariant) of the hybrid output.
  5. NUMA-aware thread placement on multi-socket hosts.

## 5b. v0.4: adaptive certificate, and where a 10x speedup could come from

**What changed.** v0.4 splits each read into as many parts as fit (5 for
150 bp), enumerates the *rarest* parts that fit a per-strand hit budget,
and certifies radius |S| − 1. It also adds a certified-repeat tier (SR:
≥ 2 exact copies found ⇒ d1 = 0, MAPQ 0) and bounded cluster
verification. Raw data: `sweep_2M_v04_*.json`, `galaxy_timing_v04.tsv`,
`full_444M_gpu_v04.stats.json`.

**Certified fraction, 2 M reads**

| k / budget | Certified | S0 | S1 | SR | Uncertified: repetitive | Uncertified: > R |
|---|---|---|---|---|---|---|
| v0.3 (k=2, cap 32/part) | 83.66 % | — | — | — | 7.78 % | 8.73 % |
| 2 / 32 | 89.38 % | 72.38 % | 15.67 % | 1.32 % | 2.68 % | 7.94 % |
| 2 / 256 | 91.31 % | 73.68 % | 16.63 % | 1.00 % | 1.13 % | 7.55 % |
| 3 / 256 | 92.51 % | 73.68 % | 17.83 % | 1.00 % | 1.13 % | 6.33 % |
| 4 / 128 | 92.44 % | 73.30 % | 18.01 % | 1.12 % | 1.55 % | 5.98 % |
| **4 / 256 (default)** | **93.03 %** | 73.67 % | 18.36 % | 1.00 % | 1.13 % | 5.79 % |

On all 444.5 M reads, v0.4 certifies **92.89 %**.

**Timing, galaxy, 20 M reads (median of 3)**

| Configuration | Wall | vs minibwa | vs BWA-MEM2 |
|---|---|---|---|
| v0.4 A100 → pipe → minibwa | 24.40 s | 1.09× | 4.09× |
| v0.4 A100 + fallback (sequential) | 31.59 s | 0.84× | 3.16× |
| v0.4 CPU + fallback (sequential) | 44.16 s | 0.60× | 2.26× |
| *v0.4 A100, certified only (92.9 %)* | *17.93 s* | — | — |

v0.4 halves the fallback's input (1.42 M reads instead of 3.3 M), but the
GPU kernel became slower: 11.0 s of mapping vs about 2 s for v0.3. More
parts and larger budgets make GPU threads diverge. End-to-end time is
unchanged.

**Anatomy of the hard reads.** These are the 31.6 M reads v0.4 leaves
uncertified out of 444.5 M, classified from minibwa's alignment of them.

| Class | Share | minibwa CPU per read |
|---|---|---|
| Repeat, MAPQ 0 | 24.2 % | ~440 µs |
| Indel/clip event shared with ≥ 3 other hard reads | 25.7 % | 500–760 µs |
| Substitutions only, NM ≤ 2, unique (radius limited by repetitive parts) | 16.1 % | ~320 µs |
| Substitutions only, NM ≥ 3 | 19.3 % | 400–550 µs |
| Private indel/clip events, unmapped, other | ~15 % | 150–1060 µs |

minibwa's average cost is 68 µs per read. The last 7 % of reads cost
≈ 52 % of its total CPU time (718 of 1368 CPU-s on 20 M reads).

**What this means for a 10× target**
- **Against minibwa:** by Amdahl's law, if the hard 7 % keep minibwa's
  per-read cost, the end-to-end speedup cannot exceed ~1.9×, however fast
  the certified path is.
- **Coverage amortization** is the most promising different idea: learn
  shared donor events once and patch the reference, so sibling reads
  become certifiable. It addresses ~26 % of hard reads; better
  certification addresses another ~16 %. Even with both, the ceiling
  against minibwa is roughly 2–3×.
- **10× over minibwa is not supported by these data.**
- **10× over BWA-MEM2**, the clinical reference aligner, is plausible:
  minibwa itself is 3.8× faster than BWA-MEM2 here, and the pipeline
  above is 4.1×.

## 6. Index builds (built concurrently in 3 groups; times are upper bounds)

| Tool | galaxy wall | a2 wall | Peak RSS (GiB) |
|---|---|---|---|
| CERTA | 0:24 | 0:32 | 12.5 |
| strobealign | 0:35 | 0:45 | 12.8 |
| minimap2 | 1:00 | 1:16 | 12.0 |
| minibwa | 4:36 | 9:08 | 55.8 |
| BWA-MEM2 | 25:13 | 54:38 | 69.3 |
| bowtie2 | 27:50 | 32:32 | 6.1 |

## Limitations

- **Single-end only**, because the prototype is single-end. Paired-end
  results may differ, since pairing changes both certification and the
  competitors' costs.
- **One sample and one run of HG002.** Timing uses a 20 M-read random
  sample and 3 repeats. This is a single-sample pilot, not the
  multi-sample protocol in `RESEARCH_ROADMAP.md` Part 4.
- **Not measured:** input decompression time is excluded, because all
  tools read uncompressed FASTQ. Sorting and duplicate marking were not
  measured.
- **Accuracy is agreement with BWA-MEM2**, not truth.
- **Shared servers:** galaxy's GPUs carried idle, resident ollama processes.
- **Version drift across rows:** the v0.3 CERTA rows were measured after
  the competitor rows, on the same idle machines.
