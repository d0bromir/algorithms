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

galaxy is ARM64 (aarch64, 128 cores) with 2× A100 80GB, CUDA and CMake.
a2 is x86-64 (64 cores) with no GPU driver, no CMake and no zlib headers.
There, `build.sh` compiles directly with g++, and gzip input is disabled.

```bash
ssh galaxy   # or a2
git clone https://github.com/d0bromir/algorithms.git && cd algorithms/certa
scripts/build.sh                   # on galaxy: CUDA_ARCHS=80 scripts/build.sh
THREADS=64 DEVICE=1 scripts/demo.sh 50000000 2000000
```

On a GPU host, `demo.sh` runs the CPU and GPU back-ends on the same reads
and fails if their outputs differ. On galaxy, `DEVICE=1` selects the second
A100; check `nvidia-smi` first, because the GPUs are shared (ollama
servers were resident on both in September 2026).

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

For the most certified reads, use `--full-cert` on an index of every
position (`certa index ... -s 1`, about 23 GB for GRCh38): pass 3, pass V
(adaptive parts) and pass L. On 20 M HG002 reads it certifies 97.6 % of
reads single-end and 97.4 % paired (default: 94.4 % and 95.4 %) in about
2.5 times the default time (galaxy, 64 cores: 38.7 s and 100.1 s; still
2.6× and 3.3× faster than BWA-MEM2); a second index (`--index2`) adds only
0.1 points more. Peak memory is the index alone, whatever the input size:
21.3 GiB single-end and 21.5 GiB paired for GRCh38 (24.3 / 25.0 GiB with an
A100), as for the default (`bench/results/hg002_2026-10/memory_*.tsv`).

### Interval output, no fallback (`--intervals`)

```bash
build-*/certa map grch38.cidx R1.fastq.gz R2.fastq.gz --gpu --intervals -o all.sam --stats stats.json
```

Every read goes to SAM; none is passed to another aligner.
- **Uncertified reads** go through tier R (`include/certa/repeat.h`).
- **Pairs with an uncertified mate** are placed from the alignments found,
  plus an exhaustive search of each mate's proper-pair window
  (`include/certa/interval.h`).

Each record carries:

| tag | meaning |
|---|---|
| `XL`, `XU` | the best local score over the whole reference (BWA-MEM scoring, both strands) lies in [XL, XU]; XL = XU means the alignment is proven optimal |
| `XQ` | a certified lower bound of BWA-MEM's MAPQ formula, from a proven upper bound on the best score at any non-overlapping locus |
| `XT` | the tier: `S0`/`S1`/`SR`/`SL`/`PR` certified; `R` exact by tier R; `IV` an interval; `MR` placed in the mate's window; `UN` unmapped |

**MAPQ.** By default MAPQ is a point estimate: BWA-MEM's formula, using the
best other locus found. A pair is trusted no more than its more confident
mate. `--iv-mapq cons` writes XQ instead.

**Tier R options:** `--iv-tau 256`, `--iv-margin 20`, `--iv-chains 8`,
`--iv-fallback 8`.

Results: [EXACT_ALIGNMENT_PROPOSAL.md §6.3](../EXACT_ALIGNMENT_PROPOSAL.md).

### Checking the output independently (`certa check`)

```bash
build-*/certa check-index grch38.cidx GRCh38.fa         # once per reference build (34 s for GRCh38)
build-*/certa map grch38.cidx R1.fq.gz R2.fq.gz --gpu --intervals --witness -o all.sam
build-*/certa check grch38.cidx all.sam                  # every record
build-*/certa check grch38.cidx all.sam --sample 100000  # or an audit sample
```

- **`--witness`** adds:
  - `XW`, the probes and support the bounds were proven from;
  - `XH`, the SUB bound behind XQ;
  - `XV`, the alignment attaining XL when it is not the record's own (for
    a PR mate: its single-end alignment);
  - `XY` (PR mates), the proper-pair distance D and each mate's parts of
    the pair certificate, which prove `XS`, the bound on every other
    proper pair.
- **`certa check`** (`src/check.cpp`) shares no alignment, bound or lookup
  code with the mapper. For each record it verifies:
  - AS, by rescoring;
  - XL, by a shown alignment;
  - XU and XH, by recomputing the bound and enumerating every probe
    occurrence itself;
  - XQ, against BWA-MEM's formula. For PR mates it uses BWA-MEM's pairing
    rule, with XS proven from both records of the pair: every chain of
    each mate's part occurrences is paired, by exhaustive DP, with every
    partner placement within D.
- **It exits** with status 2 if any record fails; `--failures f.txt` lists
  them. The two records of a pair must be adjacent (as `certa map` writes
  them); `--sample` samples pairs whole.

Results: [EXACT_ALIGNMENT_PROPOSAL.md §6.4](../EXACT_ALIGNMENT_PROPOSAL.md).

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
| `--serial-stages` | off | without a GPU, run a batch in two stages, one after the other (passes 1 and 2 over the whole batch, then the CPU passes and formatting), instead of all stages in one pass over chunks of reads (same output either way) |
| `--budget3` / `--support3` | 0 / 2 | pass 3 (CPU): reads whose parts did not all fit pass 2's budget are re-run with up to this many hits per part (≤ 4096), verifying clusters hit by ≥ `--support3` parts |
| `--local` | 0 | pass L (CPU): certified local alignment (tier SL) for reads still uncertified, verifying chains hit by ≥ t parts (2 is the cheap setting) |
| `--full-cert` | off | maximum certification: `--budget3 4096`, `--adaptive-parts 2` (unless `--witness`) and `--local 2` unless given; slower |
| `--adaptive-parts` / `--vbudget` | 0 / 4096 | pass V (CPU, s = 1 index): for reads still uncertified after pass 3, parts chosen per read: 1 = the rarest disjoint q-mers (radius \|S\| − 1); 2 = then q-mers searched with up to one edit (radius 2\|S\| − 1, at most 12); hits per strand ≤ `--vbudget`. Not with `--witness` |
| `--index2` | none | not recommended since pass V: a second index of the same reference (e.g. built with `-q 15 -s 1`, about 23 GB for GRCh38) on which reads still uncertified after every pass get every pass again (CPU). With `--full-cert` it adds 0.1 points for twice the memory and 50 % more time (HG002). Not with `--witness` or `--serial-stages` |

## Results on real data: GIAB HG002, NovaSeq X

The current numbers are in the paper ([paper/main.tex](paper/main.tex),
Table 2 for speed), every one generated from `bench/results/` by
[paper/make_numbers.py](paper/make_numbers.py). Earlier runs and their
caveats are in [bench/RESULTS.md](bench/RESULTS.md).

- **Certified fraction:** 95.4 % of all reads in paired mode (HG002), 92.9 %
  for HG001 and HG005 (NovaSeq 6000); 94.4 % single-end.
- **Speed (20 M reads or pairs, every pipeline on the same 64 cores, one
  session per host, build 77b1df2):** on CPU alone CERTA is faster than
  every other aligner timed, single-end and paired, on galaxy (ARM64) and
  a2 (x86-64). Against minibwa, the fastest of them, on galaxy / a2:
  single-end 1.88× / 2.02× in interval mode and 1.59× / 1.85× with the
  minibwa fallback; paired 1.66× / 2.09× in interval mode and 1.21× /
  1.38× with the fallback. With one A100 (galaxy): 2.32× single-end, 2.22×
  paired. BWA-MEM2 is 6.0–8.8× slower than CERTA's interval mode on CPU.
  The current build (b4c7696, same output; galaxy, one session, every
  aligner re-timed): interval mode 2.01× / 1.80× minibwa single-end /
  paired on CPU, 2.36× / 2.31× with the A100, in 21–25 GiB instead of
  29–40 GiB ([final_galaxy_b4c7696.tsv](bench/results/hg002_2026-10/final_galaxy_b4c7696.tsv)).
- **Independent checker:** `certa check` accepts all 8 M honest records
  (HG002 and HG005) and rejects all 24,306 injected faults.
- **Maximum certification against the slow, reliable mappers** (2 M HG002
  reads or pairs, galaxy, 64 cores, CPU, interval mode;
  [bench/results/hg002_2026-10/slow_galaxy.tsv](bench/results/hg002_2026-10/slow_galaxy.tsv)):

  | | certified, single-end / paired | single-end | paired |
  |---|---|---|---|
  | CERTA (default) | 94.5 % / 95.4 % | 2.7 s | 5.1 s |
  | CERTA `--full-cert --adaptive-parts 0` (no pass V) | 95.7 % / 96.3 % | 3.0 s | 6.9 s |
  | CERTA `--full-cert` (with pass V) | 97.6 % / 97.5 % | 5.6 s | 12.5 s |
  | CERTA `--full-cert --index2` (q = 15) | 97.7 % / 97.5 % | 8.5 s | 17.0 s |
  | Yara 1.0.5 `-e 4 -y full` | — | 138.2 s | 211.2 s |
  | RazerS 3 `-i 96 -rr 100` | — | est. ~2,400 s (timed: 239.6 s for the first 200 k reads, CERTA 1.6 s) | aborts (`std::bad_alloc`) |

  Yara and RazerS 3 find every alignment within 4 % errors (6 edits at
  150 bp); they give no proof. CERTA's certificates (radius ≤ 5, with pass V
  up to 12) agree with them: within their own error rate they miss no read
  CERTA certifies and place none elsewhere with more edits than the proven
  minimum; the reads they leave unmapped have their proven optimum beyond
  that rate (0.89 % of single-end reads), and CERTA certifies those too
  ([slow_reliability.txt](bench/results/hg002_2026-10/slow_reliability.txt)).
  Long reads: CERTA's HiFi pipeline (minimap2 + `tools/lr_ceiling --certify`,
  optimum proven for 93.3 % of 242,534 reads) takes 163 s, Winnowmap2 610 s.

## Results so far (synthetic data)

### Lab hosts (September 2026)

Synthetic 50 Mbp genome, 2 M × 150 bp reads, k = 2, same inputs on every host.

| Host / back-end | Map step | End to end | Certified |
|---|---|---|---|
| a2: x86-64 CPU, 64 of 64 cores (GCC 13, built without CMake) | 2.28 M reads/s | 4.3 s | 94.83 % |
| galaxy: ARM64 CPU, 64 of 128 cores (GCC 15) | 3.38 M reads/s | 3.5 s | 94.83 % |
| galaxy: NVIDIA A100 80GB PCIe (CUDA 13.1, sm_80) | 8.85 M reads/s | 3.2 s | 94.83 % |

- **Identical output on every back-end.** The SAM and uncertified-FASTQ
  checksums match across a2 (x86-64), galaxy's ARM64 CPU and the A100.
- **The certificate test passes on both hosts.**
- **Parsing now dominates.** End-to-end time is dominated by
  single-threaded FASTQ parsing and encoding (~2.3–2.5 s of the totals
  above), not by mapping. A parallel parser is therefore the next
  speed-up (Amdahl's law).

### Development machine

| Setting | Result |
|---|---|
| Certificate test (oracle), 10 seeds × 480 reads × 5 (k, cap) configs, x86-64 | all checks pass |
| Mutation tests (drop a sampling shift; narrow the band) | caught by the test |
| Synthetic 20 Mbp genome (repeats, SDs, STRs, N gaps), 500 k × 150 bp reads, 0.3 % sub + 0.01 % indel, k=2 | 95.7 % certified (60.8 % S0 exact, 34.9 % S1); 2.5 % repetitive (radius < 0), 1.7 % > R edits |
| Placement vs. simulated origin, certified reads | MAPQ 60: 100 %, MAPQ 40: 100 %, MAPQ 20: 99.94 %, MAPQ 0 (ties): 43 % |
| CPU, 8 threads (laptop, WSL2) | 357 k reads/s (map step) |
| GPU, RTX 3050 Laptop (sm_86) | 1.57 M reads/s (map step, including transfers); output byte-identical to CPU |
| Kernel resources (sm_80) | 48 registers, 27 KB stack/thread (≈6 GB local-memory reservation on an A100) |

These are synthetic-data numbers. They show that the mechanics work, that
the certificate holds, and that all back-ends give identical output. They
are **not** evidence about real genomes, where repeat content, error
profiles and the variant spectrum differ. Not yet run: GRCh38, real
NovaSeq X FASTQ (the Aim 1 measurement).

## Layout

```
include/certa/core.h       per-read algorithm (host + device)
include/certa/pair.h       pair certificate (tier PR)
include/certa/repeat.h     repeats (tier R)
include/certa/interval.h   interval output (--intervals)
include/certa/longread.h   long reads: Lemma L bounds, incl. split (jump) alignments
src/index.*                FASTA loading, sampled q-mer index, save/load
src/seqio.*                FASTQ/FASTA reading (plain or gzip)
src/parallel.h             a loop spread over threads (FASTQ reading, encoding)
src/check.*                independent checker (certa check, check-index)
src/mapper_cpu.cpp         multi-threaded CPU driver
src/mapper_gpu.cu          CUDA driver (one thread per read)
src/main.cpp               CLI: index, repeat, check, check-index; dispatch
src/cmd_map.cpp            certa map: load, open outputs, run the pipeline, report
src/map_config.*           certa map's options, parsed and checked up front
src/map_pipeline.*         batches through passes 1-3, V, L, pairs, tier R and output
src/sam_output.*           SAM and fallback FASTQ records, MAPQ
src/map_stats.*            counts, timing, the stderr summary and --stats JSON
src/cli.h, async_writer.h  argument parsing; an output file written by its own thread
tests/test_certa.cpp       brute-force oracle tests of every bound
tests/verify_lemmas.cpp    exhaustive checks of the finite lemmas
tests/test_longread.cpp    long-read certificates against a brute-force DP (planted SVs)
tests/mutants.tsv          one mutant per bound; tests/run_mutants.sh runs them
tools/simulate.py          synthetic genome/reads and placement evaluation
tools/lr_ceiling.cpp       long reads: Lemma L ceiling vs best known alignment; --certify proves OPT (Aim 6)
tools/payoff.cpp           part choices taken from the slow mappers, certified on CERTA's uncertified reads
scripts/                   build.sh, demo.sh
THEORY.md                  the proofs, with code, tests and mutants per statement
lean/                      machine-checked proofs (Lean 4; lake build): Lemmas L, E, H, S, W, C; Theorem G; band pruning, DP row
paper/                     the paper (LaTeX); make_numbers.py and make_figures.py regenerate every number and figure from bench/results
```

**Mutation tests** (from `certa/`; about an hour on 12 cores):

```bash
tests/run_mutants.sh 12
```

Each mutant must fail an oracle check, or, when THEORY.md proves it changes
no claim, pass.

## Known limitations and next steps

- **Single-end only.** Add mate-pair consistency and rescue.
- **MAPQ:** a provisional 20 × (distance gap) rule, capped at 60. Calibrate
  it on simulated reads from HPRC assemblies (roadmap §2.4).
- **Memory:** the index is memory-mapped and is nearly all of it: positions
  11.0 GiB, keys 5.5 GiB, the reference 2.9 GiB (1 byte per base) and the
  directory 1 GiB for GRCh38 with `-s 1`. A 2-bit reference (2.2 GiB less)
  and 8-bit keys (2.7 GiB less) were measured and cost too much time
  (9–20 % and 21–79 %; FINDINGS.md).
- **GPU:** one thread per read, with local-memory workspaces. Warp-cooperative
  lookups, overlapping transfers with compute, and a device-side output
  path are the obvious optimizations.
- **Clusters wider than 48 diagonals** (long tandem repeats) are left
  uncertified instead of being split.
- **FASTQ input is parsed on one thread**, which is now the end-to-end
  bottleneck. Parse (and gunzip) in parallel and overlap it with mapping.
- **Fallback merge:** combine the certified SAM with the fallback
  aligner's output (e.g., `samtools merge` + sort) in one fused pipeline
  (roadmap §2.6).
