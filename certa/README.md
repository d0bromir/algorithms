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

**Smaller indexes.** An entry's 16 key bits below the bucket are also in
the reference at its position.
- **Conversion:** `certa index --from grch38.cidx -o grch38.k8.cidx --key-bits 8`
  (or `0`) rewrites an index with 8 of them or none, in seconds and without
  the FASTA. `--key-bits` also applies when building.
- **Memory:** the GRCh38 `-s 1` index takes 17.6 or 14.9 GiB instead of
  20.4 GiB; with `--t2t`, both indexes take 36.0 or 30.4 GiB instead of
  41.7 GiB.
- **Cost:** the output is byte-identical, at 21 % or 35 % more CPU time
  (`bench/VALIDATION_RESULTS.md`, "Smaller indexes").

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
mate. `--iv-mapq cons` writes XQ instead; `--iv-mapq cert` writes XQ only
for tiers R, IV and MR, the reads the default certificate left.
- The point estimate is overconfident for those tiers: on real HG002 reads,
  21 % of the IV mates and 11 % of the MR mates at MAPQ ≥ 40 are
  misplaced, and their XQ is 0.
- With `--iv-mapq cert`, interval mode is as accurate at MAPQ ≥ 40 as the
  default mode plus minibwa with its MAPQ corrected, with no second aligner
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#mapq-of-the-uncertified-mates)).
- XQ for every record (`cons`) leaves few reads confident and drops the
  reads that carry variants ([FINDINGS.md](FINDINGS.md), Aim 3).
- `cert` costs variant calls too. With DeepVariant on 30× HG002 (chr20 and
  the CMRG genes), the same placements give 3,156 SNV errors with `cert`
  against 1,956 with `point` (minibwa: 1,745)
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#variant-calling)).
  So `point` stays the default.

**Pairs and the fitted MAPQ** (paired interval mode, opt-in;
[VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#re-paired-mates-and-a-fitted-mapq)):
- `--iv-pair all`: two certified mates that do not form a proper pair are
  re-paired as an uncertified mate is: a mate's best alignment in the
  other's proper-pair window is taken when the pair beats the two
  single-end bests minus the unpaired penalty, and becomes MR.
- `--iv-rescue local`: mate rescue by local alignment, as BWA-MEM's. The
  window search above finds a mate end to end, above about L − 54 and with
  gaps up to its band. When it finds no proper pair, a mate without a
  certificate is aligned in the proper-pair window of a confident partner
  (MAPQ ≥ 40 by `--iv-mapq cert`'s rule) by a full DP of the window, so with
  clips and gaps of any length, and becomes MR. This places chimeric reads
  (51M100S), divergent ones (11–21 mismatches) and long deletions that
  CERTA otherwise leaves unmapped or elsewhere. Certified mates are never
  moved, and a pair found end to end is never replaced. On HG002 it removes
  69 of DeepVariant's SNV errors and 157 indel errors (a third and a half
  of the gap to minibwa), for 5 % more CPU, with no rise in confident
  errors at matched N
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#a-local-mate-rescue)).
- `--mapq cal`: MAPQ is the Phred of the misplaced share of similar records
  on HG002 (`include/certa/mapq_table.h`, fitted by
  `bench/truth_mapq_fit.py`): by tier, edits, the gap to the second-best
  locus within 10 edits and the clusters verified. It ranks reads near a
  paralog lower, where the formula gives 60. The values are error rates,
  so fewer reads reach MAPQ 40 (74 % of real HG002 mates).
- `--mapq pair`: MAPQ from boosted trees over both mates' records
  (`include/certa/pair_model.h`, trees in `include/certa/mapq_model.h`,
  fitted at the locus on HG002 by `bench/truth_mapq_pair.py`). Besides a
  record's own tier, edits, gap and clusters, the trees read the other
  mate's tier, MAPQ and gap, a PR mate's pair margin and the insert size:
  whether a mate in a repeat has a unique partner decides most of the
  ranking among the less confident mates
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#a-mapq-from-both-mates)).
- `cal` and `pair` rank misplaced reads better, but they cost variant calls.
  Each edit lowers MAPQ, the sample's own variants included, so with
  DeepVariant they give 1.6–3.1× the SNV errors of `point`
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#variant-calling)).
  Use `point` (the default) for variant calling.
- `--t2t T2T.cidx --t2t-chain T2TToRef.over.chain.gz`: every pair is also
  mapped, the same way, on an index of T2T-CHM13 in the same pass, and each
  mate's T2T placement is lifted to the reference through the chain (UCSC
  `hs1ToHg38.over.chain.gz`). With `--mapq pair` the trees also read it
  (`include/certa/mapq_model_t2t.h`): whether it lifts to the record's
  locus, elsewhere or not at all, and its gap and MAPQ there. A mate from a
  copy GRCh38 lacks tends to align elsewhere on T2T-CHM13. The records,
  their placements and certificates are the reference's; only MAPQ
  changes. It needs both indexes in memory (46 GiB for GRCh38 and hs1) and
  about doubles the mapping time.
- `--feature-tags` adds `XO` (the best other candidate's distance) and `XN`
  (clusters verified), the fits' inputs not in other tags; with `--t2t`
  also `XC:Z:class,tier,MAPQ,XE,XD,XO`, the mate's T2T placement as the
  trees read it.

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

### Maximum-likelihood placement (`--bio`)

```bash
build-*/certa map grch38.cidx R1.fastq.gz R2.fastq.gz --bio --bio-insert 400,100 \
    -o ml.sam -u fallback_R1.fq -U fallback_R2.fq --stats stats.json
```

The default certificates prove the optimum of an optimization problem
(fewest edits, BWA-MEM's score), which need not be where the read came from.
`--bio` proves the optimum of the read's **likelihood** instead, under an
explicit model of recent Illumina data and of human variation
([THEORY.md §12](THEORY.md#12-maximum-likelihood-placement---bio)):

- each base call is wrong with the probability its quality states, so two
  mismatches at Q2 bases outweigh one at a Q40 base;
- the sample differs from the reference by SNVs (transitions 2.1 times as
  frequent as transversions) and indels, and indels open far more easily in
  homopolymers and short tandem repeats (gap opening 46 Phred, down to 22);
- pairs are placed jointly, with a fragment-length term, so a mate in a
  repeat follows its unique mate.

For each read it reports the **maximum-likelihood locus over the whole
reference, both strands** (proven, as the other tiers are: every alignment
cheaper than the horizon `YH` is aligned exactly), the second locus, the
ties, and a **posterior MAPQ**. Indels are left-aligned. Reads whose optimum
cannot be proven (it is not below `YH`) go to the fallback, as in the
default mode, unless their mate is certified: then the mate's proper-pair
windows are aligned exactly, and the read is placed there when that pair is
proven optimal (**pair rescue**, `YR:i:1`; THEORY.md Proposition B3). Reads
with adapters or poly-G tails go to the fallback too, since alignment is
end to end (trim them first). One CPU pass replaces passes 1–3, V and L, so
`--bio` does not combine with `--gpu`, `--intervals`, `--index2` or
`--serial-stages`.

| tag | meaning |
|---|---|
| `XT:Z:BM` | certified maximum-likelihood placement |
| `YP` | penalty of the reported alignment: −10 log₁₀ of its likelihood relative to a perfect match, in Phred; the certified least over the reference unless the mate moved the read (`YT`) |
| `YS` | the best other locus's penalty (exact), or −1: no other locus below `YH` |
| `YH` | horizon: every alignment at no listed locus costs ≥ `YH` |
| `YB` | loci tied at the least penalty (MAPQ ≤ 3 when ≥ 2) |
| `YA` | how far the alignment can shift at the same penalty within its locus (a tandem repeat's period multiples), 0 if it cannot; the record is the leftmost |
| `YT`, `YI` | placed jointly with the mate: the pair's penalty (both mates and the insert term) and its insert term |
| `YR:i:1` | placed only by the mate (pair rescue): not certified alone, so `YP` ≥ `YH`, and `YH` bounds every alignment of the read |
| MAPQ | the locus's posterior over the loci below `YH` plus one competitor at `YH`, Phred, ≤ 60 |
| `NM`, `AS`, `XE`, `XR` | as for S0/S1 (`XE`: the certified fewest edits over the reference, when the edit certificate holds) |

| option | default | meaning |
|---|---|---|
| `--bio` | off | certified maximum-likelihood placement |
| `--bio-snv` / `--bio-titv` | 0.001 / 2.1 | SNVs per base; transitions per transversion |
| `--bio-indel` / `--bio-seq-indel` / `--bio-gap-ext` | 0.0001 / 0.00001 / 0.5 | indels per base (sample, sequencing); probability an indel continues one more base |
| `--bio-str-min` / `--bio-str-slope` / `--bio-str-cap` | 4 / 3 / 24 | a gap touching a tandem run (period 1–4) of ≥ min bases opens slope Phred cheaper per base, at most cap |
| `--bio-insert` / `--bio-improper` | 400,100 / 40 | fragment length mean,sd; penalty of a pair that is not proper |
| `--bio-rescue-loci` | 8 | a certified mate's loci whose proper-pair windows the pair rescue aligns (0: no rescue) |
| `--bio-run-table` | 1 | tabulate the reference's tandem runs once (half a byte a base: 1.45 GiB for GRCh38) instead of scanning each band's window; same output, 2.65× faster on HG002 Test A (the SIMD lanes are filled from it 16 bases at a time); 0 saves the memory |
| `--bio-budget` | 4096 | hits enumerated per strand (more parts: a higher horizon) |
| `--bio-gap-cap` / `--bio-delta` / `--bio-max-band` | 50 / 60 / 1024 | band reach; clusters bounded δ above the best are not aligned; widest band aligned |

**Simulated reads with known origin** (20 Mbp synthetic genome with
diverged repeat families, segmental duplications and tandem repeats;
`tools/simulate.py reads --quality novaseqx`, the bins and 3′ decline of
GIAB HG002 NovaSeq X; sample SNVs/indels shared by both mates; 8 threads):

| | in SAM | at the truth (±10 bp) | MAPQ ≥ 20: records / misplaced |
|---|---|---|---|
| 300,000 single reads, default | 97.79 % | 94.33 % | 258,911 / 21 |
| 300,000 single reads, `--bio` | 97.83 % | 94.77 % | 263,369 / 27 |
| 150,000 pairs, default | 99.37 % | 95.40 % | 268,729 / 564 |
| 150,000 pairs, `--bio` | 97.90 % | 96.30 % | 261,665 / 23 |

- **Pairs:** 24× fewer confidently misplaced reads. The default mode's are
  mostly mates the pair certificate shifts within tandem repeats (542 of
  564), at MAPQ 6–24 values that are right only 34–42 % of the time; `--bio`
  places 95 % of reads jointly.
- **Wrong loci at MAPQ ≥ 40:** default 17 single reads and 16 mates
  (segmental-duplication copies certified at radius 0–1 and given MAPQ 60);
  `--bio` 1 and 0. Every other `--bio` record misplaced at MAPQ ≥ 40 (14
  single, 12 mates) is at its true locus, shifted by a tandem repeat's period
  where the sequence is identical, and carries `YA` > 0, which says so.
- **MAPQ:** `--bio` is conservative below 40 (single reads, MAPQ 20–29: 7
  misplaced where 34 are expected).
- **Where the locus changes:** on these Q40-dominated reads the likelihood
  optimum and the fewest-edits optimum are the same locus for every uniquely
  placed single read; they differ where qualities are low or indels compete
  with mismatches (the oracle test's planted families), and for 3.7 % of
  mates, which the mate moves.
- **Cost** (8 threads, WSL2 laptop): 300,000 single reads in 22.3 s against
  the default mode's 3.4 s (NovaSeq 6000 profile: 16.5 s against 5.7 s);
  150,000 pairs in 19.5 s against 8.8 s. The band DP drops every cell that
  already costs the horizon, which makes it 3–6 times faster with the same
  output.

**On real data** ([bench/VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md)):
HG002 reads with a known origin. The origins come from simulation from the
T2T-HG002 Q100 assembly (10 M pairs), or are read off that assembly for
19.9 M real NovaSeq X pairs; both are projected to GRCh38 through GIAB's
alignments.
- **Soundness:** none of the ~283 M alignments other aligners reported scores
  below `--bio`'s certified optimum.
- **Confident placements are not more accurate than minibwa's.** Scored at
  the locus (within 150 bp) and at the same number of confident mates (the
  95.6 % of real mates `--bio` calls confident), `--bio` misplaces 138 per
  million, minibwa 46 and CERTA's default mode 77; on simulated reads,
  1,309 against 262 and 785
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#matched-sensitivity-and-the-locus-window)).
  An earlier lead in avoidable errors came from a 20 bp window, which
  counts a read soft-clipped at an insertion as misplaced, and from calling
  fewer mates confident. Nearly all of every aligner's remaining confident
  errors are on mates whose origin GRCh38 lacks (0.05 % of real mates).
  Setting MAPQ 0 where T2T-CHM13 confidently places the read elsewhere
  removes about half of them on real reads, for every aligner alike
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#a-mapq-from-t2t-chm13)).
- **Interval mode with `--iv-pair all --mapq pair`** (below, fitted on
  HG002) is the most accurate aligner measured at every operating point,
  held out on HG005. Of the 99.75, 99.29, 98 and 95 % most confident mates,
  it misplaces 186, 60, 20 and 9 per million, against minibwa's 208, 80,
  37 and 37 and BWA-MEM2's 202, 108, 108 and 108. With `--t2t` (T2T-CHM13
  in the same pass) it is 93, 25, 7 and 3: 2–12× fewer than minibwa, and
  ahead of minibwa with a T2T second opinion from BWA-MEM2 (167, 48, 22,
  22). `--mapq cal` leads only at 98 % and below (22–33)
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#a-mapq-from-both-mates)).
- **But not for variant calling yet.** The fitted MAPQs give the reads that
  show a sample's variant a lower MAPQ than their neighbours, so DeepVariant
  misses variants:
  - with them, CERTA's records on 30× HG002 give 3,214–6,052 SNV errors,
    against minibwa's 1,745;
  - with the default MAPQ, the same placements give 1,956
    ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#variant-calling)).
- **MAPQ of the fallback:** minibwa misplaces 6–7 % of the fallback mates
  it gives MAPQ ≥ 40. Divided by 4 (`drop_certified_mates.awk -v div=4`),
  their MAPQ comes within 4× of their error rates
  ([VALIDATION_RESULTS.md](bench/VALIDATION_RESULTS.md#mapq-of-the-uncertified-mates)).
- **The cost:**
  - `--bio` leaves more reads to the fallback: with the pair rescue, 0.13 %
    of real mates with a known origin are not in the SAM file, against
    0.05 % for the default mode (98.9 % of the rescued mates are at their
    origin);
  - its MAPQ is overconfident between 30 and 49;
  - it runs 16–17× slower than the default mode on GRCh38 (up to 16 bands
    aligned at once in SIMD lanes), and needs about 23.8 GiB with the
    tandem-run table (22.3 GiB without);
  - the model has no complex indels (an insertion next to a deletion).

The resources and the plan are in [bench/TRUTH_DATA.md](bench/TRUTH_DATA.md)
and [bench/VALIDATION_PLAN.md](bench/VALIDATION_PLAN.md).

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

## QC and mapping in one pass: `certa-qc`

`certa-qc` runs fastp's QC, trimming and filtering and `certa map` in one
process. fastp reads the input, and the reads it keeps go straight into
certa map's pipeline. Nothing is compressed, written or read back in between.
This is milestone M1 (CPU) of
[FUSED_QC_MAPPING_DESIGN.md](../FUSED_QC_MAPPING_DESIGN.md).

It needs a [d0bromir/fastp](https://github.com/d0bromir/fastp) checkout with
the library built (`make lib`, see `docs/LIBRARY.md` there):

```bash
make -C ~/fastp -j && make -C ~/fastp lib
cmake -S certa -B build -DFASTP_DIR=$HOME/fastp
cmake --build build -j --target certa-qc
```

For a CUDA build of fastp, pass its libraries as a CMake list:
`-DFASTP_LIBS="$(make -s -C ~/fastp print-libs | tr ' ' ';')"`.

```bash
certa-qc ref.cidx [certa map options] --fastp [fastp options]

certa-qc GRCh38.cidx -t 32 -o certified.sam -u fallback_R1.fq -U fallback_R2.fq --stats map.json \
    --fastp -i R1.fq.gz -I R2.fq.gz -w 16 --detect_adapter_for_pe -j qc.json -h qc.html
```

- **Options.** Before `--fastp` come certa map's options, without the FASTQ
  inputs. After it come fastp's, with fastp's `-i`/`-I` as the input.
- **Refused fastp options.** `-o`, `-O` and `--stdout` are refused, because
  the reads fastp keeps are mapped. `--merge` and `--split` are not supported
  yet.
- **Dropped reads.** `--failed_out`, `--unpaired1`/`--unpaired2` and
  `--overlapped_out` are written as fastp writes them. A name ending in `.gz`
  is gzip-compressed.
- **Reports.** fastp's JSON and HTML reports (`-j`, `-h`) and certa's
  `--stats` are written as the two programs write them.
- **Threads.** fastp's work runs on its own `-w` workers, ahead of certa's
  `-t` threads, so trimming overlaps loading the index and mapping.

**The output equals the two-step pipeline** `fastp -o trimmed.fq ...` then
`certa map ref.cidx trimmed.fq ...`. `tests/qc_equivalence.py` runs both on
simulated data with adapter read-through, poly-G tails and poor 3' ends. It
compares the SAM, the uncertified FASTQ files and fastp's files for dropped
reads byte for byte, and both JSON reports (without timings).
Result on 40 000 simulated pairs: all 72 cases identical. The cases cover
single-end and paired input, six fastp option sets, four certa modes,
`-w`/`-t` 1 and 4, and plain and gzip input. fastp lowers `-w` for small
inputs. At `-w 4` it used four workers on the plain input and two on the
gzip input.

The option set with `--dedup` runs at `-w 1` only. With more than one
worker, fastp's own output changes from run to run. The copy of a duplicate
it keeps depends on the order of its packs. In a first run at `-w 4`, these
were the only 4 of 80 cases that differed.

**Speed on four cores: about the same as a plain intermediate file.**
`bench/qc_speed.py` times both on larger test data: 400 000 pairs, gzip
input, on a 4-core virtual machine. It ran three times, three rounds each.
Every run wrote the same SAM. The logs and the builds are in
[bench/results/qc_speed_2026-11/](bench/results/qc_speed_2026-11/README.md).
Medians per run, in seconds:

| Pipeline | fastp `-w` | certa `-t` | Wall time | CPU time |
|---|---|---|---|---|
| fastp, gzip file, certa map | 4 | 4 | 10.64, 10.63, 10.99 | 35.9–37.4 |
| fastp, plain file, certa map | 4 | 4 | 9.50, 9.62, 9.75 | 32.0–33.2 |
| `certa-qc` | 4 | 4 | 9.59, 9.80, 10.16 | 33.3–35.2 |
| fastp, plain file, certa map | 1 | 1 | 31.12, 31.07, 30.72 | 32.6–33.1 |
| `certa-qc` | 1 | 1 | 27.46, 28.53, 27.81 | 34.3–35.0 |

- **Four threads each.** `certa-qc` takes 1–4 % longer than the pipeline
  with a plain file, and 8–10 % less than with a gzip file.
- **Why no gain.** Mapping takes 80 % of the CPU time. In run 3, with a
  plain file and four threads each, certa map used 25.4–27.0 s and fastp
  6.2–6.8 s. All four cores are already mapping, so fastp's work takes
  cores away from the mapping.
- **One thread each.** Cores are free, so fastp's work runs beside the
  mapping. `certa-qc` takes 8–12 % less time.
- **CPU time.** `certa-qc` uses 0–9 % more CPU time than the pipeline with a
  plain file. Where that time goes is not measured yet.
- **Not tested.** GRCh38, more cores, and intermediate files on a disk
  rather than in the page cache.

**Limits (next steps in the design):**

- **Hard trim.** Reads are written as fastp trimmed them. Soft clips for the
  trimmed bases, the design's default, are not implemented yet. The library
  already reports each read's front offset and its input read for this.
- **Two thread pools.** fastp and certa map have separate threads, so `-w`
  and `-t` together can ask for more threads than there are cores. The
  design shares one pool. That removes the excess threads, not fastp's CPU
  work.
- **Copies.** Each kept read is copied from fastp's record into certa's.
  Handing device batches over without a copy is milestone M2.
- **GPU.** `--gpu` uses certa map's GPU mapper as `certa map` does. It has
  not been tested with `certa-qc` yet.

## Results on real data: GIAB HG002, NovaSeq X

The current numbers are in the paper ([paper/main.tex](paper/main.tex),
Table 2 for speed), every one generated from `bench/results/` by
[paper/make_numbers.py](paper/make_numbers.py). Earlier runs and their
caveats are in [bench/RESULTS.md](bench/RESULTS.md).

- **Certified fraction:** 95.4 % of all reads in paired mode (HG002), 92.9 %
  for HG001 and HG005 (NovaSeq 6000); 94.4 % single-end.
- **Speed (20 M reads or pairs, every pipeline on the same 64 cores, one
  session per host, build b4c7696 / 00789fe):** on CPU alone CERTA is
  faster than every other aligner timed, single-end and paired, on galaxy
  (ARM64) and a2 (x86-64). Against minibwa, the fastest of them, on galaxy
  / a2: single-end 2.01× / 2.02× in interval mode and 1.69× / 1.83× with
  the minibwa fallback; paired 1.80× / 2.14× in interval mode and 1.28× /
  1.40× with the fallback. With one A100 (galaxy): 2.36× single-end, 2.31×
  paired. BWA-MEM2 is 6.2–9.0× slower than CERTA's interval mode on CPU.
  Peak memory 21–25 GiB. The 77b1df2 session (1.88× / 2.02× and 1.66× /
  2.09× in interval mode, 29–40 GiB) is kept as `final_*_77b1df2.tsv`.
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
include/certa/bio.h        --bio: likelihood model, Lemma B, certified ML placement (Theorem B), pairs, pair rescue
src/index.*                FASTA loading, sampled q-mer index, save/load
src/seqio.*                FASTQ/FASTA reading (plain or gzip)
src/parallel.h             a loop spread over threads (FASTQ reading, encoding)
src/check.*                independent checker (certa check, check-index)
src/mapper_cpu.cpp         multi-threaded CPU driver
src/mapper_gpu.cu          CUDA driver (one thread per read)
src/main.cpp               CLI: index, repeat, check, check-index; dispatch
src/cmd_map.cpp            certa map: load, open outputs, run the pipeline, report (run_map)
src/qc_main.cpp            certa-qc: fastp (libfastp_core) as the input of certa map's pipeline
src/map_config.*           certa map's options, parsed and checked up front
src/map_pipeline.*         batches through passes 1-3, V, L, pairs, tier R and output
src/sam_output.*           SAM and fallback FASTQ records, MAPQ
src/map_stats.*            counts, timing, the stderr summary and --stats JSON
src/cli.h, async_writer.h  argument parsing; an output file written by its own thread
tests/test_certa.cpp       brute-force oracle tests of every bound
tests/verify_lemmas.cpp    exhaustive checks of the finite lemmas
tests/test_longread.cpp    long-read certificates against a brute-force DP (planted SVs)
tests/test_bio.cpp         --bio against a whole-reference likelihood oracle; Lemma B against an exact DP
tests/qc_equivalence.py    certa-qc against fastp, then certa map (byte-identical outputs)
tests/mutants.tsv          one mutant per bound; tests/run_mutants.sh runs them
tools/simulate.py          synthetic genome/reads and placement evaluation
tools/lr_ceiling.cpp       long reads: Lemma L ceiling vs best known alignment; --certify proves OPT (Aim 6)
tools/payoff.cpp           part choices taken from the slow mappers, certified on CERTA's uncertified reads
tools/bio_score.cpp        the --bio model's penalty of any aligner's SAM alignments (soundness at scale)
tools/bio_probe.cpp        per read: why --bio certifies or falls back, and the work (parts, horizon, bands, time)
bench/truth_*.py, *.sh     truth on GRCh38 for HG002 reads (Q100 haplotypes, GIAB alignments), scoring, reports
bench/qc_speed.py          certa-qc against fastp, then certa map: wall and CPU time
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
