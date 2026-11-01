# `--bio` on real data: phase 1 results (2026-10-31)

These are the pilot, Test A and Test B of [VALIDATION_PLAN.md](VALIDATION_PLAN.md),
run on galaxy (ARM64, 64 of 128 cores for each aligner).
- **Build:** `main` at 8b6c7de, with the `--bio` defaults unchanged (frozen
  before the first real run). The only option set is `--bio-insert 429,99`,
  the insert size measured in the pilot.
- **Full reports:** [results/truth_2026-10/](results/truth_2026-10/).
- **Scripts:** [truth_download.sh](truth_download.sh),
  [truth_prep.sh](truth_prep.sh), [truth_run.sh](truth_run.sh).
- **Code:** [truth_lift.py](truth_lift.py), [truth_label.py](truth_label.py),
  [truth_eval.py](truth_eval.py), [../tools/bio_score.cpp](../tools/bio_score.cpp).

## In short

- **The certificate holds at scale.** Five other aligners reported about
  283 M end-to-end alignments on reads `--bio` certified. None costs less
  under the model than `--bio`'s certified optimum. `bio_score` reproduces
  `YP` on all of `--bio`'s own records.
- **The confident placements are the most accurate measured.** On real
  NovaSeq X reads with a known origin, `--bio` has the fewest *avoidable*
  errors: misplaced mates that another aligner placed at the truth.
  - At MAPQ ≥ 40: 103 errors, against 439–988 for the others; per million
    confidently placed mates, 4–9× fewer.
  - At MAPQ ≥ 60: 31, against 325–725; 9–21× fewer per million.
  - On simulated reads: 2–4× fewer per million at MAPQ ≥ 40, 1.4–3.7× at MAPQ ≥ 60.
- **But not a gold standard yet:**
  - It leaves more reads to the fallback: 0.3 % of real mates with a
    truth (default mode 0.05 %), 1.6 % of simulated ones (1.1 %).
  - It gives MAPQ ≥ 40 to about 4 % fewer real mates (14 % fewer simulated).
  - Its MAPQ is overconfident in the 30–49 bands.
  - It is 89–96× slower than the default mode on the human genome.
  - The model excludes complex indels (an insertion next to a deletion).
- **Most confident errors cannot be fixed on GRCh38.** No aligner placed
  them at the truth. Of `--bio`'s misplaced mates with MAPQ ≥ 40, that holds
  for 98 % of real reads and 98 % of simulated ones. Their best GRCh38 match
  is not their orthologous locus.

**Since phase 1** ([below](#after-phase-1-the-fallback-and-the-cost)): a
pair rescue places mates the single-read certificate cannot. It lowers the
fallback to 3.48 % (A) and 4.11 % (B) of reads with the same avoidable errors.
A table of the reference's tandem runs and leaner data structures make
`--bio` 1.66–1.78× faster (50–58× the default mode instead of 89–96×), for
1.45 GiB more memory.

## Pilot (2 M HG002 NovaSeq X pairs, trimmed)

| | default | `--bio` |
|---|---|---|
| wall time (64 cores) | 3.3 s | 211 s (64×) |
| candidate loci per read | 12.2 | 479 |
| reads certified | 94.03 % (both mates) | 95.18 % |
| insert size (confident proper pairs) | median 429, IQR 367–501 (sd ≈ 99) | — |

## The truth sets

**Test A: simulated.**
- 10 M pairs, 2×150, with the NovaSeq X quality profile and no added
  variants (the sample's variants are in the haplotypes).
- Source: HG002 Q100 v1.1, as one male genome: maternal autosomes plus X
  and Y, paternal autosomes.
- Projected to GRCh38 through GIAB's dipcall alignments (1,593 records, 1.55 M
  blocks; 40.6 Mb of haplotype covered twice gets no truth).
- 94.09 % of mates have a GRCh38 truth. The other 5.91 % come from sequence
  with no 1:1 counterpart: satellites, acrocentric arms, SVs.

The projection check ([projection_check_A.txt](results/truth_2026-10/projection_check_A.txt))
used the 43,105 error-free mates among the first 100 k pairs:
- 85.0 % equal GRCh38 at their projected truth, 9.5 % differ by one edit,
  2.6 % by two (the sample's own variants).
- 99.71 % align within 10 edits within ±100 bp of the truth.
- The remaining 0.29 % have a GRCh38 orthologue that diverged from them.

**Test B: real.** 19,943,583 HG002 NovaSeq X pairs (SRR37356339, trimmed
with fastp: adapters, poly-G, length ≥ 30). The steps
([labels_B.txt](results/truth_2026-10/labels_B.txt)):
1. Map the pairs with CERTA's default mode to each haplotype genome. Both
   mates certify for 95.7 % of pairs.
2. A mate gets a label when it is certified unique (S0/S1, margin 2 edits)
   on one haplotype and at least 2 edits worse on the other, or when it is
   unique on both and both project to the same GRCh38 locus.
3. minibwa, mapping to the same haplotypes, disagreed with 79 labels; those
   mates were dropped.

Result: 31,931,012 mates (80.05 %) have a GRCh38 truth.

*Bias:* mates certified only by the pair certificate (PR, ~7 % of reads) or
proven repeats get no label. So the labelled set holds half of Test A's
share of segmental duplications (2.7 % against 5.5 %), and Test B
understates the difficult regions.

## Accuracy

**Test B (real reads, 31.9 M mates with a truth)**:

| Aligner | in SAM % | at truth % | avoidable errors, MAPQ ≥ 20 / 40 / 60 | MAPQ ≥ 40: placed / misplaced |
|---|---|---|---|---|
| `certa --bio` | 99.689 | 99.614 | 398 / **103** / **31** | 30,526,106 / 4,621 |
| `certa` (default) | 99.947 | 99.861 | 940 / 474 / 325 | 31,751,509 / 4,923 |
| minibwa | 100.000 | 99.901 | 848 / 587 / 389 | 31,879,502 / 8,578 |
| BWA-MEM2 | 100.000 | 99.900 | 1,458 / 988 / 725 | 31,886,149 / 9,858 |
| strobealign | 100.000 | 99.901 | 610 / 439 / 402 | 31,831,409 / 6,331 |
| minimap2 | 99.999 | 99.899 | 1,027 / 784 / 502 | 31,868,976 / 8,236 |

**Test A (simulated, 18.8 M mates with a truth)**:

| Aligner | in SAM % | at truth % | avoidable errors, MAPQ ≥ 20 / 40 / 60 | MAPQ ≥ 40: placed / misplaced |
|---|---|---|---|---|
| `certa --bio` | 98.433 | 96.277 | 3,371 / **472** / **247** | 15,656,510 / 21,465 |
| `certa` (default) | 98.883 | 97.015 | 4,511 / 2,219 / 1,218 | 17,714,036 / 23,804 |
| minibwa | 100.000 | 97.777 | **2,333** / 1,280 / 529 | 18,184,847 / 37,839 |
| BWA-MEM2 | 100.000 | 97.791 | 4,535 / 2,333 / 1,181 | 18,230,390 / 46,776 |
| strobealign | 100.000 | 97.727 | 4,129 / 1,957 / 1,419 | 18,046,383 / 27,136 |
| minimap2 | 99.988 | 97.747 | 2,892 / 1,561 / 874 | 18,103,835 / 32,685 |

- **Easy regions** (GIAB "not in any difficult region"), misplaced with
  MAPQ ≥ 20: `--bio` has 38 on real reads, against 145–368 for the others,
  and 45 on simulated reads, against 148–440.
- **Ties:** among misplaced mates with MAPQ < 20, many ties are broken
  away from the truth: 10,296 real and 112,069 simulated. `--bio` gives
  these MAPQ 0–3 by design, so they count as misplaced at MAPQ ≥ 1.
- **Moved by `--bio`:** on real reads it moved 7,960 mates to the truth
  that the default mode had wrong or absent. It moved 87,072 away, nearly
  all into the fallback (only 6,353 to a wrong locus). The default mode
  certifies more reads.

## Calibration

Real reads: misplaced mates against those expected from MAPQ (Σ 10^(−MAPQ/10)):

| MAPQ band | `--bio` | default | minibwa | BWA-MEM2 |
|---|---|---|---|---|
| 20–29 | 1,124 / 1,283 | 808 / 22 | 3,107 / 62 | 4,209 / 32 |
| 30–39 | 3,100 / 215 | 3,859 / 28 | 978 / 3 | 979 / 2 |
| 40–49 | 785 / 42 | 924 / 7 | 2,054 / 0.5 | 2,959 / 1.4 |
| 50–60 | 3,836 / 33 | 3,999 / 32 | 6,524 / 33 | 6,899 / 32 |

- **20–29:** `--bio`'s posterior MAPQ is calibrated in this band, the only
  aligner whose MAPQ is.
- **30–49:** it is overconfident there, as every aligner is.
- **50–60:** every aligner errs about 120–220× more than its MAPQ states.
  Those errors are almost all unavoidable on GRCh38.
- **Mid-band cause:** the posterior counts one competitor at the horizon
  `YH` (a convention, THEORY §12). Reads whose MAPQ is set by the horizon
  rather than by a second locus are where it errs.

## Soundness

`bio_score` rescored every end-to-end alignment the other aligners reported,
under `--bio`'s model, and compared it with `--bio`'s certified optimum
p1 for that mate:

| | alignments on certified mates | cheaper than p1 | equal | above |
|---|---|---|---|---|
| Test A (5 aligners) | 93,844,773 | **0** | 92,650,025 | 1,194,748 |
| Test B (5 aligners) | 188,731,593 | **0** | 186,633,828 | 2,097,765 |

THEORY §12 defines alignments with no insertion next to a deletion, and
`bio_band` searches exactly that space. Alignments with adjacent gaps are
therefore outside the model and counted apart:

| Aligner | adjacent I/D, Test A / Test B |
|---|---|
| strobealign | 1,760 / 4,237 |
| minibwa | 150 / 381 |
| minimap2 | 130 / 299 |
| CERTA default | 19 / 72 |
| BWA-MEM2 | 0 / 0 |

Four strobealign alignments in Test A have adjacent I/D and are cheaper
than `--bio`'s optimum, by 3–12 Phred. An example is `17M2D3I130M` at 79
against `--bio`'s `13M3D3M4I130M` at 91. They are complex indels
(deletion-insertions) that the model cannot express.

## Cost

| | Test A (20 M reads) | Test B (39.9 M reads) | peak memory |
|---|---|---|---|
| `certa --bio` | 1,206 s | 2,115 s | 22.2 GiB |
| `certa` (default) | 13.5 s | 22.0 s | 21.2 GiB |
| minibwa | 41 s | 69 s | 10.3 GiB |
| BWA-MEM2 | 210 s | 353 s | 33.4 GiB |
| strobealign | 51 s | 78 s | 13.5 GiB |
| minimap2 | 72 s | 131 s | 12.7 GiB |

The timings are indicative only: another project's job ran on cores 64–71
during the runs. `--bio` costs 89–96× the default mode on GRCh38, against
2.2–6.7× on the 20 Mbp simulated genome, and enumerates ~40× more
candidates per read.

## Issues found

- **DRAGMAP's index on galaxy is broken.** Every alignment drifts by one
  base per 60 bases, the FASTA line length, so DRAGMAP placed 0.001 % of
  Test A at the truth. It is left out of this report.
  - `idx/dragmap` was built with `dragen-os --build-hash-table` on
    `idx/ref.fa`.
  - Earlier DRAGMAP timings in this repository (`k_bench.sh`) used this
    index, with output to `/dev/null`, so they measured a broken
    reference.
- **`bio_score`'s first version scored adjacent-gap alignments.** That was
  outside the model; fixed.

## After phase 1: the fallback and the cost

These measurements come from [tools/bio_probe.cpp](../tools/bio_probe.cpp)
on the first 400,000 mate-1 reads of Test B, and from the same build on
Tests A and B. The probe runs `bio_place` per read and logs why it certified
the read or not, and what the read cost.

**Why reads fall back.** 4.35 % of the reads are not certified:
- **Not placed by the default mode either:** 3.07 % of reads. Under an
  end-to-end model these are not recoverable.
- **Placed by the default mode:** 1.28 % of reads, almost all through the
  pair: pair rescue (PR) 0.82 %, an S1 certificate 0.39 %, proven repeats
  0.07 %.
- **All parts already used:** in 69 % of the failures, every one of the 6
  parts the index allows for a 150 bp read is in S, so a larger budget would
  not help. The horizon cannot rise much above 111 Phred, because a gap may
  open for 22 anywhere (the tandem-run discount).
- **Out of reach of a single read:** where the edit certificate holds, the
  failing reads' own alignments cost 57–173 under the model (5th–95th
  percentile, median 111), at or above that horizon.

**Where the time goes** (`perf` on 3,000 of the slowest reads; output-identical
variants timed with the probe):
- **Heavy reads:** 87 % of the time goes to the 11 % of reads with ≥ 1 M DP
  cells, reads in young repeats aligning a median of 1,700 bands each, at
  about 12 µs per band.
- **Mostly confident:** 49 % of all the time goes to heavy reads with
  MAPQ ≥ 30 and p₁ = 0, whose bands are weak copies bounded below p₁ + 60.
- **Inside a band:** the tandem-run scan of the band's reference took
  38 %, the DP 35 %. The edit certificate and enumeration took about 8 %.
  The earlier note that candidate enumeration dominates was wrong.
- **What did not help:**
  - pruning the DP at p₁ + δ instead of H: 7 %;
  - narrowing the bands to the edit alignment's cost + δ: 8 %;
  - enumerating only the rarest parts needed: 0 %, because heavy reads need
    every part that fits.

**What changed** (branch `claude/certa-bio-fallback`):
- **Pair rescue** (THEORY.md Proposition B3, `bio_rescue`). When one mate is
  certified and the other is not, the other is aligned exactly in the
  certified mate's proper-pair windows. It is placed there when that pair is
  proven optimal against every pair not seen: the certified mate at an
  unlisted locus, the other mate outside the window (improper, +40), its gaps
  past the reach, or a window not aligned. Such a mate is written as BM with
  `YR:i:1`. Its YH bounds every alignment of it.
  - The oracle test checks every rescue in every setting against a
    whole-reference brute force.
  - Mutants B15–B17 drop one bound each, and all three are killed.
- **The tandem-run scan** runs on one fetched window, with precomputed gap
  opens. The output is identical on the 400,000 reads; the gain is 1.13×.
- **The tandem-run table** (`--bio-run-table`, on by default; THEORY.md §12).
  `certa map` tabulates every reference base's run length once: 1.4 s on 64
  threads, half a byte a base (1.44 GiB for GRCh38). Lengths are clamped
  where the gap-open cost stops falling (11), so every cost is the same.
  The bands read the table instead of scanning.
  - The test runs every read and rescue with the table and requires the
    scan's results; mutant B18 (clamped one too low) is killed.
- **Data structures on the hot paths**, all with the same output:
  - the band DP clears only the cells written two rows before, not three
    whole rows per row;
  - the cluster bounds are cached in a flat array by part mask, stamped
    per read, instead of a hash map;
  - nearby loci, and the span YA, are found by binary search in sorted
    lists instead of scans over every end position;
  - the rescue's buffers are reused.

**Tests A and B again** with that build ([truth_rerun_bio.sh](truth_rerun_bio.sh);
[report_A_rescue.md](results/truth_2026-10/report_A_rescue.md),
[report_B_rescue.md](results/truth_2026-10/report_B_rescue.md)), scored
against the same truth and phase 1's runs of the other aligners:

| | Test A: phase 1 | Test A: rescue | Test B: phase 1 | Test B: rescue |
|---|---|---|---|---|
| to the fallback (all reads) | 4.35 % | **3.48 %** | 4.81 % | **4.11 %** |
| mates with a truth: in SAM | 98.433 % | **99.286 %** | 99.689 % | **99.867 %** |
| at truth | 96.277 % | **97.118 %** | 99.614 % | **99.790 %** |
| placed only by the pair (YR:i:1), with a truth | – | 161,266 | – | 56,807 |
| of them at truth; misplaced at MAPQ ≥ 20 / 40 | – | 97.6 %; 874 / 12 | – | 98.9 %; 187 / 3 |
| avoidable errors, MAPQ ≥ 20 / 40 / 60 | 3,373 / 472 / 247 | 3,478 / 470 / 242 | 398 / 103 / 31 | 439 / 104 / 31 |
| wall time (64 cores) | 1,206 s | 1,079 s | 2,115 s | 1,903 s |
| peak memory | 21.54 GiB | 21.57 GiB | 22.24 GiB | 22.28 GiB |

With the run table and the data structures (the same records as the rescue
build, byte for byte, in both tests; [timings.tsv](results/truth_2026-10/timings.tsv)):

| | Test A | Test B |
|---|---|---|
| phase 1 | 1,206 s, 21.54 GiB | 2,115 s, 22.24 GiB |
| rescue, each band scanning its window | 1,079 s, 21.57 GiB | 1,903 s, 22.28 GiB |
| rescue, run table | 709 s, 23.01 GiB | 1,335 s, 23.73 GiB |
| rescue, run table and data structures | **677 s**, 23.03 GiB | **1,277 s**, 23.72 GiB |
| CERTA's default mode | 13.5 s | 22.0 s |

The table gives 1.42–1.52× and the data structures another 3–5 %. The
probe's time on the 400,000 reads falls from 1,315 s (phase 1) to 788 s.
Another job ran on cores 64–127 during the two run-table rows, so they may be
slightly slow.

- **Soundness:** 0 violations in both tests, over 94.6 M (A) and 189.9 M (B)
  alignments of the five other aligners on certified mates. For a rescued mate the check is
  against its YH, which bounds every alignment of it.
- **Rescoring:** `bio_score` reproduces YP on every one of the 19.3 M (A)
  and 38.2 M (B) `--bio` records, the 173,327 and 283,425 rescued ones
  included.
- **Remaining gap to the default mode:** 0.13 % of real mates with a truth
  are not in the SAM file, against 0.05 % for the default mode.
- **MAPQ of the anchor:** a certified mate whose partner is rescued takes
  the pair's posterior, which counts one unseen pair at each bound. That
  lowers its MAPQ where the partner costs well above its YH. Net of the
  rescued mates themselves, 7,224 fewer real mates (18,697 simulated) have
  MAPQ ≥ 40, and slightly fewer of those are misplaced (4,609 against
  4,621; 21,313 against 21,465).
- **Cost:** `--bio` is now 50–58× slower than the default mode, against
  89–96× in phase 1.

## Next

1. **Cost.** The band DP over thousands of repeat copies per heavy read now
   dominates. A DP over several bands at once (SIMD lanes, one band each)
   would cut it.
2. **MAPQ.** Replace the one-competitor-at-`YH` convention with a bound
   that does not make horizon-limited reads overconfident.
3. **Complex indels.** Allow an insertion next to a deletion in the model
   (`bio_band`, Lemma B, the oracle test).
4. **Test B labels.** Label PR-tier mates from the pair certificate, so the
   real-read truth covers segmental duplications.
5. **Phase 2.** HG005 (download 1.8 GB), the held-out verdict, with the
   build frozen after any change from 1–4.
