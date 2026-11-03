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
`--bio` 1.66–1.78× faster, for 1.45 GiB more memory. Aligning up to 16
bands at once in SIMD lanes makes it another 2.35–2.45× faster, filling
the lanes 16 bases at a time from the run table another 1.18–1.21×, and
tracing the reported alignment back over a narrowed band another
1.13–1.17×. The output is the same throughout, and `--bio` is now 16–17×
slower than the default mode instead of 89–96×.

**Complete pipelines** ([below](#complete-pipelines-against-the-truth)):
with the fallback reads mapped by minibwa, `--bio` still has the fewest
avoidable confident errors, but by 1.6–3.0× on real reads instead of 4–10×.
minibwa gives MAPQ ≥ 40 to the reads CERTA cannot certify and misplaces
4–8 % of them. On real reads, the default mode plus minibwa is the most
accurate complete pipeline at MAPQ ≥ 40.

**MAPQ of the uncertified mates** ([below](#mapq-of-the-uncertified-mates)):
divided by 4, minibwa's MAPQ on those reads comes within 4× of their error
rates. Each pipeline's confident rows are then those of CERTA alone, and
on real reads `--bio` again has 5–10× fewer avoidable confident errors than
the other aligners. Interval mode with `--iv-mapq cert` (XQ for its
uncertified tiers) matches the corrected default pipeline with no second
aligner.

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
- **Bands in SIMD lanes** (`bio_bands`; THEORY.md §12; branch
  `claude/certa-bio-simd`). Up to 16 bands of a read are aligned at once,
  one 8-bit lane each. Each cell holds its value capped at the band's
  threshold, which keeps every value below the threshold exact. The bands
  are still taken cheapest bound first; the lanes only align the next ones
  ahead of time, and the ends of bands the loop does not reach are dropped.
  - The test compares every end with the one-band DP on every read's bands
    (3.6 M ends per seed). It runs every read and rescue one band at a time
    and requires the same results. It checks both DPs against a plain DP
    without pruning on 400,000 tiny cases.
  - Mutants B19–B23 are all killed. B22, a skipped row reset, survived the
    first version of the test: it changes an end only in a rare geometry.
    That led to the tiny cases.
- **The lanes filled 16 bases at a time** (`bio_lane_setup`, `bio_opens8`;
  THEORY.md §12; branch `claude/certa-bio-window`). A `perf` profile of the
  lanes build put about a quarter of `--bio`'s time in filling the lanes,
  against 9 % in the DP itself. Each band's window was built base by base
  and written into the lanes one byte at a time, 16 bytes apart. Now:
  - the run table's 4-bit runs are unpacked 32 at a time, and each gap open
    is one 16-entry table lookup (`vqtbl1q_u8` on ARM, `pshufb` on x86);
  - each band's window is written as a row, and 16 rows at a time are moved
    into the lanes by a 16 × 16 byte transpose;
  - the lanes read the reference in place instead of a copy of each window.

  The bytes are the same as before. A new test compares them, and
  `bio_window`'s gap opens, with costs computed from scratch (67 M bytes per
  seed). The lanes and the tiny cases are now also checked with the run
  table. Mutants B24–B29 are all killed, and the test is clean under
  AddressSanitizer.
- **The traceback over a narrowed band** (`bio_align_locus`; THEORY.md §12;
  branch `claude/certa-bio-trace`). The reported locus was re-aligned with
  traceback over every cell of its cluster's band, with no threshold (about
  a tenth of the time in the profile). Its penalty V is known by then, so:
  - the band is narrowed to the diagonals within g of the end, g the most
    gap bases a penalty of V allows (Lemma W_B), and pruned at V + 1;
  - a rescued mate's traceback takes the same band in its window, which is
    otherwise over 1,000 diagonals wide.

  THEORY.md proves that this gives the whole band's traceback, ties
  included. The test traces 80,000 ends of random bands both ways, with the
  narrowed band's edges reached about 7,600 times below the end and 2,500
  above, and every locus and rescued mate of the oracle test. Mutants
  B30–B33 are all killed, and the test is clean under AddressSanitizer and
  UndefinedBehaviorSanitizer.

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
| rescue, run table and data structures | 677 s, 23.03 GiB | 1,277 s, 23.72 GiB |
| the same, bands in SIMD lanes | 288 s, 23.04 GiB | 521 s, 23.76 GiB |
| the same, lanes filled 16 bases at a time | 238 s, 23.03 GiB | 443 s, 23.73 GiB |
| the same, traceback over a narrowed band | **211 s**, 23.03 GiB | **379 s**, 23.78 GiB |
| CERTA's default mode | 13.5 s | 22.0 s |

The table gives 1.42–1.52× and the data structures another 3–5 %.
Another job ran on cores 64–127 during the two run-table rows, so they may be
slightly slow.

**The lanes** give 2.35× (A) and 2.45× (B). Both builds were run back to
back, with nothing else running. The build before the lanes took 677 s and
1,277 s again. The SAM output (by md5) and the fallback FASTQ files are
identical in both tests.

**Filling the lanes 16 bases at a time** gives 1.21× (A) and 1.18× (B).
It was run back to back with the lanes build, which took 288 s and 521 s
again. The SAM output and the fallback FASTQ files are again identical in
both tests.

**The traceback over a narrowed band** gives 1.13× (A) and 1.17× (B). It
was run back to back with the previous build, which took 239 s and 443 s
again. The SAM output and the fallback FASTQ files are identical in both
tests.

The probe's time on the 400,000 reads:

| build | all reads | heavy reads (≥ 1 M cells) | edit certificate |
|---|---|---|---|
| phase 1 | 1,315 s | 1,147 s | – |
| run table and data structures | 788 s | 650 s | 108 s |
| SIMD lanes | 315 s | 192 s | 109 s |
| SIMD lanes, `--bio-run-table 0` | 642 s | 499 s | 108 s |
| SIMD lanes, one band at a time | 773 s | 635 s | 108 s |
| SIMD lanes filled 16 bases at a time | 264 s | 145 s | 109 s |
| traceback over a narrowed band | **225 s** | **143 s** | 109 s |

The probe's output is identical in every column but the timings.

- **Heavy reads:** they align 1,991 bands each. Leaving out their edit
  certificate, a band now costs 1.04 µs: 1.06 µs before the narrowed
  traceback, 1.6 µs with the lanes alone, and 7.0 µs before them.
- **The narrowed traceback** saves time on the other reads, which align few
  bands. Their `--bio` time besides the edit certificate fell from 65 s to
  28 s.
- **The edit certificate,** which `--bio` runs first, is now 48 % of the
  time.
- **The run table** matters more with the lanes: without it, filling each
  band's window from a scan takes half the time. Test A takes 628 s without
  it, with the same SAM output, against 288 s (21.60 GiB against 23.04 GiB).
  With the lanes filled 16 bases at a time, which needs the table, it is
  632 s against 238 s (2.65×), again with the same SAM output.

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
- **Cost:** `--bio` is now 16–17× slower than the default mode, against
  89–96× in phase 1 (50–58× before the lanes, 21–24× before filling them 16
  bases at a time, 18–20× before the narrowed traceback).

## Complete pipelines against the truth

Until now CERTA's records were scored alone. Its fallback reads went to
`/dev/null` and counted as absent, while every other aligner was scored on
all its reads. [truth_complete.sh](truth_complete.sh) ran the complete
pipelines on the same reads with the build of e95e619 (`main` at c2e2cd3,
2026-11-03; galaxy, cores 0–63). They are scored against the same truth and
the same phase 1 runs of the other aligners
([report_A_complete.md](results/truth_2026-10/report_A_complete.md),
[report_B_complete.md](results/truth_2026-10/report_B_complete.md)):
- **`certa + minibwa`, `certa --bio + minibwa`:** the default mode or
  `--bio` writes its uncertified pairs to FASTQ, and minibwa maps them as in
  the timed pipeline. That means `-K 10m`, with `drop_certified_mates.awk`
  keeping CERTA's record of a certified mate. That is 536,624 and 447,096
  pairs in Test A, and 1,200,514 and 1,185,464 in Test B.
- **`certa` interval mode:** CERTA places every read itself, with no
  fallback.
- **Checks:**
  - Every mate has exactly one primary record, and no mate is in both
    outputs.
  - The default mode's records are identical to phase 1's.
  - `--bio`'s fallback FASTQ files match the narrowed-traceback benchmark's
    (md5).

Mates with a truth. *Confident:* MAPQ ≥ 40. *Avoidable:* placed elsewhere,
while another of the nine runs in the report placed them at the truth.

| Test B (real) | at truth % | confident % | misplaced per M, MAPQ ≥ 40 / 60 | avoidable per M, MAPQ ≥ 40 / 60 | easy regions, misplaced at MAPQ ≥ 20 | wall (s) |
|---|---|---|---|---|---|---|
| `certa --bio` + minibwa | 99.900 | 95.69 | 219 / 147 | **11.6 / 6.2** | **244** | 379 + 24 |
| `certa` + minibwa | **99.902** | 99.48 | **187** / 120 | 20.2 / 13.4 | 253 | 21 + 26 |
| `certa` interval mode | 99.899 | 99.48 | 211 / 149 | 40.4 / 32.0 | 464 | 32 |
| minibwa | 99.901 | 99.84 | 269 / **107** | 22.5 / 16.1 | 285 | 69 |
| BWA-MEM2 | 99.900 | **99.86** | 309 / 189 | 35.2 / 26.5 | 368 | 353 |
| strobealign | 99.901 | 99.69 | 199 / 173 | 18.1 / 16.9 | 280 | 78 |
| minimap2 | 99.899 | 99.81 | 258 / 188 | 29.0 / 19.9 | 301 | 131 |
| *`certa --bio` alone* | *99.790* | *95.58* | *151 / 112* | *3.5 / 1.1* | *47* | *379* |
| *`certa` alone* | *99.861* | *99.44* | *155 / 106* | *15.7 / 10.9* | *145* | *21* |

| Test A (simulated) | at truth % | confident % | misplaced per M, MAPQ ≥ 40 / 60 | avoidable per M, MAPQ ≥ 40 / 60 | easy regions, misplaced at MAPQ ≥ 20 | wall (s) |
|---|---|---|---|---|---|---|
| `certa --bio` + minibwa | 97.630 | 83.59 | 1,741 / 1,020 | **69 / 48** | 307 | 210 + 15 |
| `certa` + minibwa | 97.699 | 94.70 | 1,552 / 915 | 152 / 85 | **297** | 13 + 24 |
| `certa` interval mode | 97.690 | 94.68 | 1,671 / 1,055 | 234 / 160 | 497 | 20 |
| minibwa | 97.777 | 96.63 | 2,081 / **421** | 80 / 39 | 345 | 41 |
| BWA-MEM2 | **97.791** | **96.87** | 2,566 / 981 | 139 / 75 | 440 | 210 |
| strobealign | 97.727 | 95.90 | **1,504** / 1,182 | 121 / 91 | 348 | 51 |
| minimap2 | 97.747 | 96.20 | 1,805 / 1,043 | 96 / 59 | 334 | 72 |
| *`certa --bio` alone* | *97.118* | *83.10* | *1,363 / 837* | *30 / 21* | *55* | *210* |
| *`certa` alone* | *97.015* | *94.13* | *1,344 / 841* | *127 / 72* | *148* | *13* |

The wall times of the complete pipelines add the two steps run one after
the other: CERTA, then minibwa on its fallback. The other aligners' times
come from phase 1's session.

- **Real reads, default mode plus minibwa.** It is the most accurate
  complete pipeline at MAPQ ≥ 40, with 187 misplaced per million against
  199–309 for the others. It also places the most mates at the truth,
  though the spread is only 0.003 points. At MAPQ ≥ 60, minibwa is lower
  (107 against 120).
- **`--bio` plus minibwa.** It has the fewest avoidable confident errors on
  both tests, but its lead shrinks:
  - real reads, MAPQ ≥ 40: 11.6 per million against 18.1–35.2 (1.6–3.0×),
    where `--bio` alone had 3.5;
  - simulated reads: 69 against 80–234 at MAPQ ≥ 40, but at MAPQ ≥ 60
    minibwa is lower per million (39 against 48);
  - its edge in easy regions mostly goes: 244 misplaced at MAPQ ≥ 20,
    against 280–368 for the others and 47 for `--bio` alone.
- **The fallback's MAPQ is the cause.**
  - minibwa gives MAPQ ≥ 40 to 35,219 of `--bio`'s fallback mates on real
    reads and misplaces 5.9 % of them, where MAPQ 40 claims 0.01 %. On
    simulated reads the rate is 6.6 %. For the default mode's fallback it
    is 8.1 % (real) and 3.6 % (simulated).
  - These mates make up 31 % of `--bio` + minibwa's misplacements at
    MAPQ ≥ 40 on real reads (22 % simulated), and 14–17 % of the default
    pipeline's.
  - So a read CERTA cannot certify is one whose heuristic MAPQ cannot be
    trusted.
  - Correcting that MAPQ: [below](#mapq-of-the-uncertified-mates).
- **Interval mode** needs no second aligner. But the reads it places
  without a certificate carry its confident errors: 40 avoidable per
  million at MAPQ ≥ 40 on real reads, the most of any run.
- **Simulated reads:** every CERTA pipeline places fewer mates at the truth
  than every other aligner (97.63–97.70 % against 97.73–97.79 %).
  strobealign has the lowest misplacement rate at MAPQ ≥ 40 (1,504 per
  million, against 1,552 for CERTA + minibwa).
- **Avoidable counts** depend on which runs are in the report. With these
  nine, `--bio` alone has 106 at MAPQ ≥ 40 on real reads; with phase 1's
  six, it had 103.

## MAPQ of the uncertified mates

A MAPQ rule changes MAPQ only, never a placement. So
[truth_caps.py](truth_caps.py) scores each rule by rewriting the arrays of
the runs above, with no remapping
([caps_A.md](results/truth_2026-10/caps_A.md),
[caps_B.md](results/truth_2026-10/caps_B.md)). For XQ, interval mode was
also run with `--iv-mapq cons` (same build): the same placements, with
lower MAPQ for 36.1 M of the 39.9 M real mates.

The rule for interval mode is now built in as `--iv-mapq cert`: XQ for
tiers R, IV and MR, the point estimate otherwise. Its rows below come from
its own runs (the build of this change):
- they match the rewrite of the other two runs exactly, on both tests;
- the same build writes the same records as e95e619 in the default mode
  and in interval mode with `--iv-mapq point`, on both tests.

- **minibwa's MAPQ on the fallback mates is inflated at every level**, not
  only at 40. On real reads its MAPQ 60 is worth Phred 13–15 (the misplaced
  fraction), 40–59 is worth 5–9, and 20–39 is worth 2–3. On simulated
  reads: 15–17, 6–13 and 3–9.
  - Uncorrected, its MAPQ 40–60 band expects at most one misplaced mate
    per test, where there are 1,011–6,079.
  - Divided by 4, it is at most 3.6× too confident in any band (0, 1–4,
    5–9, 10–14, 15–19) on either test. Capped at 39, it only leaves the
    confident band.
  - `drop_certified_mates.awk -v div=4` divides it in the pipeline. Scoring
    its output gives the same arrays as the rewrite, on both tests.
- **Interval mode's uncertified tiers have the same problem.** On real
  reads it misplaces 21 % of the IV mates it gives MAPQ ≥ 40 (1,535 of
  7,337) and 11 % of the MR mates (166 of 1,506).
  - Their XQ (Proposition Q) is 0 for every IV and MR record.
  - Tier R aligns uncertified reads exactly. Its XQ is not always 0, and at
    XQ ≥ 20 it misplaces 0.4 % of real mates (4 of 1,089) and 0.9 % of
    simulated ones.
  - `--iv-mapq cons` cannot be used as MAPQ by itself: it gives XQ to the
    certified tiers too, which leaves 0.56 % (real) and 3.7 % (simulated)
    of mates at MAPQ ≥ 40.

Mates with a truth; columns as above. *Composed:* `--bio`'s records, plus
interval mode's for the mates `--bio` leaves to the fallback. The other
aligners' rows are in the tables above.

| Test B (real) | at truth % | confident % | misplaced / avoidable per M, MAPQ ≥ 40 | ≥ 60 | easy regions, misplaced at MAPQ ≥ 20 | wall (s) |
|---|---|---|---|---|---|---|
| `certa` + minibwa | 99.902 | 99.48 | 187 / 20.2 | 120 / 13.4 | 253 | 21 + 26 |
| `certa` + minibwa, fallback MAPQ / 4 | 99.902 | 99.44 | 155 / 15.7 | 106 / 10.9 | 145 | 21 + 26 |
| `certa` interval mode | 99.899 | 99.48 | 211 / 40.4 | 149 / 32.0 | 464 | 32 |
| `certa` interval mode, XQ for IV, MR | 99.899 | 99.45 | 158 / 15.9 | 107 / 11.1 | 151 | 32 |
| `certa` interval mode, `--iv-mapq cert` | 99.899 | 99.44 | 155 / 15.4 | 106 / 10.8 | 145 | 32 |
| `certa --bio` + minibwa | 99.900 | 95.69 | 219 / 11.6 | 147 / 6.2 | 244 | 379 + 24 |
| `certa --bio` + minibwa, fallback MAPQ / 4 | 99.900 | 95.58 | 151 / 3.5 | 112 / 1.1 | 47 | 379 + 24 |
| `certa --bio` + interval mode, MAPQ / 4 (composed) | 99.897 | 95.58 | 151 / 3.5 | 112 / 1.1 | 47 | – |
| minibwa | 99.901 | 99.84 | 269 / 22.5 | 107 / 16.1 | 285 | 69 |
| strobealign | 99.901 | 99.69 | 199 / 18.1 | 173 / 16.9 | 280 | 78 |

| Test A (simulated) | at truth % | confident % | misplaced / avoidable per M, MAPQ ≥ 40 | ≥ 60 | easy regions, misplaced at MAPQ ≥ 20 | wall (s) |
|---|---|---|---|---|---|---|
| `certa` + minibwa | 97.699 | 94.70 | 1,552 / 151.8 | 915 / 84.6 | 297 | 13 + 24 |
| `certa` + minibwa, fallback MAPQ / 4 | 97.699 | 94.13 | 1,344 / 126.5 | 841 / 71.8 | 148 | 13 + 24 |
| `certa` interval mode | 97.690 | 94.68 | 1,671 / 233.8 | 1,055 / 159.9 | 497 | 20 |
| `certa` interval mode, XQ for IV, MR | 97.690 | 94.34 | 1,339 / 115.1 | 834 / 65.3 | 153 | 20 |
| `certa` interval mode, `--iv-mapq cert` | 97.690 | 94.13 | 1,316 / 112.7 | 818 / 63.8 | 148 | 20 |
| `certa --bio` + minibwa | 97.630 | 83.59 | 1,741 / 69.2 | 1,020 / 48.2 | 307 | 210 + 15 |
| `certa --bio` + minibwa, fallback MAPQ / 4 | 97.630 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 55 | 210 + 15 |
| `certa --bio` + interval mode, MAPQ / 4 (composed) | 97.611 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 55 | – |
| minibwa | 97.777 | 96.63 | 2,081 / 80.0 | 421 / 38.5 | 345 | 41 |
| strobealign | 97.727 | 95.90 | 1,504 / 121.0 | 1,182 / 90.7 | 348 | 51 |

- **The correction works.** With the fallback's MAPQ divided by 4, or
  capped at 39, each complete pipeline's rows at MAPQ ≥ 40 and ≥ 60 are
  CERTA alone's, exactly.
  - Division also brings the easy-region errors at MAPQ ≥ 20 back to
    CERTA alone's. On real reads that is 253 → 145 for the default mode and
    244 → 47 for `--bio`.
  - The cost is fewer confident mates: 0.04 points (default) and 0.11
    (`--bio`) on real reads, and 0.57 and 0.49 points on simulated ones.
- **`--bio`'s lead is back.** On real reads at MAPQ ≥ 40, it has 3.5
  avoidable errors per million confident mates against 18.1–35.2 for the
  other aligners (5–10×). It has 47 easy-region errors at MAPQ ≥ 20 against
  280–368. On simulated reads it has 30 against 80–139 at MAPQ ≥ 40, and 21
  against 39–91 at MAPQ ≥ 60.
- **The default mode needs no second aligner.**
  - Interval mode with `--iv-mapq cert` matches the confident rows of the
    corrected default + minibwa pipeline, or beats them. On real reads at
    MAPQ ≥ 40 it has 155 misplaced and 15.4 avoidable per million, against
    155 and 15.7. On simulated reads: 1,316 and 112.7, against 1,344 and
    126.5.
  - minibwa places 0.003 points more of the real mates at the truth (about
    1,000 of 31.9 M) and 0.009 points more of the simulated ones. Once its
    MAPQ is divided, all of its records are below MAPQ 20.
  - Interval mode is faster: 32 s against 21 + 26 s on Test B, and 20
    against 13 + 24 on Test A. It needs one index, not two; minibwa's step
    peaks at 7.8–8.6 GiB.
- **Interval mode with `--iv-mapq cert` against the other aligners:**
  - It has the fewest misplaced mates at MAPQ ≥ 40 on both tests: on real
    reads 155 per million against 199–309, and on simulated reads 1,316
    against 1,504–2,566.
  - It calls fewer mates confident: 99.44 % against 99.69–99.86 % on real
    reads, and 94.13 % against 95.90–96.87 % on simulated ones.
  - At MAPQ ≥ 60 it ties minibwa on real reads (106 against 107). On
    simulated reads, minibwa is far lower (421 against 818).
  - Its avoidable errors at MAPQ ≥ 40 are the fewest on real reads (15.4
    against 18.1–35.2), but not on simulated ones (112.7 against 80–139).
- **`--bio` still needs a second aligner**, since it does not combine with
  `--intervals`.
  - The composed rows estimate interval mode as its fallback. They match
    `--bio` + minibwa at MAPQ ≥ 20 and above, and place 0.003 points fewer
    mates at the truth on real reads (0.019 on simulated ones).
  - With XQ instead of the division, the composed run has 154 / 4.1 per
    million at MAPQ ≥ 40 on real reads, and 97 easy-region errors.
  - This is an estimate. The records come from a separate run, and nothing
    was timed. Placing `--bio`'s uncertified reads with interval mode's
    tiers in the same process would remove minibwa from `--bio` too.
- **What was fitted.** The division by 4 was fitted on these two HG002
  tests, the calibration sample; the held-out check is phase 2. XQ is a
  proven bound, with nothing fitted.

## Next

1. **Cost.** The edit certificate that `--bio` runs first is now 48 % of
   its time; its band DP (`band_distance_vec8`) is also the largest cost of
   `--full-cert`.
2. **MAPQ of uncertified records.**
   - Test `--iv-mapq cert` in variant calling before making it the
     default: XQ for every record dropped the reads that carry variants
     (FINDINGS.md, Aim 3).
   - Place the reads `--bio` cannot certify with interval mode's tiers in
     the same process, so that `--bio` needs no second aligner either.
3. **MAPQ.** Replace the one-competitor-at-`YH` convention with a bound
   that does not make horizon-limited reads overconfident.
4. **Complex indels.** Allow an insertion next to a deletion in the model
   (`bio_band`, Lemma B, the oracle test).
5. **Test B labels.** Label PR-tier mates from the pair certificate, so the
   real-read truth covers segmental duplications.
6. **Phase 2.** HG005 (download 1.8 GB), the held-out verdict, with the
   build frozen after any change from 1–5.
