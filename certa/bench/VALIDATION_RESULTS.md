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
- **The confident placements are not more accurate than minibwa's**
  ([below](#matched-sensitivity-and-the-locus-window), 2026-11-04). Phase 1
  reported 4–21× fewer *avoidable* errors for `--bio` on real reads (1.4–4×
  simulated): misplaced mates that another aligner placed at the truth.
  Two choices made that:
  - a 20 bp window, which counts a read soft-clipped at an insertion as
    misplaced;
  - counting each aligner at its own MAPQ ≥ 40, which credits `--bio` for
    calling about 4 % fewer real mates confident (14 % simulated).

  At the locus and the same number of confident mates, minibwa misplaces
  3× fewer real mates than `--bio` (46 against 138 per million) and 5×
  fewer simulated ones (262 against 1,309).
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
  is not their orthologous locus. At the locus this bounds every aligner:
  of minibwa's 89 misplaced per million confident real mates, 0.6 are
  mates another GRCh38 aligner places at the truth.

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
CERTA with its fallback reads mapped by minibwa, and interval mode, which
places every read itself. minibwa gives MAPQ ≥ 40 to the reads CERTA cannot
certify and misplaces 4–8 % of them.

**MAPQ of the uncertified mates** ([below](#mapq-of-the-uncertified-mates)):
divided by 4, minibwa's MAPQ on those reads comes within 4× of their error
rates, and each pipeline's confident rows are then those of CERTA alone.
Interval mode with `--iv-mapq cert` (XQ for its uncertified tiers) matches
the corrected default pipeline with no second aligner.

**Matched sensitivity and the locus window**
([below](#matched-sensitivity-and-the-locus-window)): at the locus and at
the same number of confident mates, minibwa is the most accurate aligner
measured on both tests. Interval mode misplaces 1.4× as many confident
real mates (2.9× simulated), and `--bio` 3.0× (5.0×). Nearly all the
remaining errors are on mates no GRCh38 run places at the truth.
T2T-CHM13 places 38 % of those real mates there; a pangenome mapper (vg
giraffe on HPRC v2.1, HG002 left out) is not more accurate than minibwa.

**A MAPQ from T2T-CHM13** ([below](#a-mapq-from-t2t-chm13)): setting MAPQ
0 where BWA-MEM2 on T2T-CHM13 confidently places the read elsewhere
removes 42–59 % of every aligner's confident errors on real reads (minibwa
89 → 52 per million, interval mode 122 → 65), and 17–22 % on simulated
ones. A second GRCh38 aligner in its place removes 0–3 %. It helps every
aligner alike, so minibwa stays ahead.

**Re-paired mates and a fitted MAPQ**
([below](#re-paired-mates-and-a-fitted-mapq)): two opt-in options of paired
interval mode. `--iv-pair all` re-pairs certified mates of improper pairs;
`--mapq cal` writes a MAPQ fitted on HG002 from the gap to the second-best
locus, the edits and the clusters verified. Held out on HG005
([below](#held-out-hg005)), they make CERTA the most accurate aligner
measured at operating points of 98 % of mates and below: 22–33 misplaced per
million against minibwa's 37 and 108–110 for the others (14–19 against 22
with the T2T veto on both). At 99.3 % and above minibwa stays ahead (80
against 86).

**A MAPQ from both mates** ([below](#a-mapq-from-both-mates)): two more
opt-in options. `--mapq pair` writes a MAPQ from boosted trees over both
mates' records, fitted on HG002. `--t2t` maps every pair on T2T-CHM13 too,
in the same pass, and the trees read where that placement lifts to. Held
out on HG005, CERTA is the most accurate aligner measured at every
operating point. With one index it misplaces 186, 60, 20 and 9 per million
of the 99.75, 99.29, 98 and 95 % most confident mates, against minibwa's
208, 80, 37 and 37. With `--t2t` it is 93, 25, 7 and 3, ahead of minibwa
with BWA-MEM2's T2T veto too (167, 48, 22, 22). Placements are unchanged,
and CERTA's MAPQ equals the fit's replay on every mate.

**Variant calling** ([below](#variant-calling)): DeepVariant on 30× HG002
and HG005 (chr20 and the CMRG genes) makes 1.8–5.7× more SNV errors on
CERTA's records with the fitted MAPQs than on minibwa's (HG002:
3,214–6,052 against 1,745; HG005: 1,253–3,099 against 543).
- **The cause is the MAPQ, not the placements:** the fitted MAPQ is lower
  on reads that show the sample's variants, and refitting does not remove
  that.
- **With BWA-MEM's MAPQ formula** (the default), the same placements give
  12 % (HG002) and 17 % (HG005) more SNV errors than minibwa.
- **That remaining gap** is mostly reads CERTA does not place: chimeric,
  very divergent, or with long deletions.

So the fitted MAPQs are not fit for variant calling. *(2026-11-19: they
are, written on BWA-MEM's scale; see below.)*

**Smaller indexes** ([below](#smaller-indexes)): the key bits stored per
entry can be read from the reference instead. `certa index --from X.cidx
--key-bits 8` (or `0`) gives byte-identical output with a 17.6 GiB (or
14.9 GiB) GRCh38 index instead of 20.4 GiB, at 21 % (35 %) more CPU time.

**How far from the limit** ([below](#how-far-from-the-limit), 2026-11-19):
- **On GRCh38, minibwa's placement is at the limit.** An oracle that scores
  each misplaced mate at its truth finds that minibwa's search misses
  0.05 % of simulated mates and 0.001 % of real ones. Everything else it
  misplaces scores better somewhere else in GRCh38 (or ties), so no
  aligner that maximizes the alignment score there can place it right.
- **The limit itself is far away.** Mapped to the sample's own genome,
  3 per million of the 96 % most confident simulated mates are misplaced,
  against 1,702 for minibwa on GRCh38. That information is not in GRCh38.
- **CERTA's search was short of the GRCh38 limit** (2.2× minibwa's misses
  on simulated reads), and its confident errors are a matter of MAPQ.

**The fitted MAPQ on BWA's scale** ([below](#the-fitted-mapq-on-bwas-scale)):
the fitted MAPQ is a Phred, while DeepVariant was trained on BWA-MEM's
MAPQ, which means far less by the same numbers. Written on BWA-MEM's scale
(`--mapq-scale bwa`, now the default with `--mapq pair`), it gives fewer
DeepVariant errors than CERTA's default MAPQ, and with `--t2t` as many SNV
errors as minibwa on HG002 (1,748 against 1,745; CMRG 791 against 800).

**Joint seeds of a pair** ([below](#joint-seeds-of-a-pair)): both mates'
seeds find the copy of a repeat that fits the pair when each mate alone
fits other copies better. CERTA's search misses fall below minibwa's on
both tests, and it places more mates at their locus than minibwa and
BWA-MEM2, for 45 % more CPU (about minibwa's).

**The default build, all together** ([below](#the-default-build-all-together)):
on HG005 (held out) CERTA leaves the fewest mates away from their locus,
and with `--mapq pair --t2t` misplaces 2–12× fewer confident mates than
minibwa. In DeepVariant it makes fewer SNV errors than minibwa and
BWA-MEM2 on HG002 (1,696 against 1,745 and 1,761) and 2 % more than
minibwa on HG005 (555 against 543); indels are level on HG002 and behind
on HG005 (77 against 63). That configuration costs twice minibwa's CPU
time and 4.6× its memory.

**The pair model refit** ([below](#the-pair-model-refit-on-the-default-builds-placements), 2026-11-20):
both trees fitted again on the default build's records.
- **Held out on HG005:** `--mapq pair` misplaces fewer of its confident
  mates than minibwa at every operating point (185 against 208 per million
  at 99.75 %), at about minibwa's CPU time. With `--t2t` it misplaces
  2.2–14× fewer.
- **DeepVariant with `--t2t`:** fewer SNV errors than with minibwa or
  BWA-MEM2 on both genomes (HG002 1,671 against 1,745 and 1,761; HG005 534
  against 543 and 553). Indels stay behind minibwa (HG005 74 against 63).
  Most of that indel gap comes from BWA-MEM's gap costs, which CERTA
  shares, and a few cases are alignments its search misses.
- **XM** (the margin of a pair placed from the alignments found) did not
  help the trees.

**Tier R's wide band** ([below](#tier-rs-wide-band), 2026-11-20, `--iv-wide 1`):
in a repeat, tier R aligned the two sides of an indel longer than 16 bases
in separate bands and clipped one side. One band of 95 diagonals around a
clipped alignment aligns such reads through the indel, for 0.3 % more CPU.
- **Variant calls:** DeepVariant's indel errors fall by 13–25 on HG002 and
  by 1–3 on HG005, and its SNV errors by 8–15 on HG002.
- **With `--t2t`:** CERTA's total errors are below minibwa's on HG002 (3,446
  against 3,520) and equal on HG005 (606).
- **Truth sets:** they barely change.

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

*At truth* in this section and the next ones means within 20 bp of the
truth, and each aligner is counted at its own MAPQ thresholds. Both
overstate `--bio`'s lead, which disappears at the locus and at the same
number of confident mates ([below](#matched-sensitivity-and-the-locus-window)).

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
*At truth:* within 20 bp, each run at its own MAPQ; at the locus and at a
matched number of confident mates the ranking changes
([below](#matched-sensitivity-and-the-locus-window)).

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

- **Real reads, default mode plus minibwa.** At its own MAPQ ≥ 40 it has
  the fewest misplaced mates per million, 187 against 199–309, but it calls
  fewer mates confident (99.48 % against 99.69–99.86 %). It places the most
  mates at the truth, by only 0.003 points. At MAPQ ≥ 60, minibwa is lower
  (107 against 120).
- **`--bio` plus minibwa.** In these tables it has the fewest avoidable
  confident errors on both tests, though more than `--bio` alone. That is
  an artifact of the 20 bp window and of calling 4 % fewer mates confident
  ([below](#matched-sensitivity-and-the-locus-window)):
  - real reads, MAPQ ≥ 40: 11.6 per million against 18.1–35.2,
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

Mates with a truth; columns as above (20 bp, each run at its own MAPQ).
*Composed:* `--bio`'s records, plus
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
- **`--bio` at its own MAPQ ≥ 40** has 3.5 avoidable errors per million
  confident mates on real reads, against 18.1–35.2 for the other aligners,
  and 47 easy-region errors at MAPQ ≥ 20 against 280–368. On simulated
  reads it has 30 against 80–139 at MAPQ ≥ 40, and 21 against 39–91 at
  MAPQ ≥ 60. This is not a lead in accuracy: at the locus and at the same
  number of confident mates, minibwa misplaces fewer
  ([below](#matched-sensitivity-and-the-locus-window)).
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
  - At its own MAPQ ≥ 40 it has fewer misplaced mates per million than
    the others on both tests: on real reads 155 against 199–309, and on
    simulated reads 1,316 against 1,504–2,566.
  - But it calls fewer mates confident: 99.44 % against 99.69–99.86 % on
    real reads, and 94.13 % against 95.90–96.87 % on simulated ones. At the
    same number of confident mates it is not ahead
    ([below](#matched-sensitivity-and-the-locus-window)).
  - At MAPQ ≥ 60 it ties minibwa on real reads (106 against 107). On
    simulated reads, minibwa is far lower (421 against 818).
  - In this table its avoidable errors at MAPQ ≥ 40 are the fewest on real
    reads (15.4 against 18.1–35.2), but not on simulated ones (112.7
    against 80–139).
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

## Matched sensitivity and the locus window

Two choices in the comparisons above favour CERTA, and `--bio` most
(2026-11-04):
- **The 20 bp window.** A mate is at the truth when its unclipped leftmost
  base is within 20 bp of the truth's. That counts mates placed at their
  locus but written differently as misplaced. On real reads, each other
  aligner has 2,022–2,504 confident mates (MAPQ ≥ 40) 21–150 bp from the
  truth; interval mode has 1,055 and `--bio` 399.
  - About half of the other aligners' are soft-clipped by 10 bases or more,
    and their offset is the clip's length (median offset less the clip −1
    to 6 bases, on both tests). The read starts in sequence GRCh38 lacks:
    they clip it, where CERTA writes an indel or `--bio` aligns it end to
    end.
  - The rest have no long clip or indel; their median offset is 30–33 bp.
  - *At the locus* below means within 150 bp (`truth_eval.py score
    --window 150`). On Test B, 1,000 bp gives the same ranking.
- **Each run at its own MAPQ.** Errors per million confident mates credit a
  run that calls fewer mates confident: `--bio` calls 95.58 % of real mates
  confident, the others 99.44–99.86 %. Below, every run is taken at the
  same N: its N most confident mates with a truth, ties at the boundary
  MAPQ counted pro rata.

[truth_headroom.sh](truth_headroom.sh) rescored every run of the tables
above at 150 bp, with no remapping ([truth_headroom.py](truth_headroom.py),
[truth_nearmiss.py](truth_nearmiss.py); reports
[headroom_A.md](results/truth_2026-10/headroom_A.md),
[headroom_B.md](results/truth_2026-10/headroom_B.md)). The minibwa
pipelines have their fallback MAPQ divided by 4. Misplaced / avoidable per
million among each run's N most confident mates with a truth; N is the
share interval mode with `--iv-mapq cert`, or `--bio` + minibwa, calls
confident:

| Test B (real) | 20 bp, N = 99.44 % | locus, N = 99.44 % | 20 bp, N = 95.58 % | locus, N = 95.58 % |
|---|---|---|---|---|
| `certa` interval mode, `--iv-mapq cert` | 155 / 15.4 | 122 / 8.9 | 106 / 10.8 | 77 / 5.5 |
| `certa` + minibwa | 156 / 15.9 | 123 / 9.5 | 106 / 10.9 | 77 / 5.7 |
| `certa --bio` + minibwa | 394 / 35.7 | 332 / 24.5 | 151 / 3.5 | 138 / 1.3 |
| minibwa | 154 / 17.0 | **89 / 0.6** | 107 / 16.1 | **46 / 0.3** |
| BWA-MEM2 | 189 / 26.5 | 114 / 0.7 | 189 / 26.5 | 114 / 0.7 |
| strobealign | 173 / 16.9 | 111 / 3.6 | 173 / 16.9 | 111 / 3.6 |
| minimap2 | 198 / 20.6 | 130 / 4.5 | 188 / 19.9 | 121 / 4.0 |

| Test A (simulated) | 20 bp, N = 94.13 % | locus, N = 94.13 % | 20 bp, N = 83.10 % | locus, N = 83.10 % |
|---|---|---|---|---|
| `certa` interval mode, `--iv-mapq cert` | 1,316 / 112.7 | 1,232 / 95.5 | 818 / 63.8 | 762 / 57.4 |
| `certa` + minibwa | 1,349 / 127.7 | 1,266 / 111.0 | 841 / 71.8 | 785 / 65.6 |
| `certa --bio` + minibwa | 4,349 / 234.5 | 4,245 / 207.4 | 1,363 / 30.4 | 1,309 / 13.2 |
| minibwa | 607 / 42.7 | **430 / 3.3** | 421 / 38.5 | **262 / 2.2** |
| BWA-MEM2 | 981 / 74.8 | 767 / 3.4 | 981 / 74.8 | 767 / 3.4 |
| strobealign | 1,182 / 90.7 | 995 / 52.5 | 1,182 / 90.7 | 995 / 52.5 |
| minimap2 | 1,043 / 59.0 | 831 / 13.1 | 1,043 / 59.0 | 831 / 13.1 |

- **minibwa is the most accurate aligner measured**, at both N, on both
  tests, at the locus.
  - Real reads: 89 misplaced per million against 122 for interval mode
    (1.4×) and 111–130 for the others. At `--bio`'s N it has 46 against
    `--bio`'s 138 (3.0×) and interval mode's 77.
  - Simulated reads: 430 against 1,232 for interval mode (2.9×), and 262
    against `--bio`'s 1,309 (5.0×).
- **CERTA's earlier lead came from the two choices.** At 20 bp, interval
  mode ties minibwa on real reads (155 against 154). minibwa's avoidable
  errors fall from 17.0 to 0.6 per million at the locus: they were near
  misses that CERTA wrote with an indel. `--bio`'s avoidable errors (1.3)
  are no longer the fewest: minibwa has 0.3.
- **Nearly every confident error at the locus is caused by the
  reference.** No GRCh38 run places 16,639 real mates (0.052 %) and
  296,586 simulated ones (1.58 %) at the truth: their origin is missing
  from GRCh38 or has diverged there. Of minibwa's 89 misplaced per million
  on real reads, 0.6 are mates another GRCh38 run places at the truth (the
  other heuristic aligners 0.7–4.5, CERTA 8.9–9.5). So a better algorithm
  on GRCh38 cannot be much more accurate than minibwa. What separates the
  aligners is which of the other mates they call confident, and minibwa's
  MAPQ ranks more of them low than CERTA's does. CERTA's certificate
  proves a read unique in GRCh38, not that it came from there.
- **Interval mode with `--iv-mapq cert` still matches the corrected
  default pipeline** (122 against 123 per million on real reads, 1,232
  against 1,266 simulated), with no second aligner.

**Other references** (Test B, at the locus). Of the 16,639 mates no GRCh38
run places at the truth:

| run | at the truth | at MAPQ ≥ 40 | misplaced at MAPQ ≥ 40 | misplaced / avoidable per M, N = 99.44 % |
|---|---|---|---|---|
| BWA-MEM2 on T2T-CHM13, lifted to GRCh38 | 6,372 (38.3 %) | 5,984 | 7,329 | 782 / 578.7 |
| vg giraffe, HPRC v2.1 graph without HG002 | 4,636 (27.9 %) | 1,420 | 2,529 | 94 / 14.7 |

- **T2T-CHM13 (hs1).** BWA-MEM2 on hs1, lifted through UCSC's
  `hs1ToHg38` chain by [lift_sam.py](lift_sam.py). It places 38 % of those
  mates at the truth, nearly all at MAPQ ≥ 40, so these reads do identify
  their origin given a reference that has it. As a GRCh38 aligner it is
  not usable: the lift leaves 2.09 M of its 39.8 M mapped records (5.3 %)
  without a GRCh38 position and misplaces others (698 easy-region errors
  at MAPQ ≥ 20, against 8–35 for the GRCh38 runs). That limit is the lift's, not T2T's; a truth
  in hs1 coordinates would measure T2T itself.
- **A pangenome.** vg giraffe 1.78.0 on `hprc-v2.1-mc-grch38-eval`, the
  evaluation graph that leaves out HG002, HG005 and NA19240 (d46
  filtered), projected onto GRCh38 with `--ref-name GRCh38`, without
  haplotype sampling. The
  published distance index is format v3, so `vg autoindex -w sr-giraffe`
  rebuilt it (58 min, 104 GiB peak, on a2). Mapping took 641 s on a2's 64
  cores (x86-64; not comparable to galaxy's times), at 69.7 GiB peak.
  - At N = 99.44 % it misplaces 94 per million, against minibwa's 89; at
    N = 98.65 % (minibwa's share at MAPQ ≥ 60), 73 against 46.
  - It leaves 20,998 mates with a truth unmapped (0.066 %), of which
    GRCh38 aligners place at least 19,622 at the truth, and has 169
    easy-region errors at MAPQ ≥ 20 (GRCh38 runs: 8–35).
  - It places 28 % of the mates GRCh38 cannot, but most at low MAPQ.

So on GRCh38 the remaining confident errors are the reference's. Much
better accuracy needs information GRCh38 lacks:
for instance, mapping to both GRCh38 and T2T and lowering the MAPQ of a
read whose confident T2T placement lifts to another locus (next section).

## A MAPQ from T2T-CHM13

The rule was fixed before the first evaluation (2026-11-04;
[truth_t2t.py](truth_t2t.py), [truth_t2t.sh](truth_t2t.sh); reports
[t2t_A.md](results/truth_2026-10/t2t_A.md),
[t2t_B.md](results/truth_2026-10/t2t_B.md)):
- BWA-MEM2 2.2.1 maps the same reads to T2T-CHM13 (hs1), and its records
  are lifted to GRCh38 as above.
- When the T2T record has MAPQ ≥ 40 and lifts more than 150 bp from a
  GRCh38 aligner's placement, or to another contig or strand, that
  placement's MAPQ becomes 0. No placement changes.

Misplaced per million, at the locus and at the same N as in the section
above, without → with the rule:

| run | Test B, N = 99.44 % | Test B, N = 95.58 % | Test A, N = 94.13 % | Test A, N = 83.10 % |
|---|---|---|---|---|
| `certa` interval mode, `--iv-mapq cert` | 122 → 65 | 77 → 41 | 1,232 → 1,023 | 762 → 615 |
| `certa` + minibwa | 123 → 65 | 77 → 41 | 1,266 → 1,050 | 785 → 634 |
| `certa --bio` + minibwa | 332 → 137 | 138 → 72 | 4,245 → 3,536 | 1,309 → 1,056 |
| minibwa | 89 → **52** | 46 → **27** | 430 → **354** | 262 → **209** |
| BWA-MEM2 | 114 → 63 | 114 → 63 | 767 → 597 | 767 → 597 |
| strobealign | 111 → 60 | 111 → 60 | 995 → 812 | 995 → 812 |
| minimap2 | 130 → 71 | 121 → 64 | 831 → 648 | 831 → 648 |

- **Real reads: about half the confident errors go** (42–59 % at
  N = 99.44 %), for 0.04–0.06 points of confident mates. On simulated
  reads it is 17–22 % (N = 94.13 %), for 0.09–0.15 points.
- **The gain comes from T2T's sequence, not from a second aligner.** With
  BWA-MEM2 on GRCh38 as the second opinion (`truth_t2t.py eval --second
  bwa-mem2`), the same rule removes 0–2 % on real reads (minibwa 89 → 89,
  interval mode 122 → 121) and 0–3 % on simulated ones (1,232 → 1,200).
- **Most of the mates it sets to 0 were placed right.** Of the confident
  ones, 14–23 % were misplaced on real reads and 24–34 % on simulated
  ones. They are few (0.04–0.06 % of real mates), so the cost is small,
  but the flags are far from precise. How many of the false ones come from
  the lift, where the chain and the truth's projection (dipcall) align
  repeats differently, was not separated.
- **It does not change the order.** With the rule, minibwa (52 per million
  on real reads) stays ahead of every CERTA pipeline (65), and at `--bio`'s
  N it has 27 against `--bio`'s 72.
- **Variants** (the grids in the reports) are mixed:
  - a T2T threshold of 20 removes a few more errors (minibwa 51 on real
    reads, 324 on simulated);
  - counting T2T records that do not lift as disagreeing costs about 0.3
    points more of confident mates and helps less on real reads;
  - raising the MAPQ of placements T2T confirms adds errors for every
    aligner on real reads.
- **What was fitted:** only the threshold and the window, set beforehand;
  CHM13 is another individual. The truth is still HG002's, so a held-out
  sample needs truth labels for HG001 or HG005.
- **Cost:** BWA-MEM2 on hs1 took 237 s (A) and 411 s (B) on galaxy's 64
  cores, at 38.2 and 37.4 GiB: more than any GRCh38 aligner here, BWA-MEM2
  included. CERTA on hs1 would be the fast second opinion; not tested.

So a reference with the copies GRCh38 lacks is worth about a factor of 2
in confident errors on real reads, for every aligner alike. It is not an
advantage of CERTA's. Larger gains need flags that are more precise and
catch more: a certified T2T placement (CERTA on hs1), or more haplotypes
than CHM13's.

## Re-paired mates and a fitted MAPQ

Two changes aimed at what the matched comparison found (2026-11-04;
[truth_cal.sh](truth_cal.sh), [truth_mapq_fit.py](truth_mapq_fit.py);
reports [cal_A.md](results/truth_2026-10/cal_A.md),
[cal_B.md](results/truth_2026-10/cal_B.md)). Both are opt-in options of
paired interval mode; the default mode is unchanged (same records on
test A), and so is interval mode without them (same records as
`--iv-mapq cert` on both tests).

**`--iv-pair all`: certified mates of an improper pair are re-paired.**
On simulated reads 82 % of interval mode's avoidable confident errors
(1,360 of 1,651) were mates of pairs that were not proper, both mates
certified: each sat at its single-end optimum, while the pair pointed at
the truth. Interval mode already re-paired pairs with an uncertified mate
the way BWA-MEM does (`place_pair`); `--iv-pair all` sends these pairs
through it too. A mate's best alignment in the other's proper-pair window
is taken when the pair beats the two single-end bests minus the unpaired
penalty; the moved mate becomes MR (XQ 0, and its single-end certificate
stays in XL/XU). Otherwise the records are the certified ones.
- Avoidable errors per million at interval mode's N fall from 8.7 to 4.3
  on real reads and from 93 to 34 on simulated ones. The mates at the
  locus rise from 99.908 % to 99.911 % (real) and from 97.728 % to 97.788 %
  (simulated).
- The test re-pairs 14–20 improper certified pairs per seed and checks that
  every re-paired pair is proper, every moved mate is MR with XQ 0, the
  anchor keeps its locus, tier and MAPQ, and every claim passes the checker.

**`--mapq cal`: MAPQ fitted at the locus.** BWA-MEM's formula, which
interval mode uses, reaches 60 once the second locus is about three
mismatches away. On real reads the error rate keeps falling with the gap
over its whole range. For S0 mates at MAPQ ≥ 40 it is:

| gap to the second-best locus (edits) | 2 | 3 | 4 | 5 | 6 | 8 | 10 | none within 10 |
|---|---|---|---|---|---|---|---|---|
| misplaced per million | 10,244 | 5,413 | 3,654 | 2,322 | 1,474 | 587 | 423 | 25 |

A nearby paralog in GRCh38 marks the reads whose own copy may be missing
from it. The edits of the best alignment and the clusters verified add to
that: among S1 mates with no second locus within 10 edits, 47 per million
are misplaced with one edit and one cluster, against 1,242 with two edits
and two to four clusters. `--mapq cal` writes as MAPQ the Phred of the
misplaced share of such records on HG002 (`cal_mapq`,
[mapq_table.h](../include/certa/mapq_table.h)):
- records with MAPQ ≥ 40 in tiers S0, S1 and PR, by tier, edits (0–5), the
  gap to the second-best locus (0–10 edits, or none; `--mapq cal` sets
  `--mapq-limit 10`) and clusters verified (1, 2–4, 5–16, 17 or more);
- every other record by tier and MAPQ;
- cells with fewer than 1,000 mates take their parent's value.

The table was fitted on tests A and B pooled (HG002, the calibration
sample). CERTA's output equals the fit's replay on every mate of both
tests, and the placements are those of `--iv-pair all`. Its values are
Phred-scaled error rates, so they are lower than the formula's: 74.5 % of
real mates get MAPQ ≥ 40 (most S1 reads get 39), against 99.4 % before.

Misplaced per million at matched sensitivity, at the locus; *veto*: the
T2T rule of the previous section, with BWA-MEM2's or CERTA's placements on
hs1:

| run | B, N = 99.44 % | B, N = 95.58 % | A, N = 94.13 % | A, N = 83.10 % |
|---|---|---|---|---|
| `certa` interval mode, `--iv-mapq cert` | 122 | 77 | 1,232 | 762 |
| the same, `--iv-pair all --mapq cal` | 100 | **34** | 1,063 | 267 |
| the same, with BWA-MEM2's T2T veto | 58 | **21** | 894 | 224 |
| the same, with CERTA's T2T veto | 70 | 23 | 966 | 239 |
| minibwa | 89 | 46 | 430 | 262 |
| minibwa, with BWA-MEM2's T2T veto | **52** | 27 | **354** | **209** |
| BWA-MEM2 / strobealign / minimap2 | 114 / 111 / 130 | 114 / 111 / 121 | 767 / 995 / 831 | 767 / 995 / 831 |

- **Real reads:** with both changes and BWA-MEM2's T2T veto, CERTA has
  fewer confident errors than every aligner as shipped at both N (58 against
  89–130; 21 against 46–121). minibwa with the same veto is ahead at the
  looser N (52 against 58) and behind at the stricter one (27 against 21).
- **Simulated reads:** at the stricter N, CERTA with the veto is ahead of
  every aligner as shipped (224 against 262–995) but behind minibwa with
  the veto (209). At the looser N it is behind minibwa as shipped (894
  against 430). The simulated truth keeps the segmental duplications that
  test B's labels under-sample (5.5 % against 2.7 % of mates).
- **CERTA on hs1 as the second opinion** (interval mode, `--iv-pair all`)
  takes 29 s for test B against BWA-MEM2's 411 s, and keeps 70–80 % of its
  gain. Its flags are not more precise: 13–22 % of the confident mates it
  sets to 0 on real reads were misplaced, against 14–23 % for BWA-MEM2's
  (for the fitted MAPQ, whose confident mates are fewer, 6 % and 5 %).
  Counting only certified hs1 placements changes nothing, since with
  `--iv-mapq cert` every hs1 record at MAPQ ≥ 40 is S0, S1 or PR.
- **What was fitted:** the table, on HG002 (both tests). HG005 is the
  held-out check ([below](#held-out-hg005)).

## Held out: HG005

The configurations above were fixed on HG002 and then run once on HG005
(2026-11-04; [truth_hg005.sh](truth_hg005.sh); reports
[t2t_H.md](results/truth_2026-11_hg005/t2t_H.md),
[curve_H.md](results/truth_2026-11_hg005/curve_H.md),
[labels_H.txt](results/truth_2026-11_hg005/labels_H.txt)). The build is
this section's commit
([frozen_source.sha](results/truth_2026-11_hg005/frozen_source.sha)).
- **Truth:**
  - HG005's HPRC release 2 haplotypes (trio-phased hifiasm,
    chromosome-level), aligned to GRCh38 with dipcall's minimap2 settings;
  - the first 20 M NovaSeq 6000 pairs (2×151) of SRR14724528, trimmed as
    test B's;
  - labels made as test B's: certified unique on one haplotype, with
    minibwa's check (52 disagreements).

  30,127,192 mates (75.4 %) have a truth, in the same strata as test B
  (2.7 % in segmental duplications, 34.9 % in difficult regions).

Misplaced per million at the locus among each run's N most confident
mates with a truth (N = 99.29 % is interval mode's confident share), and
each run's cost on galaxy's 64 cores (CPU seconds as user + sys; peak
memory). The veto adds BWA-MEM2 on hs1, a second run (21,819 s at
38.3 GiB), so its rows add the times and take the larger memory:

| run | 99.75 % | 99.29 % | 98 % | 95 % | CPU | memory |
|---|---|---|---|---|---|---|
| `certa` interval mode, `--iv-mapq cert` | 363 | 111 | 65 | 65 | **2,124 s** | 24.6 GiB |
| the same, `--iv-pair all --mapq cal` | 319 | 86 | **33** | **22** | 2,149 s | 24.6 GiB |
| minibwa | 208 | **80** | 37 | 37 | 3,741 s | **10.1 GiB** |
| BWA-MEM2 | **202** | 108 | 108 | 108 | 17,300 s | 32.9 GiB |
| strobealign | 323 | 108 | 108 | 108 | 4,000 s | 24.2 GiB |
| minimap2 | 252 | 152 | 110 | 110 | 6,477 s | 12.7 GiB |
| *with BWA-MEM2's T2T veto:* `--iv-pair all --mapq cal` | 185 | 50 | **19** | **14** | 23,968 s | 38.3 GiB |
| *with the veto:* minibwa | 167 | **48** | 22 | 22 | 25,560 s | 38.3 GiB |
| *with the veto:* BWA-MEM2 | **160** | 60 | 60 | 60 | 39,119 s | 38.3 GiB |

- **The HG002 result holds.** The two changes cut interval mode's
  confident errors by 23 % at N = 99.29 % (111 → 86; on test B 122 → 100
  at 99.44 %) and by 66 % at 95 % (65 → 22).
- **At 98 % and below, CERTA is the most accurate aligner measured:** 22–33
  per million against minibwa's 37 and 108–110 for the others, so 1.1–1.7×
  fewer than minibwa and 3–5× fewer than the rest. With the T2T veto on
  both, it is 14–19 against minibwa's 22.
- **At the loosest operating points minibwa is ahead:** 80 against 86 at
  99.29 %, and 208 against 319 at 99.75 %. CERTA places 99.904 % of the
  mates at the locus, minibwa 99.906 %.
- **Why the curves flatten:** minibwa, BWA-MEM2 and strobealign give MAPQ 60
  to more than 98 % of the mates, so below that their rate is the MAPQ-60
  group's; the fitted MAPQ keeps ranking.
- **The T2T veto** removes 41–51 % of every aligner's confident errors at
  N = 99.29 %, as on HG002. CERTA on hs1 as the second opinion again keeps
  70–80 % of BWA-MEM2's gain (`--mapq cal`: 86 → 59, against 50). A second
  GRCh38 aligner in its place removes 0–7 %.

So the lead is real but confined to the stricter operating points, and it
is 1.1–1.7× over the best other aligner, not a large factor everywhere.

## A MAPQ from both mates

Two more opt-in options of paired interval mode (2026-11-05;
[truth_pair.sh](truth_pair.sh), [truth_mapq_pair.py](truth_mapq_pair.py);
reports in [results/truth_2026-11_pair/](results/truth_2026-11_pair/)).

**What `--mapq cal` missed.** At the looser operating points the mates
that decide the ranking are the less confident ones, and `--mapq cal`
ranks those by tier and MAPQ alone. On test B at N = 99.44 %, minibwa's N
most confident mates hold 76,245 correctly placed mates that interval
mode's do not. 45,608 of them are PR mates (a mate in a repeat placed by
the pair certificate), which `--mapq cal` puts at 14 and minibwa at 60. On
test A it is 166,825 of 239,660. What ranks such a mate is the other mate:
minibwa, as BWA-MEM, raises a repeat mate's MAPQ by its unique mate's.

**`--mapq pair`: boosted trees over both mates' records**
([pair_model.h](../include/certa/pair_model.h),
[mapq_model.h](../include/certa/mapq_model.h)). For each mapped mate, the
trees read:
- its own record: tier, MAPQ (`--iv-mapq cert`), the gap to the
  second-best locus (XD, else XO, minus XE), XE, the clusters verified
  (log₂), XB, XQ, NM, soft clips, a PR mate's pair margin (AS plus the
  other mate's AS minus XS, the bound on every other proper pair), the
  proper and mate-unmapped flags, and |TLEN| in classes;
- the other mate's tier, MAPQ, gap, XE, XQ, NM, soft clips and margin.

The trees are fitted by Newton boosting on the log loss of "placed more
than 150 bp from the truth", over HG002 tests A and B (every misplaced
mate, and one in eight of the others with weight 8). MAPQ is
−10 log₁₀ of the predicted probability, rounded, at most 60. The header
holds the score threshold of every MAPQ value, so CERTA and the script
compute MAPQ with the same additions and comparisons, and it is
deterministic. The default is 300 trees of depth 4. Deeper trees (400 of
depth 6) were compared on HG002 only, under a rule fixed beforehand: adopt
them for a model if they lower test B's errors by at least 5 % at two of
its three operating points and raise none by more than 5 %. They did for
the T2T model (63 / 32 / 9 per million against 72 / 33 / 10) and not for
the other (115 / 70 / 18 against 118 / 71 / 18).

**`--t2t`: T2T-CHM13 in the same pass.** With `--t2t hs1.cidx --t2t-chain
hs1ToHg38.over.chain.gz`, every chunk of reads is also mapped on T2T-CHM13
the way it is mapped on GRCh38: passes 1 to L, the pair certificate, tier
R and the pair placement ([map_pipeline.cpp](../src/map_pipeline.cpp)
`t2t_chunk`). Each mate's T2T placement is lifted through the chain
([chain.cpp](../src/chain.cpp), lift_sam.py's rule) and compared with its
GRCh38 placement. The comparison gives one of four classes: lifted within
150 bp (same contig and strand), lifted elsewhere, not lifted, or
unmapped there. With `--mapq pair` a second set of trees
([mapq_model_t2t.h](../include/certa/mapq_model_t2t.h)) also reads that
class, the gap to the second-best locus there, XE minus XE there, the MAPQ
and tier there, and the other mate's class and gap. Only MAPQ changes: on
test A, every other field of every GRCh38 record is the same with and
without `--t2t`. Mapping takes 1.9–2.2× as long (test B: 60 s against
32 s on galaxy's 64 cores), at 46.1 GiB for both indexes.

Checks:
- the T2T MAPQ of every mate equals that of a separate CERTA run on hs1
  (test A), and the C++ lift agrees with the Python one (lift_sam.py, with
  chrY's PARs moved to chrX) on 99.975 % of placed mates; the rest are PAR
  placements;
- CERTA's MAPQ equals the script's replay of the trees on every mapped
  mate of tests A and B and of HG005, with and without `--t2t` (the frozen
  build);
- unit tests of the inputs (worked by hand), the tree walk, the MAPQ
  thresholds and the chain lift; mutants Q1–Q8 and N1–N4 are killed (201
  in all, with P18 and K1–K7 below).

**Held-out half of HG002.** Fitted on the even pairs of tests A and B,
scored on the odd ones (misplaced per million among each run's N most
confident mates with a truth, at the locus; [cv_A.md](results/truth_2026-11_pair/cv_A.md),
[cv_B.md](results/truth_2026-11_pair/cv_B.md)), with each run's cost on
the whole of test B on galaxy (CPU seconds as user + sys, the mean of three
runs for CERTA; peak memory):

| run | B, 99.75 % | B, 99.44 % | B, 95.58 % | A, 94.13 % | A, 83.10 % | CPU (B) | memory (B) |
|---|---|---|---|---|---|---|---|
| interval mode, `--iv-pair all --mapq cal` | 172 | 103 | 34 | 1,069 | 267 | **1,879 s** | 24.5 GiB |
| the same, `--mapq pair` | **118** | **71** | **18** | **405** | **70** | 1,979 s | 24.5 GiB |
| the same, `--mapq pair --t2t` (depth 6) | **63** | **32** | **9** | **111** | **23** | 3,876 s | 42.6 GiB |
| minibwa | 137 | 93 | 48 | 433 | 262 | 3,367 s | **9.7 GiB** |
| BWA-MEM2 | 143 | 119 | 119 | 778 | 778 | 16,111 s | 31.5 GiB |
| strobealign | 169 | 115 | 115 | 1,005 | 1,005 | 3,753 s | 13.5 GiB |
| minimap2 | 168 | 134 | 125 | 844 | 844 | 6,677 s | 12.6 GiB |

CERTA's peak memory includes up to 2 GiB of SAM output queued for a reader
slower than the mapper (here `grep | md5sum`); its index is 20.4 GiB
(GRCh38) and 21.3 GiB (hs1).

The halves share loci, so this checks the fit, not its transfer to
another genome; HG005 does that.

**Held out: HG005.** The models were then fitted on every pair of tests A
and B, frozen ([frozen_source.sha](results/truth_2026-11_pair/frozen_source.sha)),
and run once on HG005 ([holdout_H.md](results/truth_2026-11_pair/holdout_H.md);
the truth and the other runs are those of the section above):

| run | 99.75 % | 99.29 % | 98 % | 95 % | MAPQ ≥ 40 | CPU | memory |
|---|---|---|---|---|---|---|---|
| interval mode, `--iv-pair all --mapq cal` | 319 | 86 | 33 | 22 | 72.8 % | **2,149 s** | 24.6 GiB |
| the same, `--mapq pair` | **186** | **60** | **20** | **9** | 89.3 % | 2,281 s | 24.6 GiB |
| the same, `--mapq pair --t2t` | **93** | **25** | **7** | **3** | 95.6 % | 4,630 s | 45.9 GiB |
| minibwa | 208 | 80 | 37 | 37 | 99.8 % | 3,741 s | **10.1 GiB** |
| BWA-MEM2 | 202 | 108 | 108 | 108 | 99.8 % | 17,300 s | 32.9 GiB |
| strobealign | 323 | 108 | 108 | 108 | 99.6 % | 4,000 s | 24.2 GiB |
| minimap2 | 252 | 152 | 110 | 110 | 99.7 % | 6,477 s | 12.7 GiB |
| *with BWA-MEM2's T2T veto:* minibwa | 167 | 48 | 22 | 22 | 99.7 % | 25,560 s | 38.3 GiB |
| *with BWA-MEM2's T2T veto:* BWA-MEM2 | 160 | 60 | 60 | 60 | 99.7 % | 39,119 s | 38.3 GiB |

- **With one index, CERTA is now the most accurate aligner measured at
  every operating point.** `--mapq pair` misplaces 186, 60, 20 and 9 per
  million against minibwa's 208, 80, 37 and 37: 1.1×, 1.3×, 1.9× and 4×
  fewer. At the loosest point it also beats BWA-MEM2 (202).
- **With T2T-CHM13 in the same pass, the lead is 2–12×:** 93, 25, 7 and 3
  per million. It is also ahead of every other aligner given BWA-MEM2's
  T2T veto, which needs a second aligner and a second index: 1.7× at
  99.75 % (93 against 160), 1.9× at 99.29 % (25 against 48), 3.1× at 98 %
  and 7× at 95 %.
- **The HG002 gains held.** Between HG002's held-out half and HG005 the
  ratios to minibwa are similar (test B at 99.44 %: 71 and 32 against 93;
  HG005 at 99.29 %: 60 and 25 against 80), on another genome sequenced on
  another instrument (NovaSeq 6000 against X).
- **Placements are unchanged;** only the order of confidence changes. The
  fitted MAPQ is an error rate, so fewer mates reach MAPQ 40 (89.3 % and
  95.6 % against minibwa's 99.8 %). A pipeline that filters at a fixed
  MAPQ keeps fewer reads than with the other aligners; matched sensitivity
  compares at the same number of reads.
- **The rows of the other aligners equal the section above**, as does
  `--mapq cal`'s, so the scoring is the one used there.
- **Cost** (HG005, galaxy, CPU seconds as user + sys, peak memory):
  `--mapq pair` 2,281 s at 24.6 GiB, with `--t2t` 4,630 s at 45.9 GiB;
  minibwa 3,741 s at 10.1 GiB, BWA-MEM2 17,300 s at 32.9 GiB. The veto
  adds BWA-MEM2 on hs1 (21,819 s at 38.3 GiB). So `--t2t` costs 1.2×
  minibwa's CPU time and 4.5× its memory. CERTA's mapping stage took
  43.8 s, and 76.2 s with `--t2t` (the wall times here are those of the
  scoring pipeline).

What remains unproven: the model was fitted on HG002 with a truth from
its Q100 assembly, so any bias of that truth's labels (certified unique on
one haplotype) carries into the model; HG005 was labelled the same way.
Neither option has been tested in variant calling.

## Variant calling

The MAPQ options of the last three sections were scored by how they rank
misplaced reads. This section tests them where it matters downstream, in
DeepVariant's calls (2026-11-05; scripts in [vc/](vc/)).

**Setup.**
- **Reads:** Google's 30× NovaSeq 6000 PCR-free 2×151 FASTQs of HG002 and
  HG005 (`gs://brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/`).
  Each genome is streamed once to five mappings at the same time
  ([vc/fanout.py](vc/fanout.py); nothing is stored).
- **Regions:** each mapping keeps the records overlapping chr20 and the
  CMRG genes, padded by 1 kb (81.3 Mbp; [vc/vc_map.sh](vc/vc_map.sh)).
- **Aligners:**
  - CERTA is PR #19's frozen build (`953e1f819c018e62`): interval mode with
    `--iv-mapq cert --iv-pair all`, plus `--mapq cal`, `--mapq pair`, or
    `--mapq pair --t2t`;
  - minibwa and BWA-MEM2 run as in the sections above.
- **Calling:** DeepVariant 1.10.0 (WGS model) over the same regions
  ([vc/call.sh](vc/call.sh)).
- **Scoring:** hap.py 0.3.15 with vcfeval, on PASS calls:
  - HG002 against GIAB v5.0q (within the regions) and CMRG v1.00;
  - HG005 against v4.2.1.

  Errors are false negatives plus false positives.

**HG002 against v5.0q.** The mapping cost is that of the 30× run on
galaxy (CPU seconds as user + sys; peak memory). The five ran together, so
their wall times, about 2.7 h, are the download's.

| aligner | SNV errors | SNV recall | SNV precision | indel errors | indel recall | indel precision | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA `--mapq cal` | 6,052 | 94.22 % | 99.32 % | 6,127 | 80.07 % | 92.85 % | **43,555 s** | 21.8 GiB |
| CERTA `--mapq pair` | 4,040 | 96.24 % | 99.45 % | 4,253 | 86.37 % | 95.11 % | 45,126 s | 21.9 GiB |
| CERTA `--mapq pair --t2t` | 3,214 | 97.07 % | 99.51 % | 3,280 | 89.58 % | 96.25 % | 90,488 s | 43.0 GiB |
| minibwa | **1,745** | 98.48 % | **99.67 %** | **1,775** | **94.47 %** | **97.95 %** | 75,512 s | **9.5 GiB** |
| BWA-MEM2 | 1,761 | **98.49 %** | 99.63 % | 1,848 | 94.26 % | 97.84 % | 333,333 s | 34.6 GiB |

Against CMRG v1.00 the order is the same:
- **SNV errors:** 1,719, 1,286 and 1,077 for the three CERTA runs, against
  minibwa's 800 and BWA-MEM2's 805;
- **indel errors:** 868, 620 and 500, against 278 and 291.

**HG005 against v4.2.1** (a genome the models were not fitted on):

| aligner | SNV errors | SNV recall | SNV precision | indel errors | indel recall | indel precision | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA `--mapq cal` | 3,099 | 96.51 % | 99.72 % | 862 | 92.70 % | 98.78 % | **44,470 s** | 21.8 GiB |
| CERTA `--mapq pair` | 1,654 | 98.13 % | 99.87 % | 373 | 96.85 % | 99.48 % | 46,058 s | 21.8 GiB |
| CERTA `--mapq pair --t2t` | 1,253 | 98.65 % | 99.83 % | 251 | 97.95 % | 99.59 % | 93,154 s | 43.2 GiB |
| minibwa | **543** | 99.43 % | **99.91 %** | **63** | **99.46 %** | **99.92 %** | 75,787 s | **9.0 GiB** |
| BWA-MEM2 | 553 | **99.45 %** | 99.88 % | 73 | 99.37 % | 99.91 % | 335,649 s | 34.7 GiB |

v4.2.1 leaves out more of the difficult regions than v5.0q, so every
aligner makes fewer errors. CERTA's fitted MAPQs make 2.3–5.7× minibwa's
SNV errors here, and 4–14× its indel errors.

With the default MAPQ, the same reads mapped again (as for HG002, check 1
below) give 633 SNV and 85 indel errors: 17 % and 35 % more than minibwa.
`--iv-mapq cert` gives 1,133 and 207.

**The MAPQ, not the placements.** Three checks on HG002:

1. *The same reads, mapped again with each MAPQ mode*
   ([vc/vc_subset.sh](vc/vc_subset.sh)): the 9.95 M read pairs in the five
   runs' records, with the placements the same in every mode. This subset
   reproduces the full runs: `--mapq cal` gives 6,047 SNV errors against
   6,052, `--mapq pair` 4,036 against 4,040.

   | CERTA's MAPQ | SNV errors | indel errors | SNV errors outside difficult regions |
   |---|---|---|---|
   | `--iv-mapq point` (BWA-MEM's formula; the default) | **1,956** | **2,060** | **108** |
   | `--iv-mapq cert` (XQ for tiers R, IV and MR) | 3,156 | 3,590 | 340 |
   | `--mapq cal` | 6,047 | 6,128 | 1,912 |
   | `--mapq pair` | 4,036 | 4,254 | 866 |
   | *minibwa (the full run)* | *1,745* | *1,775* | *67* |

2. *MAPQ ≥ 30 set to 60* ([vc/remap60.sh](vc/remap60.sh)). SNV errors fall
   from 6,052 to 3,761 (`--mapq cal`), from 4,040 to 2,689 (`--mapq pair`)
   and from 3,214 to 2,198 (`--mapq pair --t2t`).
3. *The reads at the misses.* The sites are 134 SNVs outside difficult
   regions that CERTA `--mapq cal` misses and minibwa calls: one in ten of
   chr20's 1,334 ([vc/allele_mapq.py](vc/allele_mapq.py)).

   | run | reads showing the variant: mean MAPQ | showing the reference | variant reads below MAPQ 5 |
   |---|---|---|---|
   | CERTA `--mapq cal` | 22.9 | 33.3 | 9.3 % |
   | CERTA `--mapq pair` | 28.7 | 38.2 | 2.7 % |
   | CERTA `--mapq pair --t2t` | 34.8 | 41.4 | 0.7 % |
   | minibwa | 58.3 | 59.5 | 0.3 % |

What the checks show:
- **The fitted MAPQs are biased against the sample's variants.**
  - `--mapq cal`'s table and `--mapq pair`'s trees lower MAPQ with each edit
    a read carries, because misplaced reads carry more edits.
  - A read that shows the sample's own variant carries one edit more than
    its neighbours that show the reference. So at the same locus it gets a
    lower MAPQ: 22.9 against 33.3 with `--mapq cal`.
  - 9.3 % of such reads fall below MAPQ 5, the threshold below which
    DeepVariant drops reads.
  - XQ (`--iv-mapq cert`) is biased the same way, less strongly.
- **The scale matters too.**
  - DeepVariant's model was trained on BWA-MEM, which gives a unique read 60.
  - CERTA's fitted MAPQ is an error rate and tops out at 43–47.
  - Raising MAPQ ≥ 30 to 60 removes a third to a half of the extra errors,
    but not the bias.
- **With BWA-MEM's MAPQ formula, CERTA's placements are close to minibwa's:**
  - 1,956 SNV and 2,060 indel errors, 12 % and 16 % more than minibwa.
  - The gap is in difficult regions, mostly indels: 2,044 indel errors there
    against minibwa's 1,761.
  - Outside difficult regions CERTA has 108 SNV errors against minibwa's 67.
- **The accuracy score of the sections above cannot see this.**
  - It ranks reads by their chance of being misplaced anywhere in the
    genome. It does not check that reads at one locus are ranked alike,
    whatever allele they carry.
  - Those results stand as rankings of misplaced reads.
  - For variant calling, use the default MAPQ.

**Refitting does not remove the bias.** `truth_mapq_pair.py fit --drop`
holds inputs at 0 so the trees cannot split on them; the header stays
compatible, so CERTA runs the refit unchanged
([vc/vc_model.sh](vc/vc_model.sh), [vc/neutral.sh](vc/neutral.sh),
[vc/neutral3.sh](vc/neutral3.sh)). Each refit was fitted on HG002's even
pairs, scored on the odd ones at matched N, and run through DeepVariant
on the HG002 subset:

| `--mapq pair` model | B, 99.75 % | B, 99.44 % | B, 95.58 % | A, 94.13 % | SNV errors | indel errors | variant reads' mean MAPQ (reference reads') |
|---|---|---|---|---|---|---|---|
| all inputs (as committed) | **118** | **71** | **18** | **405** | 4,036 | 4,254 | 28.7 (38.2) |
| N1: no XE, NM, soft clips or XQ, either mate | 130 | 80 | 18 | 466 | 4,128 | 4,814 | 28.7 (38.8) |
| N2: N1, no record MAPQ either | 162 | 97 | 18 | 694 | 4,353 | 5,037 | 27.3 (38.6) |
| N3: none of the record's own inputs (only the gap, ties, the pair flags, \|TLEN\| and the mate) | 194 | 113 | 19 | 659 | 3,901 | 3,921 | 30.5 (38.1) |
| *the default MAPQ* | | | | | *1,956* | *2,060* | *52.1 (57.7)* |
| *minibwa* | 137 | 93 | 48 | 433 | *1,745* | *1,775* | *58.3 (59.5)* |

The bias is in the tiers themselves ([vc/allele_tags.py](vc/allele_tags.py)).
At the 134 sites above:
- **S0:** 61.2 % of the reads that show the reference are S0 (exact);
  none of the reads that show the variant can be.
- **Weaker tiers:** the variant reads fall to PR, IV, R or MR 2.5× as
  often (46.3 % against 16.1 %).

Any MAPQ fitted to rank misplaced reads learns those tiers and their
correlates (the mate's edits, the gap) and ranks variant reads lower.
DeepVariant, trained on BWA-MEM, reads that pattern as a mapping artifact.

**Where the default MAPQ loses: reads CERTA does not place.**
- **How many:** of the 18.07 M reads minibwa places with MAPQ ≥ 20 in the
  regions, CERTA puts 179,002 (0.99 %) more than 1 kb away, on another
  contig, or nowhere ([vc/lost_reads.py](vc/lost_reads.py)).
- **Which reads:** the more a read differs from GRCh38, the likelier this
  is. By minibwa's NM, CERTA loses 970 per million of the reads with NM 0,
  17,289 with NM 3, 60,829 with NM 6, and 317,308 with NM ≥ 10.
- **Where they go:** 90,236 are unmapped; most of the rest are IV records
  elsewhere, often clipped and at low MAPQ.
- **Why:** 10 % are adapter read-through (the reads are untrimmed); 47 %
  are reads minibwa soft-clips (chimeras, SV breakpoints); 43 % are
  full-length but divergent.
- **The cost:** at one missed cluster of four SNVs within 13 bp
  (16:81,086,429–442), the reads that carry it have 7–13 mismatches. CERTA
  places them elsewhere as local hits that clip most of the read, so that
  haplotype loses coverage there.
- **The existing options do not fix it:**
  - `--full-cert` recovers 6.7 % of these reads, at about 4× the time;
  - a 1,200 bp mate window (`--max-dist`; HG002's inserts: p99 954 bp)
    recovers 0.4 % (1,950 SNV errors against 1,956).
- **What could:** 57 % of the lost reads have a mate that CERTA placed
  with MAPQ ≥ 20 within 2 kb. CERTA's mate rescue (`rescue_mate`) cannot
  recover them: it aligns end to end, above L − 54 (about 10 mismatches),
  in a band of 47 gaps.
  - A sample of 40 of them shows why: they are mostly chimeric reads that
    minibwa aligns in part (51M100S, 110S41M, 33M118S), full-length reads
    with 11–21 mismatches, and reads with long deletions
    (46M16D2M1I86M16S, 110M25D3M3D38M).
  - A local alignment in the mate's window, with clipping and long gaps
    allowed, as in BWA-MEM's mate rescue, could place them. On HG002 it
    places most of them, and closes about a third of the SNV gap and half
    of the indel gap; the rest is the MAPQ of reads placed where minibwa
    places them ([below](#a-local-mate-rescue)).

## Smaller indexes

What CERTA's index holds (2026-11-05):

| part | GRCh38 (`-q 22 -s 1`) | T2T-CHM13 (hs1) | per entry |
|---|---|---|---|
| positions | 10.98 GiB | 11.61 GiB | 4 bytes |
| key bits below the bucket | 5.49 GiB | 5.81 GiB | 2 bytes |
| reference sequence | 2.89 GiB | 2.90 GiB | 1 byte per base |
| bucket directory (2²⁸ buckets) | 1.00 GiB | 1.00 GiB | |
| **in memory** | **20.36 GiB** | **21.32 GiB** | |

Beyond the index, `certa map` holds the per-thread workspaces (0.9 MiB
each), the read and record batches in flight (about 0.9 GiB), and up to
2 GiB of SAM queued for a reader slower than the mapper.

**Filtering repeats would save little.** The mapper enumerates a 22-mer's
positions only when they fit a budget (at most 4,096 per strand) or tier
R's anchor cap (4,096); otherwise it reads at most the first 64 (tier SR's
sample). So the positions past the 64th of the 22-mers seen more than
4,096 times are never read. They are 3.2 % of GRCh38's entries (9,797
22-mers, 95 M entries; 5.2 % of hs1's). Dropping them, with each 22-mer's
count kept, would save 0.53 GiB (0.89 GiB on hs1) without changing any
output. It would also take from `certa check` the complete occurrence
lists it verifies against. Not done.

**The key bits are redundant.** An entry's position holds its 22-mer in
the reference, so the 16 key bits below the bucket can be read there
instead of stored. A v4 index stores the high 8 of them, or none:
`certa index --from grch38.cidx -o grch38.k8.cidx --key-bits 8` (or `0`).
The conversion needs no FASTA and no sort (5–11 s on galaxy). The lookup
narrows the bucket by the stored byte, then reads the remaining bits from
the reference at each candidate entry
([core.h](../include/certa/core.h) `lookup_bucket_ref`). The ranges are
those of the full keys, so every output is the same.

Tests:
- every key's range, and those of absent keys, equal the full-key index's
  for 12, 13 and 14 bits below the bucket;
- the mapping and tier R give the same results;
- `check-index` validates v4 files and rejects a wrong stored byte;
- mutants K1–K7 are killed.

Measured on galaxy (64 threads on cores 0–63; `--intervals --iv-mapq
cert --iv-pair all --mapq pair`; CPU seconds as user + sys, the mean of
three runs on test B and one on A; the highest peak memory; the SAM body's
md5):

| index | in memory | B: CPU | B: peak memory | A: CPU | A: peak memory | SAM |
|---|---|---|---|---|---|---|
| v3, 16 key bits | 20.36 GiB | **1,979 s** | 24.5 GiB | **1,211 s** | 21.5 GiB | |
| v4, 8 key bits | 17.62 GiB | 2,387 s (+21 %) | 19.7 GiB | 1,443 s (+19 %) | 18.5 GiB | the same |
| v4, none | **14.87 GiB** | 2,665 s (+35 %) | **15.8 GiB** | 1,572 s (+30 %) | **15.8 GiB** | the same |
| *with `--t2t`, both indexes:* v3 | 41.68 GiB | **3,876 s** | 42.6 GiB | **2,411 s** | 42.6 GiB | |
| *with `--t2t`:* 8 key bits | 36.04 GiB | 4,796 s (+24 %) | 36.9 GiB | 2,965 s (+23 %) | 36.9 GiB | the same |
| *with `--t2t`:* none | **30.38 GiB** | 5,283 s (+36 %) | **31.3 GiB** | 3,190 s (+32 %) | **31.3 GiB** | the same |

- **Memory falls by the index's saving:** with `--t2t`, by 5.65 and
  11.30 GiB, exactly the stored bytes dropped. GRCh38 alone falls by more
  (4.8 and 8.6 GiB): the faster run there outpaces the `grep | md5sum`
  reader, so up to 2 GiB of SAM waits in the output queue.
- **The cost is CPU time:** each lookup waits on one more random read (the
  reference at a candidate's position), +19–36 %. The wall time on test B
  is the same with 8 key bits (41.1 s against 41.2 s) because the reader
  limits it there; with `--t2t` it is 66 s, 79 s and 87 s.
- **A prefetch of those reference bases does not help.** The test added a
  third stage to the mapper's prefetch pipeline: the bytes a lookup will
  read at its candidates, fetched one read ahead. It saved under 1 % (test
  B, 8 key bits: 2,367 s against 2,387 s; none: 2,665 s, unchanged) and
  was not kept. So the cost is not the latency of that read.

## A pairing fix

CI's seed 12 on PR #19 found a pair that `place_pair` reported as proper
although its mates were not. The bug predates PR #19. When the best proper
pair put a mate on a locus that overlapped the mate's single-end
placement, `place_pair` treated the single-end record as that placement:
it checked only that the two overlapped and that the scores agreed. So it
kept the single-end record under the pair's flags. Now the single-end
record is kept only when it is the same alignment (contig, strand,
position and CIGAR; [interval.h](../include/certa/interval.h)
`same_alignment`); otherwise the mate is placed by the pair (tier MR).

A deterministic test reproduces the case in a tandem array (period 37), and
mutant P18, the old rule, is killed on every seed.

On real reads, with `--mapq pair`, the fix changes the placement or mate
fields of 2,128 records of test A (106 per million) and 3,620 of test B
(91 per million). It changes no FLAG; the moved mates' MAPQ changes.
Against the truth it is neutral:
- **Test A:** 597 mates move, 83 with a truth. 46 of those are at the
  locus before and after: 3 move to it and 3 away from it.
- **Test B:** 775 mates move, 49 with a truth. 47 are at the locus before
  and after.
- **Confident errors at matched N:** test A 393.3 and 73.3 per million
  before, 393.0 and 73.3 after; test B 110.6, 66.8 and 17.5 before, 110.8,
  66.9 and 17.5 after (the model fitted on all of tests A and B).

So the fix is one of consistency: a pair flagged proper is now one.

## A local mate rescue

`--iv-rescue local` (opt-in; 2026-11-05) places the reads of
[Where the default MAPQ loses](#variant-calling) with a local alignment
in the confident mate's window, as BWA-MEM's mate rescue does
([interval.h](../include/certa/interval.h) `rescue_mate_local`,
[pair.h](../include/certa/pair.h) `window_local`):
- **When:** the window searches find no proper pair, and a mate without a
  certificate (unmapped, tier R or IV) has a partner whose single-end
  placement has MAPQ ≥ 40 by `--iv-mapq cert`'s rule (its XQ if it has no
  certificate). The read's score must be at least 30 (BWA-MEM's -T) and
  beat what the pair needs (its single-end score − 17).
- **How:** the full DP of the proper-pair window (850 bp of starts, within
  the partner's contig), with CERTA's objective (BWA-MEM's scores, 5 per
  clipped end) and no band, so with clips at either end and gaps of any
  length. It is exhaustive in the window. Its fill runs 8 cells at a time
  in 16-bit lanes and equals the scalar fill cell for cell (tested). The
  rescued mate is written as MR (XQ 0, [XL, XU] the single-end one); its
  MAPQ is the pair's, at most 40.
- **What it never does:** move a certified mate, or replace a pair found
  end to end. The window search proves there is nothing above its
  threshold, so the DP runs only below it.

**Two rules that did not hold.** The first version ran the local alignment
whenever a window search failed, next to any partner with a point MAPQ
≥ 20.
- It placed more reads (lost reads 96,538 against the final 127,511) and
  gave about the same calls (1,872 SNV and 1,892 indel errors against
  1,880 and 1,906).
- But it raised confident errors at matched N on test B (`point`: 179.7 →
  180.5 per million at 99.75 %). Two causes:
  - **Competing with end-to-end pairs:** its pair could outbid one found
    end to end that had moved a certified mate. That mate then stayed at
    its single-end locus with MAPQ 60, wrong in about 14 of 36 cases on
    test A.
  - **The partner's point MAPQ:** a partner with a point MAPQ below 60 was
    usually misplaced, and its mate followed it. With the first rule fixed
    (no competing), 29 % (test A) and 44 % (test B) of the rescued mates
    with a truth were still misplaced. With the partner's MAPQ by
    `--iv-mapq cert`'s rule ≥ 40, 7 % and 20 % are.
- Gating on that MAPQ instead (20, 40 or 60, replayed on the placements)
  gives no rise at any N; 40 is the confidence threshold of the tables
  above.

**Placement against the truth** (HG002 tests A and B on galaxy, 64
threads; [truth_rescue.sh](truth_rescue.sh),
[truth_rescue.py](truth_rescue.py)). Misplaced mates per million among each
run's N most confident, at truth_pair.sh's shares:

| test, MAPQ | `--iv-rescue` | N₁ | N₂ | N₃ | mates moved with a truth: wrong → right, right → wrong | CPU | memory |
|---|---|---|---|---|---|---|---|
| A, `point` | e2e (default) | 1,320.8 | 796.0 | | | 1,171 s | 24.82 GiB |
| A, `point` | local | **1,318.1** | **793.7** | | 806, 12 | 1,177 s | 24.89 GiB |
| A, `cert` | e2e | 1,151.0 | 702.9 | | | 1,169 s | 24.82 GiB |
| A, `cert` | local | 1,151.0 | 702.9 | | 806, 12 | 1,174 s | 24.89 GiB |
| B, `point` | e2e | 179.7 | 143.0 | 88.3 | | 1,882 s | 24.74 GiB |
| B, `point` | local | **179.6** | **142.9** | **88.0** | 114, 3 | 1,952 s | 24.75 GiB |
| B, `cert` | e2e | 219.3 | 117.6 | 73.9 | | 1,879 s | 24.72 GiB |
| B, `cert` | local | 219.3 | 117.6 | 73.9 | 114, 3 | 1,951 s | 24.77 GiB |

N₁, N₂, N₃: 94.13 % and 83.10 % of the 18.82 M mates with a truth (test
A); 99.75 %, 99.44 % and 95.58 % of the 31.93 M (test B).
- **No rise:** with `point` the errors fall at every N; with `cert` the
  rescued mates have MAPQ 0, so nothing changes.
- **The moved mates:** test A moves 871 mates with a truth. 806 go to their
  locus and 12 leave it; 94 of them were unmapped.

**Variant calling** (the subsets of
[Variant calling](#variant-calling), DeepVariant 1.10.0 and hap.py on a2;
mapping cost of the subset, 64 threads;
[vc/mate_rescue.sh](vc/mate_rescue.sh)). Errors are FN + FP:

| HG002 (9.95 M pairs) | v5.0q SNV | v5.0q indel | CMRG SNV | CMRG indel | reads lost to minibwa | CPU | memory |
|---|---|---|---|---|---|---|---|
| PR #19's frozen build, default MAPQ | 1,956 | 2,060 | 861 | 302 | 179,002 | 1,430 s | 24.82 GiB |
| this build, `--iv-rescue e2e` (default) | 1,949 | 2,063 | 858 | 305 | 179,002 | 1,440 s | 24.85 GiB |
| `--iv-rescue local` | **1,880** | **1,906** | 866 | 308 | 127,511 | 1,516 s | 24.97 GiB |
| *its first version (above)* | *1,872* | *1,892* | *863* | *304* | *96,538* | *1,605 s* | *24.97 GiB* |
| *diagnostic: the lost reads given minibwa's records* | *1,849* | *1,889* | *854* | *305* | *0* | | |
| *diagnostic: also minibwa's MAPQ where CERTA places a read at its locus* | *1,769* | *1,785* | *810* | *281* | *0* | | |
| minibwa (the 30× run) | 1,745 | 1,775 | 800 | 278 | | 75,512 s (30×) | 9.5 GiB |

| HG005 (v4.2.1) | SNV | indel | reads lost to minibwa | CPU | memory |
|---|---|---|---|---|---|
| PR #19's frozen build, default MAPQ | 633 | 85 | | 1,354 s | 24.93 GiB |
| this build, `--iv-rescue e2e` (default) | 633 | 85 | 156,602 | 1,361 s | 24.92 GiB |
| `--iv-rescue local` | 633 | 81 | 107,266 | 1,426 s | 25.02 GiB |
| *its first version* | *628* | *81* | *79,960* | *1,520 s* | *24.97 GiB* |
| minibwa (the 30× run) | 543 | 63 | | 75,787 s (30×) | 9.0 GiB |

- **The gain:** on HG002 v5.0q, 69 fewer SNV and 157 fewer indel errors.
  That is a third of the SNV gap to minibwa and over half of the indel
  gap. CMRG does not change (866 and 308 against 858 and 305). On HG005,
  4 fewer indel errors and the same SNV errors (the first version: 5 and
  4 fewer).
- **The cost:** 5 % more CPU on both subsets, and 0.1 GiB more memory.
- **Not at minibwa's level:** CERTA still makes 135 more SNV and 131 more
  indel errors on HG002, and 90 and 18 more on HG005. The two diagnostic
  rows show where they come from.
- **The ceiling of a mate rescue:** give every read CERTA still loses
  minibwa's record ([vc/swap_lost.py](vc/swap_lost.py)) and CERTA makes
  1,849 SNV and 1,889 indel errors.
  - The lost reads are worth about 100 SNV and 174 indel errors. The
    local rescue gets 69 and 157 of those.
  - The rest of the gap to minibwa (104 SNV, 114 indel errors) is in the
    MAPQ of reads CERTA places where minibwa does. Give the 1.57 M of them
    whose MAPQ differs minibwa's MAPQ as well (`swap_lost.py --mapq`),
    and CERTA makes 1,769 SNV and 1,785 indel errors: 24 and 10 more than
    minibwa (CMRG: 810 and 281 against 800 and 278).
  - So the remaining gap is not placement. It is the default MAPQ of reads
    placed where minibwa places them: BWA-MEM's formula with CERTA's
    second-best locus, lower on the reads that show the sample's variants
    ([above](#variant-calling)).
- **The reads still lost:** 127,511 on HG002. Of those with a mate that
  CERTA places within 2 kb with MAPQ ≥ 20 (by `lost_reads.py`'s count),
  half are recovered (101,669 → 50,910; the first version: → 20,779).
  Of the 50,910 left:
  - 14,552 score below 30 under CERTA's scoring: short chimeric pieces
    such as 31M120S, which BWA-MEM's AS counts without the clip penalty;
  - 4,378 have mates starting more than 850 bp apart;
  - the other 31,980 are left by the rules above: a partner below the MAPQ
    gate, a pair found end to end elsewhere, or a window alignment that
    does not beat the single-end one by the unpaired penalty (not counted
    apart).

## The MAPQ of mates placed by the pair

`--pr-mapq found` (opt-in; 2026-11-07). [A local mate rescue](#a-local-mate-rescue)
left most of the gap to minibwa in the default MAPQ of reads placed where
minibwa places them, mostly properly paired (PR) mates:
- **The default, `bound`:** a PR mate's pair term is BWA-MEM's raw MAPQ of
  the margin over S2, the certified upper bound on every other proper pair
  ([pair.h](../include/certa/pair.h) `certify_pair`). S2 is at least the
  pair's score − 10 whenever the certificate needed other anchors, so a PR
  mate in a repeat gets a low MAPQ even when no other pair was found.
- **`found`:** the margin over the best other proper pair actually found
  (`Result::sub_found`: proper, non-overlapping in a mate, within 10 of the
  pair) or over the unpaired alternative (the mates' best single-end scores
  − 17), whichever is higher, as `mem_sam_pe` takes it. The records, their
  certificates and XQ are unchanged.

| run | HG002 SNV | HG002 indel | HG002 segdups | HG005 SNV | HG005 indel | HG005 segdups | test A, `point` (N₁, N₂) | test B, `point` (N₁, N₂, N₃) | CPU (test B) | memory |
|---|---|---|---|---|---|---|---|---|---|---|
| `bound` (default) | 1,880 | 1,906 | 842 | 633 | 81 | 367 | 1,318.1, 793.7 | 179.6, 142.9, 88.0 | 1,952 s | 24.75 GiB |
| `found` | **1,827** | **1,851** | 833 | 626 | 78 | 391 | 2,246.1, 1,020.5 | 222.5, 177.9, 102.8 | 2,013 s | 24.76 GiB |
| minibwa | 1,745 | 1,775 | 784 | 543 | 63 | 308 | | | | |

DeepVariant's errors (FN + FP) on the subsets of [Variant calling](#variant-calling),
and the misplaced mates per million among each run's N most confident
([truth_prq.sh](truth_prq.sh)); CPU and memory on galaxy, 64 threads. With
`--iv-mapq cert`, test B: 219.3 / 117.6 / 73.9 → 195.2 / 158.8 / 88.7.
- **The calls:** 53 fewer SNV and 55 fewer indel errors on HG002, 7 and 3
  on HG005. SNV misses drop below minibwa's on both (HG002 1,424 against
  1,436, HG005 453 against 472), but the false positives rise (366 → 403,
  124 → 173).
- **The ranking:** the misplaced mates among the confident ones nearly
  double on test A. The raised PR mates are 8.5–13 % misplaced, so the
  bound S2 carries information that the pairs found do not.
- **T2T-CHM13 separates them.** Of the raised PR mates, those whose
  `--t2t` placement lifts to the same locus (XC class 1, T2T MAPQ ≥ 40) are
  1.8–2.2 % misplaced; those it places elsewhere, 65–75 %
  ([truth_t2t_gate.py](truth_t2t_gate.py)). Raising only where T2T agrees,
  and MAPQ 0 where T2T places the mate elsewhere, gives 116.3 / 88.1 / 54.8
  on test B (`bound`: 179.6 / 142.9 / 88.0).

**The T2T-gated MAPQ in variant calling** (2026-11-07). The subsets mapped
with `--t2t --feature-tags` on galaxy ([vc/t2t_feat.sh](vc/t2t_feat.sh)).
[vc/t2t_gate.py](vc/t2t_gate.py) rewrote the MAPQ of the `bound` BAM; the
placement and MAPQ of every record matched the `--t2t` run.

| run | HG002 errors | HG002 segdups | CMRG | HG005 errors | HG005 segdups | CPU (30× HG002) | memory |
|---|---|---|---|---|---|---|---|
| `bound` (default) | 3,786 | 842 | 1,174 | 714 | 367 | 43,555 s | 21.8 GiB |
| `found` | 3,678 | 833 | 1,161 | 704 | 391 | 43,555 s | 21.8 GiB |
| `found` only where T2T agrees | 3,672 | 819 | 1,146 | **677** | 355 | 90,488 s | 43.0 GiB |
| MAPQ 0 where T2T places it elsewhere | 3,767 | 825 | 1,159 | 711 | 361 | 90,488 s | 43.0 GiB |
| both | **3,657** | **804** | 1,132 | 680 | 356 | 90,488 s | 43.0 GiB |
| minibwa | 3,520 | 784 | 1,078 | 606 | 308 | 75,512 s | 9.5 GiB |
| BWA-MEM2 | 3,609 | 789 | 1,096 | 626 | 326 | 333,333 s | 34.6 GiB |

Errors are SNV and indel FN + FP (HG002 v5.0q and CMRG, HG005 v4.2.1); the
cost is that of the 30× HG002 runs without and with the T2T pass (`--mapq
cal`, and `--mapq pair --t2t`) and of minibwa's and BWA-MEM2's.
- The gated rule closes a third of the gap to minibwa on HG005 (held out),
  at twice the CPU and memory. The raise does it; the veto adds nothing.
- **Not adopted:** `bound` stays the default and `found` is opt-in.
- **The rest of the gap is not MAPQ.** At the false positives only CERTA
  calls, minibwa places the reads carrying the false allele at the same
  locus with the same or a higher MAPQ. About half of those calls are
  wrong genotypes at real variants. Outside the difficult regions CERTA
  misses 81 SNVs against minibwa's 48: reads carrying the sample's alleles
  are missing from its pileup ([below](#reads-certa-lost)).

## Reads CERTA lost

Of the reads minibwa places with MAPQ ≥ 20 in HG002's targets, 94,300 were
placed by CERTA more than 1 kb away, on another contig, or not at all (PR
20's build with `--iv-rescue local`; [vc/lost_run.sh](vc/lost_run.sh),
[vc/lost_diag.py](vc/lost_diag.py)). Re-scored under CERTA's objective
against CERTA's own record:
- **73 % (68,713): the search did not reach minibwa's locus.** There
  minibwa's alignment scores above CERTA's XL and below its XU: CERTA's
  interval was honest, and its search missed the alignment. Nearly all of
  them have an exact run of ≥ 22 bases there (a q-mer seed exists), but
  none of tier R's parts occurs exactly there.
- **27 % (25,587): CERTA's placement scores at least as well,** mostly
  clipped or chimeric reads; the pair decides. About 5 k have a mate CERTA
  places within 2 kb of minibwa's locus with MAPQ ≥ 40.
- **None score above XU:** no bound was wrong.

**Why the search missed.** Tier R's parts are disjoint, placed greedily,
and grown until rare; in a repeat they grow long (up to 61 bases).
- Example: chr16:81,086,437, where HG002 has four heterozygous SNVs in
  13 bp inside a repeat. minibwa places 8 reads end to end with MAPQ 60;
  CERTA put them at unrelated loci (four records at chr1:211,861,927) with
  MAPQ 0–24.
- At their locus those reads differ in about 7 places, one in every part,
  yet keep exact runs of 54–66 bases. So no chain reaches the locus, and
  the bound is the span bound, 5 below L per part (XU 121–136).
- The lower bound then comes from sampled copies. At each part's start
  these are the 8 copies whose exact match extends furthest. In a repeat
  other copies match the sample's alleles there longer, so the read's own
  copy is never aligned.

**Fix 1: seed votes** ([repeat.h](../include/certa/repeat.h) `repeat_tier`,
`RepeatParams::vote_occ`; on by default, `--iv-vote-occ 0` turns them off).
- **When:** a read tier R did not prove exact (lo < U_eff).
- **How:** every q-mer of the read with 1–256 hits votes for its diagonal,
  on each strand. Runs of diagonals at most 16 apart, within one band, form
  groups. The 8 groups with the most votes per strand are aligned in a band
  of ±16 (`band_score`); a better score raises lo, and the band is kept for
  A*.
- **Why it is sound:** only lo rises, by alignments actually computed. hi =
  max(lo, U_eff) and the chain evaluation are unchanged, and the parts
  (the witness) are the same. A read becomes exact only when lo reaches
  U_eff, which every chain above U_eff was already evaluated against.

**Fix 2: alternatives for the pair** (`RepeatOut::alt`, [interval.h](../include/certa/interval.h)
`place_pair`).
- Tier R keeps up to 4 bands at loci away from A*, scoring ≥ 30 and more
  than lo − 17 (a lower one can never beat the unpaired alternative).
- The pair stage aligns them again and pairs each with the other mate's
  placement and alternatives, as BWA-MEM pairs every alignment it finds.
  They are added after the local rescue, so they never suppress it.
- They also count as other pairs for the pair's MAPQ, and a mate they move
  is written as MR (XQ 0, MAPQ at most 40).

**The votes' settings, on the 63,093 lost pairs** ([vc/vote_sweep.sh](vc/vote_sweep.sh);
a2, 64 threads):

| `--iv-vote-occ`, `--iv-vote-top` | reads still lost | CPU | memory |
|---|---|---|---|
| 0 (off; PR 20) | 94,300 | 61.6 s | 20.61 GiB |
| 64, 8 | 59,212 | 81.0 s | 20.62 GiB |
| **256, 8 (default)** | 58,051 | 94.8 s | 20.63 GiB |
| 1,024, 8 | 57,286 | 130.3 s | 20.65 GiB |
| 256, 16 | 55,928 | 102.8 s | 20.63 GiB |
| 4,096, 8 | 56,623 | 192.5 s | 20.77 GiB |
| 256, 8 with the alternatives | **51,787** | 95.6 s | 20.63 GiB |

**Placement against the truth** (HG002 tests A and B on galaxy, 64
threads; [truth_fix.sh](truth_fix.sh)). Misplaced mates per million among
each run's N most confident:

| test, MAPQ | run | N₁ | N₂ | N₃ | mates moved with a truth: wrong → right, right → wrong | unplaced with a truth | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| A, `point` | PR 20 | 1,318.1 | 793.7 | | | 554 | 1,177 s | 24.89 GiB |
| A, `point` | votes | 1,266.2 | 754.0 | | 2,578, 535 | 0 | 1,398 s | 24.98 GiB |
| A, `point` | votes + alternatives | **1,265.4** | **753.9** | | **2,969, 635** | 0 | 1,431 s | 24.98 GiB |
| A, `cert` | PR 20 | 1,151.0 | 702.9 | | | 554 | 1,174 s | 24.89 GiB |
| A, `cert` | votes + alternatives | 1,148.7 | 700.8 | | 2,969, 635 | 0 | 1,430 s | 24.98 GiB |
| B, `point` | PR 20 | 179.6 | 142.9 | 88.0 | | 168 | 1,952 s | 24.75 GiB |
| B, `point` | votes | 170.8 | 134.8 | 81.0 | 315, 148 | 0 | 2,293 s | 24.85 GiB |
| B, `point` | votes + alternatives | **170.9** | **134.9** | **81.0** | **426, 185** | 0 | 2,340 s | 24.86 GiB |
| B, `cert` | PR 20 | 219.3 | 117.6 | 73.9 | | 168 | 1,951 s | 24.77 GiB |
| B, `cert` | votes + alternatives | 218.9 | 117.2 | 73.7 | 426, 185 | 0 | 2,340 s | 24.85 GiB |

- **Fewer confident errors at every N:** 4–8 % with `point`, fewer with
  `cert`, where the moved mates mostly have MAPQ 0.
- **More reads at their locus:** net 2,334 mates with a truth on test A and
  241 on test B; the alternatives add 291 and 74 of those. No mate with a
  truth is left unplaced.
- **What moves wrong:** 635 mates on test A go from their locus to an
  alignment elsewhere that scores better. The reference causes those, not
  the search; their MAPQ is mostly low, since the confident errors fall.
- **The cost:** 22 % more CPU on test A and 20 % on test B, and 0.1 GiB.

**Variant calling** (the subsets; the BAMs mapped on galaxy with
[vc/dv_map.sh](vc/dv_map.sh), called on a2 with call.sh; errors FN + FP):

| run | HG002 SNV | HG002 indel | HG002 segdups | CMRG | HG005 SNV | HG005 indel | HG005 segdups | reads still lost (HG002) | CPU (HG002, HG005) | memory |
|---|---|---|---|---|---|---|---|---|---|---|
| PR 20 | 1,880 | 1,906 | 842 | 1,174 | 633 | 81 | 367 | 94,300 | 1,497 s, 1,407 s | 21.36 GiB |
| votes + alternatives | **1,863** | **1,890** | **835** | **1,151** | **624** | **79** | **364** | 51,787 | 1,964 s, 1,837 s | 21.41 GiB |
| minibwa | 1,745 | 1,775 | 784 | 1,078 | 543 | 63 | 308 | | 75,512 s, 75,787 s (30×) | 9.5 GiB |
| BWA-MEM2 | 1,761 | 1,848 | 789 | 1,096 | 553 | 73 | 326 | | 333,333 s, 335,649 s (30×) | 34.6 GiB |

CERTA's cost is that of mapping the subsets on galaxy (64 threads; PR 20's
build mapped again there, [vc/dv_map.sh](vc/dv_map.sh)). Its records on
galaxy (ARM64) and a2 (x86-64) are the same (md5 of the SAM body of HG002's
BAM).
- 33 fewer errors on HG002 and 11 on HG005 (held out); CMRG 23 fewer.
- **The cost:** 31 % more CPU on both subsets, more than on the truth tests
  (20–22 %): the subsets hold every read any of five aligners placed in
  the targets, so their share of repeat reads is higher.
- CERTA is still behind minibwa: by 233 errors on HG002 and 97 on HG005.
- **The reads still lost:** 51,787. 24 k are search misses, 14 k of them
  with ≥ 10 differences at minibwa's locus. 28 k are clipped or chimeric
  reads that score as well where CERTA places them, so the pair decides.
  Truth decides between those, not minibwa: on tests A and B the fixes
  lowered the confident errors, and the reads they moved went to their
  locus 4.7 and 2.3 times as often as away from it.

**Tests.** [test_certa.cpp](../tests/test_certa.cpp):
- `test_tier_r_votes`: 60 copies of a 300-bp element, each with its own
  substitutions at 5 % of sites, and a 6-bp tandem repeat. Reads carry the
  consensus base at 70 % of their copy's own sites, plus 2 substitutions.
  - With votes, every bound stays sound against the local oracle: lo ≤ OPT
    ≤ hi, exact ⇒ lo = OPT, and the CIGAR rescores to lo.
  - lo never falls, hi = max(lo, hi without votes), and the parts are the
    same.
  - lo reaches OPT on 178–180 of 180 reads (139–164 without; seeds
    20260927, 1, 2, 13); the test requires 97 %.
- `test_pair_alts`:
  - A mate whose best alignment (95) is at copy B, and whose alignment
    beside its partner (85, clipped below the window search's threshold)
    is at copy A. It stays at B without the alternatives, and moves to A as
    MR with them (XQ 0, MAPQ 40, XL 95).
  - Two repeat mates, each best elsewhere (95), whose alternatives (90 and
    90) make the only proper pair (180 against the unpaired 173). They pair
    only through their alternatives, both as MR with MAPQ 0.
- Votes are also on in the interval-mode configurations of the tier R,
  pair-placement and `certa check` tests. All checks pass on four seeds.
- Mutants RV1–RV5 (votes, alternatives) and RA1–RA3 (pairing them): all
  killed on seeds 1 and 2.
- Mutants F1–F12 of `--pr-mapq found` ([above](#the-mapq-of-mates-placed-by-the-pair)):
  F1, F3, F5, F6 and F8 survive (a replaced pair counted as found, the
  overlap judgment, the properness checks and the second window search).
  The option is opt-in, and those paths need their own threshold families
  before it could become a default.

## How far from the limit

The question (2026-11-19): has minibwa reached the limit of biologically
correct placement, and what is left for CERTA? Two measurements on HG002's
truth sets ([truth_oracle.py](truth_oracle.py),
[../tools/truth_oracle.cpp](../tools/truth_oracle.cpp),
[truth_personal.py](truth_personal.py); reports in
[results/truth_2026-11_limit/](results/truth_2026-11_limit/)).

**What GRCh38 lets an aligner know.** For every mate a run leaves away from
its truth, the oracle aligns the read at its truth and at the run's
placement (the full DP of a 300-bp window, under CERTA's objective:
BWA-MEM's scores). It then compares the best pair with that mate at its
truth against the pair the run wrote; a pair that is not proper costs
BWA-MEM's unpaired penalty of 17. Each such mate is one of three kinds:
- **missed:** the pair at the truth scores more. An aligner that found
  every alignment would have placed the mate there, so the search failed.
- **tie:** the same score. Neither the read nor GRCh38 tells the loci apart.
- **reference:** the truth scores less. No aligner that maximizes the score
  on GRCh38 places the mate there. Its origin differs from GRCh38 more than
  another locus does: a copy GRCh38 lacks or collapses, or the sample's own
  variants.

Test A (simulated, 18.82 M mates with a truth) and test B (real, 31.93 M),
with each run's misplaced mates per million among its N most confident (its
own MAPQ; in brackets, how many of them are missed or ties) and its cost on
galaxy's 64 cores (CPU seconds as user + sys; peak memory). CERTA is interval
mode with `--iv-pair all --iv-rescue local` and its default (point) MAPQ:

| test A | at the locus | missed | tie | reference | 94.13 % | 83.10 % | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| minibwa | 97.811 % | 9,668 | 139,864 | 262,385 | **429.1** (3.9) | **261.7** (2.0) | 2,167 s | **10.3 GiB** |
| BWA-MEM2 | 97.828 % | **2,797** | 138,979 | 266,949 | 766.6 (3.4) | 766.6 (3.4) | 9,759 s | 33.4 GiB |
| CERTA (the seed votes, [above](#reads-certa-lost)) | 97.805 % | 21,699 | 134,586 | 256,856 | 1,265.4 (55.0) | 753.9 (34.9) | **1,497 s** | 25.1 GiB |

| test B | at the locus | missed | tie | reference | 99.75 % | 99.44 % | 95.58 % | CPU | memory |
|---|---|---|---|---|---|---|---|---|---|
| minibwa | 99.909 % | 382 | 9,435 | 19,155 | **131.7** (0.8) | **88.7** (0.8) | **45.8** (0.1) | 3,367 s | **9.7 GiB** |
| BWA-MEM2 | 99.908 % | **287** | 9,623 | 19,312 | 137.6 (1.0) | 114.3 (0.8) | 114.3 (0.8) | 16,111 s | 31.5 GiB |
| CERTA (the seed votes) | **99.912 %** | 312 | 8,617 | 19,151 | 170.9 (3.0) | 134.9 (2.7) | 81.0 (1.5) | **2,433 s** | 25.1 GiB |

The other aligners' runs are those of the tables above (October 31); the
reference class is split by how much less the truth scores in the reports.

**The sample's own genome.** The reads of test A mapped by minibwa to each
of HG002's haplotype genomes; for each pair, the haplotype where it scores
more, projected to GRCh38 as the truth is:
- 98.608 % of the mates land at their truth, against 97.81–97.83 % for the
  GRCh38 runs. 0.95 % land elsewhere (repeats within the sample's own
  genome) and 0.44 % do not project.
- The mates of pairs that both haplotypes place alike (or that one
  haplotype scores 10 more), with MAPQ ≥ 40 there, are 96.43 % of all
  mates; 3.0 per million of them are misplaced. On GRCh38, minibwa's
  96.51 % most confident misplace 1,702 per million, and its 94.13 % 429.
- Of minibwa's reference-class mates, 46 % land at their truth on the
  sample's genome (13 % confidently); of its ties, 51 % (10 %).

(The haplotype mappings were not timed: they are a diagnostic that needs the
sample's own assembly.)

What this says:
- **On GRCh38, placement is at its limit.** minibwa misses 0.05 % of
  simulated mates and 0.001 % of real ones. Everything else it misplaces,
  no aligner that maximizes this score on GRCh38 can place right. Among its
  confident mates the share is 99 %: 425.1 of its 429.1 misplaced per
  million on test A, 88.0 of 88.7 on test B are reference-caused.
- **CERTA's search was short of that limit**: on simulated reads it missed
  2.2× minibwa's mates, nearly all pairs whose two mates lie in one repeat.
  The joint seeds of the [next section](#joint-seeds-of-a-pair) close that
  gap.
- **What separates the aligners' confident errors is MAPQ**: which of the
  reference-class mates each calls confident. There the limit is far away.
  On the sample's own genome, 3 per million at 96 % of the mates, two
  orders of magnitude below every GRCh38 aligner. CERTA's fitted MAPQ with
  T2T-CHM13 recovers part of it (2–12× fewer than minibwa, held out on
  HG005: [A MAPQ from both mates](#a-mapq-from-both-mates)). The rest
  needs sequence that neither GRCh38 nor CHM13 holds: a pangenome, or the
  sample's own assembly.
- Ties are 0.7 % of simulated mates and 0.03 % of real ones. Every aligner
  gives most of them MAPQ 0, which is the right answer.

## The fitted MAPQ on BWA's scale

`--mapq-scale bwa` (2026-11-19; the default with `--mapq pair`;
[pair_model.h](../include/certa/pair_model.h) `kBwaScale`).

**Why.** [Variant calling](#variant-calling) found that DeepVariant makes
1.8–5.7× minibwa's SNV errors on records with the fitted MAPQs, and named
two causes: the inputs that count a read's own edits, and the scale. A refit
without those inputs did not help. The scale alone was never tested.

The model's MAPQ m is a Phred: 10^(−m/10) of such records are misplaced.
BWA-MEM's MAPQ uses the same numbers for something else. On HG002's truth
sets (tests A and B, [How far from the limit](#how-far-from-the-limit)), the
mates minibwa and BWA-MEM2 give 60 are misplaced at 0.5–8 in 10,000 (Phred
31–43), those they give 50–54 at 1–5 % (Phred 13–20), and those they give
20–29 at 27–41 % (Phred 4–6). A caller trained on BWA-MEM reads a MAPQ of
35 as a doubtful alignment; to the model it means 3 errors in 10,000.

**The scale.** `kBwaScale` maps the model's Phred to the MAPQ minibwa gives
records misplaced as often: 60 from Phred 40, four Phred per level down to
55 at Phred 20, then minibwa's measured levels (Phred 15 → 51, 10 → 45,
5 → 22). It is monotone, so the order of the records is the model's, with
ties where levels merge; `XP` keeps the Phred. It was written into the code
from minibwa's levels before any calls; two tables that reach 60 sooner
(from Phred 25 and from Phred 33) were tried beside it.

**DeepVariant** (1.10.0; hap.py vcfeval; errors = FN + FP; PR #19's
placements, the MAPQ rewritten by [vc/bwscale.sh](vc/bwscale.sh)). The
mapping cost is that of the 30× runs of [Variant calling](#variant-calling)
(galaxy, CPU seconds as user + sys; peak memory); the scale costs nothing:

| HG002, the subset of 9.95 M pairs | v5.0q SNV | v5.0q indel | CMRG SNV | CMRG indel | CPU (30×) | memory |
|---|---|---|---|---|---|---|
| default MAPQ (`point`) | 1,956 | 2,060 | 861 | 302 | 43,555 s | 21.8 GiB |
| `--mapq pair`, Phred | 4,036 | 4,254 | 1,284 | 620 | 45,126 s | 21.9 GiB |
| `--mapq pair`, 60 from Phred 25 | **1,905** | 1,956 | 856 | 301 | 45,126 s | 21.9 GiB |
| `--mapq pair`, 60 from Phred 33 | 1,910 | **1,954** | 860 | **300** | 45,126 s | 21.9 GiB |
| `--mapq pair`, `kBwaScale` (60 from Phred 40) | 1,911 | 1,964 | **850** | 302 | 45,126 s | 21.9 GiB |
| *minibwa* | *1,745* | *1,775* | *800* | *278* | *75,512 s* | *9.5 GiB* |

| the 30× runs | HG002 SNV | HG002 indel | CMRG SNV | CMRG indel | HG005 SNV | HG005 indel | CPU (HG002) | memory |
|---|---|---|---|---|---|---|---|---|
| `--mapq pair`, Phred | 4,040 | 4,253 | 1,286 | 620 | 1,654 | 373 | 45,126 s | 21.9 GiB |
| `--mapq pair`, 60 from Phred 25 | | | | | 609 | 73 | 45,126 s | 21.9 GiB |
| `--mapq pair`, `kBwaScale` | | | | | 600 | 73 | 45,126 s | 21.9 GiB |
| `--mapq pair --t2t`, Phred | 3,214 | 3,280 | 1,077 | 500 | 1,253 | 251 | 90,488 s | 43.0 GiB |
| `--mapq pair --t2t`, 60 from Phred 25 | **1,735** | 1,946 | 801 | 295 | 574 | 76 | 90,488 s | 43.0 GiB |
| `--mapq pair --t2t`, `kBwaScale` | 1,748 | 1,947 | **791** | 293 | 555 | 74 | 90,488 s | 43.0 GiB |
| minibwa | 1,745 | **1,775** | 800 | **278** | **543** | **63** | 75,512 s | **9.5 GiB** |
| BWA-MEM2 | 1,761 | 1,848 | 805 | 291 | 553 | 73 | 333,333 s | 34.6 GiB |

HG005's default MAPQ (the subset): 633 SNV and 85 indel errors. HG005 was
not used to choose anything here.

- **The scale was the cause.** On the same records, the BWA-MEM scale
  halves DeepVariant's errors (HG002: 4,036 → 1,911 SNV; HG005: 1,654 →
  600), and the fitted MAPQ then does better than the default MAPQ (HG002:
  1,911 against 1,956 SNV and 1,964 against 2,060 indel errors; HG005: 600
  against 633 and 73 against 85).
- **With `--t2t`** CERTA's SNV calls are level with minibwa's on HG002
  (1,748 against 1,745) and ahead on CMRG (791 against 800); on HG005, 555
  against minibwa's 543 and BWA-MEM2's 553.
- **Indels stay behind minibwa:** 1,947 against 1,775 on HG002, 74 against
  63 on HG005.
- **The tables:** on HG002 the three give the same calls within 0.5 %
  (all errors with `--t2t`: 4,777 with 60 from Phred 25, 4,779 with
  `kBwaScale`). `kBwaScale` keeps the model's order among more of the
  confident records, and on HG005 it is the better of the two tried there.

**What the scale costs the ranking** (the default build's `--mapq pair`
runs, [The default build](#the-default-build-all-together); misplaced per
million among the N most confident, ranked by `XP` or by the MAPQ written;
[scale_cost.md](results/truth_2026-11_limit/scale_cost.md)):

| HG005 (held out) | 99.75 % | 99.29 % | 98 % | 95 % | MAPQ 60 |
|---|---|---|---|---|---|
| `--mapq pair`, by `XP` | 207.3 | 61.7 | 20.7 | 9.4 | |
| `--mapq pair`, by MAPQ | 207.3 | 61.7 | 21.2 | 9.4 | 89.3 % |
| `--mapq pair --t2t`, by `XP` | 99.2 | 25.2 | 7.2 | 2.6 | |
| `--mapq pair --t2t`, by MAPQ | 99.2 | 26.1 | 7.4 | 3.0 | 95.5 % |
| *minibwa* | *208.0* | *79.9* | *37.4* | *37.4* | |

The scale merges the levels above Phred 40 into 60, and below 20 maps
several Phreds to one level; that moves a few mates across the boundaries
of N: at most 0.9 per million here, 1.3 on test B and 2.9 on test A.

## Joint seeds of a pair

`--iv-pair-votes 8` (on by default; `0` turns it off), with
`--iv-pair-vote-occ 4096` and `--iv-pair-vote-budget 8192` (2026-11-19;
[interval.h](../include/certa/interval.h) `pair_votes`, `place_pair`).

**Why.** The missed mates of [How far from the limit](#how-far-from-the-limit)
are mostly pairs whose two mates lie in one repeat: 87 % of CERTA's on test
A have MAPQ 0, and in a sample of 4,000 such pairs about four in five had
neither mate certified unique (both in tiers SR, R or IV, or moved by the
pair, MR). Each mate's best single-end alignment lies at another
copy, for instance 150 for one mate at copy A and 150 for the other at copy
B, while the copy they come from has each one mismatch down (145 + 145 =
290 > 150 + 150 − 17). The pair stage searched beside each mate's own best
placement and its proven copy, where the other mate scores 130; tier R's
alternatives (at most 4 per mate) rarely hold the right copy when many
copies score near the best. BWA-MEM pairs every hit of both mates; here
the mates' seeds find the copies that fit both.

**How.** When no pair found so far scores the sum of the mates' best
single-end scores (so a better pair may exist):
- the other mate is searched in the proper-pair window of each of a mate's
  alternatives, as in the window of its single-end placement;
- every 22-mer of each mate with 1–`occ` hits votes for its diagonal,
  rarest first, until `budget` diagonals per mate and strand. Rare q-mers
  are the ones that tell copies apart: they cover the sites where copies
  differ. Runs of diagonals within ±16 form groups;
- a joint locus is a group of one mate on the forward strand and a group of
  the other on the reverse strand that fit a proper pair. The 8 joint loci
  with the most votes are aligned, each mate in its group's band, and every
  pair above max(the unpaired alternative, the best pair found − 10) is a
  candidate. `place_pair` chooses among all candidates as before.
- Joint loci with the same votes are ordered by a hash of the pair's bases
  and the loci, so the pairs from the copies of a repeat spread over those
  copies instead of all taking the first. The order depends on the reads
  only.
- A mate certified unique (S0, S1, SL) with a point MAPQ ≥ 20 anchors the
  pair: the joint seeds are skipped, and that mate is not searched for
  beside the other mate's alternatives. The window searches around it are
  the search that matters there.

Nothing it adds is a certificate: a mate it moves is tier MR with XQ 0, as
before. With it off (`--iv-pair-votes 0`) the SAM output is byte-identical
to the build before it (test B, md5).

**Against the truth** (HG002 tests A and B on galaxy, 64 threads; the
oracle's classes as in [How far from the limit](#how-far-from-the-limit);
misplaced per million among each run's N most confident with the default
MAPQ; CPU seconds as user + sys, peak memory; reports
[oracle_A.md](results/truth_2026-11_limit/oracle_A.md),
[oracle_B.md](results/truth_2026-11_limit/oracle_B.md),
[final_A.md](results/truth_2026-11_limit/final_A.md),
[final_B.md](results/truth_2026-11_limit/final_B.md)):

| test A | at the locus | missed | moved wrong → right, right → wrong | 94.13 % | 83.10 % | CPU | memory |
|---|---|---|---|---|---|---|---|
| CERTA, joint seeds off | 97.805 % | 21,699 | | 1,265.4 | 753.9 | **1,497 s** | 25.1 GiB |
| occ 512, no budget | 97.834 % | 9,754 | 6,300, 729 | 1,253.6 | 748.4 | 2,103 s | 25.1 GiB |
| occ 2048, no budget | **97.847 %** | **4,242** | 8,906, 884 | **1,252.9** | **748.2** | 2,580 s | 24.6 GiB |
| occ 4096, budget 8192 (the default) | 97.846 % | 6,165 | 8,556, 779 | 1,262.3 | 752.4 | 2,205 s | 24.8 GiB |
| *minibwa* | *97.811 %* | *9,668* | | *429.1* | *261.7* | *2,167 s* | *10.3 GiB* |
| *BWA-MEM2* | *97.828 %* | *2,797* | | *766.6* | *766.6* | *9,759 s* | *33.4 GiB* |

| test B | at the locus | missed | moved wrong → right, right → wrong | 99.75 % | 99.44 % | 95.58 % | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA, joint seeds off | 99.912 % | 312 | | 170.9 | 134.9 | 81.0 | **2,433 s** | 25.1 GiB |
| occ 512, no budget | 99.912 % | 143 | 243, 103 | 170.7 | 134.6 | 80.9 | 3,328 s | 25.3 GiB |
| occ 2048, no budget | 99.913 % | **142** | 256, 112 | **170.6** | **134.5** | 80.9 | 4,249 s | 25.5 GiB |
| occ 4096, budget 8192 (the default) | **99.913 %** | 147 | 255, 100 | 170.8 | 134.7 | 80.9 | 3,521 s | 25.2 GiB |
| *minibwa* | *99.909 %* | *382* | | *131.7* | *88.7* | *45.8* | *3,367 s* | *9.7 GiB* |
| *BWA-MEM2* | *99.908 %* | *287* | | *137.6* | *114.3* | *114.3* | *16,111 s* | *31.5 GiB* |

- **CERTA now places more mates at their locus than minibwa and BWA-MEM2
  on both tests**, and its search misses fewer mates than minibwa on both
  (6,165 against 9,668 simulated; 147 against 382 real). BWA-MEM2 still
  misses fewer simulated mates (2,797).
- **Moved mates go to their locus far more often than away from it**: 11
  times as often on test A (8,556 against 779), 2.6 times on test B (255
  against 100). The joint seeds only add candidates, and a tie keeps the
  pair found first, so a mate that leaves its locus goes to a pair that
  scores more: GRCh38's reference class.
- **The confident errors barely change** (test A 1,265 → 1,262 per
  million): most moved mates have MAPQ 0, being repeats. The gain is in
  where the reads are, which coverage and paralog-aware callers use.
- **The cost:** 47 % more CPU on test A and 45 % on test B, which brings
  CERTA's interval mode to minibwa's CPU time (test B: 3,521 s against
  3,367 s); memory is unchanged. Where it goes, at occ 512 on test A (609 s
  added): the window searches beside the alternatives 307 s and the votes
  304 s (each removed in turn from the build). With the budget, occ 4,096
  costs about what occ 512 does (2,201 s against 2,107 s); occ 2,048
  without it costs 2,580 s. Skipping anchored mates in the window searches
  saved nothing measurable.
- **The settings:** occ 4,096 with the budget was chosen on the 4,000
  missed pairs of test A (3,146 of 7,967 mates placed at their locus,
  against 2,376 with occ 512) and on its cost there, before the default's
  runs above. `--iv-pair-votes 0` gives back the time.

**Tests** ([test_certa.cpp](../tests/test_certa.cpp) `test_pair_votes`):
- three copies of a 500-bp fragment: copy A holds one mate exactly and the
  other with 4 substitutions, copy B the reverse, copy T each with one.
  Without the joint seeds the mates stay unpaired at A and B; with them
  (top 8 or top 1) they pair at T, 7 above the unpaired alternative, both
  MR with MAPQ 0 and XQ 0, and the CIGARs rescore;
- with only unique q-mers allowed to vote nothing pairs, unless a mate's
  alternative at T is given: the other mate is then found in its window;
- a pair already at its mates' bests is the same with and without them;
- with two copies of T, 16 shifted pairs split between them (10 and 6 on
  seed 20260927), and each takes the same copy again.

All checks pass on seeds 20260927, 1 and 2. Mutants PV1–PV6 (the joint
seeds, the window searches beside the alternatives, the order of the joint
loci, the pruned scan), BS1–BS2 (the scale) and, again on this build,
RV1–RV5 and RA1–RA3 (tier R's votes and alternatives): all 16 killed on
seeds 1 and 2 ([tests/mutants.tsv](../tests/mutants.tsv)). CERTA's
`--mapq pair` on the BWA scale equals the fit's replay through `kBwaScale`
(`truth_mapq_pair.py replay --bwa`) on every one of test A's 19,940,557
mapped mates.

## The default build, all together

Paired interval mode as it now runs by default (`--intervals --iv-pair all
--iv-rescue local`: the seed votes, tier R's alternatives and the joint
seeds), with its three MAPQs, against minibwa and BWA-MEM2 (2026-11-20;
build `001a8c2f8ec5fbd9`; [ranking_H.md](results/truth_2026-11_limit/ranking_H.md)).

**HG005, held out** (30.13 M mates with a truth): mates away from their
locus, misplaced per million among each run's N most confident (CERTA's
`--mapq pair` ranked by the MAPQ it writes, on BWA-MEM's scale), and the
cost of the 20 M pairs on galaxy's 64 cores:

| run | away from the locus | 99.75 % | 99.29 % | 98 % | 95 % | CPU | memory |
|---|---|---|---|---|---|---|---|
| CERTA, default MAPQ | **28,309** | 317.6 | 124.5 | 71.2 | 71.2 | 3,915 s | 24.8 GiB |
| CERTA, `--mapq pair` | **28,309** | 207.3 | 61.7 | 21.2 | 9.4 | 3,995 s | 24.8 GiB |
| CERTA, `--mapq pair --t2t` | **28,309** | **99.2** | **26.1** | **7.4** | **3.0** | 7,602 s | 46.1 GiB |
| minibwa | 28,344 | 208.0 | 79.9 | 37.4 | 37.4 | **3,741 s** | **10.1 GiB** |
| BWA-MEM2 | 28,410 | 202.2 | 107.9 | 107.9 | 107.9 | 17,300 s | 32.9 GiB |

**DeepVariant** (the subsets of [Variant calling](#variant-calling): 9.95 M
pairs of HG002, the same extraction for HG005; errors = FN + FP; the CPU
and memory columns repeat the HG005 truth runs above, the one cost all five
were measured on with the same reads):

| run | HG002 SNV | HG002 indel | CMRG SNV | CMRG indel | HG005 SNV | HG005 indel | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,870 | 1,894 | 843 | 299 | 637 | 79 | 3,915 s | 24.8 GiB |
| CERTA, `--mapq pair` | 1,841 | 1,826 | 851 | 295 | 592 | 77 | 3,995 s | 24.8 GiB |
| CERTA, `--mapq pair --t2t` | **1,696** | 1,807 | **789** | 280 | 555 | 77 | 7,602 s | 46.1 GiB |
| minibwa | 1,745 | **1,775** | 800 | **278** | **543** | **63** | **3,741 s** | **10.1 GiB** |
| BWA-MEM2 | 1,761 | 1,848 | 805 | 291 | 553 | 73 | 17,300 s | 32.9 GiB |

- **Placement:** CERTA leaves the fewest mates away from their locus on all
  three truth sets (HG005: 28,309 against 28,344 and 28,410; tests A and B
  in [Joint seeds of a pair](#joint-seeds-of-a-pair)). The differences on
  real reads are small, as the limit says they must be.
- **Confidence, held out:** `--mapq pair` misplaces fewer of its confident
  mates than minibwa at 99.29 % and below (1.3–4×) and as many at 99.75 %;
  with `--t2t`, 2–12× fewer at every point.
- **Variant calls:** with `--t2t`, fewer SNV errors than minibwa and
  BWA-MEM2 on HG002 (v5.0q 1,696 against 1,745 and 1,761; CMRG 789 against
  800 and 805; segmental duplications 596 against 627 and 630) and the
  highest SNV recall on both genomes. On HG005, held out, 2 % more SNV
  errors than minibwa (555 against 543; BWA-MEM2 553). Indels: level on
  HG002 (1,807 against 1,775 and 1,848), behind on HG005 (77 against 63).
  Outside the difficult regions CERTA still makes more SNV errors than
  minibwa (HG002 86 against 67; HG005 19 against 15): reads it does not
  place where minibwa does ([Reads CERTA lost](#reads-certa-lost)).
- **Cost:** without `--t2t`, about minibwa's CPU time (+7 % on HG005) at
  2.5× its memory; with `--t2t`, twice its CPU time and 4.6× its memory.
  On the variant-calling subsets, whose reads favour repeats, the joint
  seeds cost more than on the truth sets: HG002 3,291 s against 1,964 s
  without them (+68 %; HG005 +64 %).

## The pair model refit on the default build's placements

`--mapq pair` with both models fitted again, on the records of the default
build (2026-11-20; [truth_pair.sh](truth_pair.sh) `features`, `fit`, `check`
and `holdout`; reports in [results/truth_2026-11_refit/](results/truth_2026-11_refit/)).

**Why.** The trees of [A MAPQ from both mates](#a-mapq-from-both-mates) were
fitted on the records of PR #19's build. Since then the local mate rescue,
tier R's seed votes and alternatives and the joint seeds have moved mates
and changed what their records hold (more MR mates, new margins). On HG005,
`--mapq pair` misplaced 186 per million of its 99.75 % most confident mates
on PR #19's placements and 207.3 on the default build's.

**What.**
- **The records:** tests A, B and H mapped by the default build (the source
  of commit 320a506) with `--intervals --iv-mapq cert --iv-pair all
  --iv-rescue local --t2t --feature-tags`.
- **The trees:** fitted as on 2026-11-05, 300 of depth 4 for the pair model
  and 400 of depth 6 for the T2T model, on every misplaced mate of tests A
  and B and one in eight of the others.
- **One new input tried:** XM, the margin of a pair placed from the
  alignments found (`truth_mapq_pair.py fit --xm`). The joint seeds' pairs
  carry it, and the trees did not read it.
- **Two rules, fixed before the runs they decide:**
  - XM becomes an input of a model if, fitted on HG002's even pairs, it
    lowers test B's misplaced mates per million on the odd pairs by at least
    5 % at two of the three N and raises none by more than 5 % (the rule of
    the depth choice).
  - The refit replaces a model of 2026-11-05 unless, on HG005 (held out, run
    once), it misplaces more than 5 % more mates at two or more of the four N.
- **The fits:** `fit --workers 7` computes each level's histograms in seven
  forked processes. The trees are those of a serial fit, bit for bit
  (checked on synthetic data), and the eight fits took 26 minutes on 64
  cores; serially, the T2T fits on every pair had not reached round 50 of
  400 after 16 minutes.

**Cross-validated on HG002** (fitted on the even pairs, scored on the odd
pairs at matched N; [cv_A.md](results/truth_2026-11_refit/cv_A.md),
[cv_B.md](results/truth_2026-11_refit/cv_B.md)). The cost is that of mapping
test B with each configuration (galaxy, 64 cores; CPU seconds as user + sys,
peak memory); XM costs nothing to compute:

| model | B, 99.75 % | B, 99.44 % | B, 95.58 % | A, 94.13 % | A, 83.10 % | CPU (B) | memory (B) |
|---|---|---|---|---|---|---|---|
| `--mapq pair` | 117.1 | 71.2 | 17.9 | 403.5 | **67.8** | 3,473 s | 24.3 GiB |
| the same, with XM | **116.5** | **71.0** | **17.8** | **383.0** | 70.9 | 3,473 s | 24.3 GiB |
| `--mapq pair --t2t` | **60.8** | 31.3 | **8.7** | 110.2 | 23.4 | 6,350 s | 45.4 GiB |
| the same, with XM | 61.8 | **31.0** | 8.8 | **103.9** | **22.8** | 6,350 s | 45.4 GiB |

- **XM is not an input.** On test B it changes the misplaced mates by 0.3–0.6 %
  for the pair model and 1.0–1.6 % for the T2T model, against the 5 % the
  rule asks. On test A it lowers them by 5–6 % at 94.13 % and raises the pair
  model's by 4.6 % at 83.10 %. With XM the pair model gives MAPQ ≥ 40 to
  fewer mates (test A 56 % against 67 %).
- **The refit holds the level of the 2026-11-05 fit** on its own placements
  (test B 118 / 71 / 18 for the pair model, 63 / 32 / 9 for the T2T model).

**Held out: HG005** (run once with the refit;
[holdout_H.md](results/truth_2026-11_refit/holdout_H.md); misplaced per
million among each run's N most confident of the 30.13 M mates with a
truth; ranked by `XP`, the model's Phred, or by the MAPQ written on
BWA-MEM's scale; the other rows as in [The default build](#the-default-build-all-together);
cost of the 20 M pairs on galaxy's 64 cores):

| run | 99.75 % | 99.29 % | 98 % | 95 % | CPU | memory |
|---|---|---|---|---|---|---|
| CERTA, default MAPQ | 317.6 | 124.5 | 71.2 | 71.2 | 3,915 s | 24.8 GiB |
| CERTA, `--mapq pair`, the trees of 2026-11-05, by `XP` | 207.3 | 61.7 | 20.7 | 9.4 | 3,995 s | 24.8 GiB |
| CERTA, `--mapq pair`, refit, by `XP` | 185.2 | 60.8 | 20.9 | 9.9 | 3,981 s | 24.8 GiB |
| CERTA, `--mapq pair`, refit, by MAPQ | 185.2 | 60.8 | 21.3 | 10.0 | 3,981 s | 24.8 GiB |
| CERTA, `--mapq pair --t2t`, the trees of 2026-11-05, by `XP` | 99.2 | 25.2 | 7.2 | **2.6** | 7,602 s | 46.1 GiB |
| CERTA, `--mapq pair --t2t`, refit, by `XP` | **95.1** | **24.8** | **7.1** | **2.6** | 7,544 s | 46.1 GiB |
| CERTA, `--mapq pair --t2t`, refit, by MAPQ | **95.1** | 25.9 | 7.4 | 2.7 | 7,544 s | 46.1 GiB |
| minibwa | 208.0 | 79.9 | 37.4 | 37.4 | **3,741 s** | **10.1 GiB** |
| BWA-MEM2 | 202.2 | 107.9 | 107.9 | 107.9 | 17,300 s | 32.9 GiB |
| minibwa with BWA-MEM2's T2T veto | 167.4 | 47.5 | 22.3 | 22.3 | 25,560 s | 38.3 GiB |

- **The refit passes its rule.** The pair model misplaces 11 % fewer of its
  confident mates at 99.75 % (185.2 against 207.3 per million) and 1.5 %
  fewer at 99.29 %; it is 1 % worse at 98 % and 5.3 % worse at 95 % (9.9
  against 9.4), one N past 5 %. The T2T model is better or level at every N
  (95.1 against 99.2 at 99.75 %).
- **With one index, `--mapq pair` is again ahead of minibwa at every
  point:** 185.2, 60.8, 21.3 and 10.0 per million against 208.0, 79.9, 37.4
  and 37.4, by the MAPQ written (1.1–3.7× fewer). Before the refit it was
  level with minibwa at 99.75 %.
- **With `--t2t`, 2.2–14× fewer than minibwa**, and fewer than minibwa with
  BWA-MEM2's T2T veto at every point, which needs a second aligner and index.
- **The cost is unchanged:** the trees have the same number and depth.

**DeepVariant** (1.10.0; the variant-calling subsets mapped again by the
refit build, which places every read as the default build does, so only
MAPQ differs; errors = FN + FP; the CPU and memory columns repeat the HG005
runs above):

| run | HG002 SNV | HG002 indel | CMRG SNV | CMRG indel | HG005 SNV | HG005 indel | CPU | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,870 | 1,894 | 843 | 299 | 637 | 79 | 3,915 s | 24.8 GiB |
| CERTA, `--mapq pair`, the trees of 2026-11-05 | 1,841 | 1,826 | 851 | 295 | 592 | 77 | 3,995 s | 24.8 GiB |
| CERTA, `--mapq pair`, refit | 1,830 | 1,831 | 839 | 294 | 597 | 75 | 3,981 s | 24.8 GiB |
| CERTA, `--mapq pair --t2t`, the trees of 2026-11-05 | 1,696 | 1,807 | 789 | 280 | 555 | 77 | 7,602 s | 46.1 GiB |
| CERTA, `--mapq pair --t2t`, refit | **1,671** | 1,806 | **783** | 285 | **534** | 74 | 7,544 s | 46.1 GiB |
| minibwa | 1,745 | **1,775** | 800 | **278** | 543 | **63** | **3,741 s** | **10.1 GiB** |
| BWA-MEM2 | 1,761 | 1,848 | 805 | 291 | 553 | 73 | 17,300 s | 32.9 GiB |

- **With `--t2t`, the refit makes fewer SNV errors than minibwa and
  BWA-MEM2 on both genomes:** HG002 1,671 against 1,745 and 1,761, CMRG 783
  against 800 and 805, and HG005 (held out) 534 against 543 and 553. With
  the trees of 2026-11-05 HG005 was 555.
- **Indels stay behind minibwa:** HG002 1,806 against 1,775, HG005 74
  against 63. They are level with BWA-MEM2 on HG005 (73) and ahead of it on
  HG002 (1,848). Why: the indel bullets of [Next](#next).
- **Without `--t2t` the refit's calls are level with the trees of
  2026-11-05:** HG002 1,830 against 1,841 SNV errors and 1,831 against 1,826
  indel errors, HG005 597 against 592 and 75 against 77. They stay below the
  default MAPQ's (1,870 / 1,894 and 637 / 79).
- The subsets were mapped in 3,271 s (HG002) and 3,082 s (HG005) at
  21.5 GiB, and with `--t2t` in 5,663 s and 5,532 s at 42.9 GiB (galaxy, 64
  cores).

**Checks.**
- CERTA's MAPQ equals the script's replay of the refit trees on every mapped
  mate, with and without `--t2t`: tests A and B (19,940,557 and 39,667,148
  mates) and HG005 (39,568,207 mates), on the build of these headers
  (`truth_pair.sh check` and `holdout`).
- `test_certa` passes with the refit headers on seeds 20260927, 1, 2 and
  11–15 (WSL, g++ 13.3).

## Tier R's wide band

`--iv-wide 1` (on by default; `0` turns it off; 2026-11-20;
[repeat.h](../include/certa/repeat.h) `repeat_tier`, `RepeatParams::wide`).

**Why.** In [hg005_indel_sites.md](results/truth_2026-11_refit/hg005_indel_sites.md),
at 5:1,430,797, a 24-bp insertion in HG005:
- minibwa aligns 24 reads through the insertion (63M24I64M);
- CERTA writes 28 of its 30 reads there clipped (68M83S, AS 63, tier IV)
  and loses the reads with the insertion near their start from the site.

Under CERTA's own scores, the alignment through the insertion scores 97.
Tier R aligns a read in bands that each hold one run of diagonals: a
sampled copy ±8, a chain of rare parts, or a vote group ±16. In a repeat no
part is rare, so no chain spans the indel. The two sides of an indel longer
than 16 bases then fall into different bands.

**How.** Take a read that tier R did not prove exact and whose best
alignment so far scores at least 20 below a perfect match. Tier R also
aligns one band of 95 diagonals centred on that alignment's band, which
holds an indel of up to 47 bases either way. As with the seed votes, only lo
rises, and only by an alignment actually computed. hi, the parts and every
certificate are unchanged.

**Against the truth.** These runs compare the default build with and without
the wide band (`--intervals --iv-pair all --iv-rescue local`, default MAPQ;
reports moved7_A.md, moved7_B.md, moved7_H.md, oracle7_A.md and oracle7_B.md
in [results/truth_2026-11_refit/](results/truth_2026-11_refit/)). The
columns give misplaced mates per million among each run's N most confident.
The cost is that of mapping test B back to back, twice each with SAM to
/dev/null, on galaxy's 64 cores:

| test | run | away from the locus | search misses | N₁ | N₂ | N₃ | N₄ | moved: wrong → right, right → wrong | CPU (B) | memory (B) |
|---|---|---|---|---|---|---|---|---|---|---|
| A | without | 405,364 | **6,165** | **1,262.3** | **752.4** | | | | **3,443 s** | 21.6 GiB |
| A | with | **405,361** | 6,169 | 1,264.8 | 753.5 | | | 20, 17 | 3,453 s | 21.6 GiB |
| B | without | 27,925 | 147 | 170.8 | 134.7 | 80.9 | | | **3,443 s** | 21.6 GiB |
| B | with | **27,922** | **145** | 170.8 | 134.7 | 80.9 | | 4, 1 | 3,453 s | 21.6 GiB |
| H | without | 28,309 | | **317.6** | **124.5** | **71.2** | **71.2** | | **3,443 s** | 21.6 GiB |
| H | with | 28,309 | | 317.8 | 124.7 | 71.3 | 71.3 | 0, 0 | 3,453 s | 21.6 GiB |

N: 94.13 % and 83.10 % (test A); 99.75 %, 99.44 % and 95.58 % (test B);
99.75 %, 99.29 %, 98 % and 95 % (HG005). Away from the locus counts the
misplaced and unplaced mates with a truth.
- **Few reads change:** 4,041 mates on test A, 5,846 on test B and 8,398 on
  HG005, of which 879, 316 and 321 have a truth.
- **The truth sets barely move.** Their labels are certified unique on one
  haplotype, so they hold few reads across a long indel in a repeat. The
  confident errors rise by at most 0.2 %: a misplaced read whose alignment
  improves where it is gains MAPQ.
- **The cost is 0.3 % more CPU** (3,453 s against 3,443 s), with the same
  memory.

**DeepVariant** (1.10.0; the variant-calling subsets mapped by the default
build with the refit models, without and with the wide band; errors are FN
+ FP; [vc_errors_wide.md](results/truth_2026-11_refit/vc_errors_wide.md)):

| run | HG002 SNV | HG002 indel | CMRG SNV | CMRG indel | HG005 SNV | HG005 indel | CPU (HG002) | memory |
|---|---|---|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,870 | 1,894 | 843 | 299 | 637 | 79 | 3,291 s | 21.6 GiB |
| the same, wide band | 1,855 | 1,869 | 836 | 291 | 637 | 78 | 3,266 s | 21.6 GiB |
| CERTA, `--mapq pair` | 1,830 | 1,831 | 839 | 294 | 597 | 75 | 3,271 s | 21.5 GiB |
| the same, wide band | 1,822 | 1,818 | 838 | 292 | 597 | 72 | 3,311 s | 21.6 GiB |
| CERTA, `--mapq pair --t2t` | 1,671 | 1,806 | 783 | 285 | **534** | 74 | 5,663 s | 42.8 GiB |
| the same, wide band | **1,656** | 1,790 | **776** | 279 | **534** | 72 | 5,836 s | 43.0 GiB |
| minibwa | 1,745 | **1,775** | 800 | **278** | 543 | **63** | 75,512 s (30×) | **9.5 GiB** |
| BWA-MEM2 | 1,761 | 1,848 | 805 | 291 | 553 | 73 | 333,333 s (30×) | 34.6 GiB |

CERTA's CPU is that of mapping the HG002 subset on galaxy's 64 cores.
minibwa's and BWA-MEM2's is that of their 30× runs. The wide band's runs
shared the host with another project's jobs; back to back it costs 0.3 %
more (above).
- **HG002:** with every MAPQ, indel errors fall by 13–25 and SNV errors by
  8–15. With `--t2t`, CERTA's total errors are 3,446, against minibwa's 3,520
  and BWA-MEM2's 3,609, and its indel gap to minibwa halves (1,790 against
  1,775; it was 1,806).
- **HG005, held out:** indel errors fall from 75 to 72 with `--mapq pair` and
  from 74 to 72 with `--t2t`, below BWA-MEM2's 73; SNV errors are unchanged.
  Of the four HG005 indel errors traced to the search, three are gone
  (5:1,430,797, 20:38,081,181 and 20:50,366,199). At 5:1,430,797 CERTA now
  aligns 21 of its 38 reads through the 24-bp insertion (before: none of
  30; minibwa: 24 of 39). With `--t2t`, CERTA's total errors on HG005 equal
  minibwa's (606; BWA-MEM2 626).
- **Indels stay behind minibwa:** HG005 72 against 63. The rest comes from
  BWA-MEM's gap costs, which CERTA shares (the indel bullets of [Next](#next)).

**Tests.** [test_certa.cpp](../tests/test_certa.cpp) `test_tier_r_wide`:
- **Reads:** 160 reads with an insertion or a deletion of 9–40 bases at
  least 30 bases from either end, 120 of them from a 400-bp element with 20
  identical copies.
- **Without the wide band**, lo reaches OPT on 32–43 of the 120 element reads,
  about those whose indel is at most 16 bases.
- **With it**, lo reaches OPT on all 120, and every read it raises is aligned
  through its indel.
- **Unique sequence:** for the reads from unique sequence, the chains already
  span the indel.
- **Invariants:** the bounds stay sound against the oracle on every read, hi
  is max(lo, hi without it), and the parts are the same.
- **Default configurations:** the wide band is also on in the interval-mode
  configurations of the tier R, pair-alternative, joint-seed,
  pair-placement and `certa check` tests.

All checks pass on seeds 20260927, 1, 2 and 11–15 (CI's), and so do the
lemma checks. Mutants WB1–WB4 (the band never aligned, not centred on the
band of A*, lo raised without its band, the threshold never met) are all
killed on seeds 1 and 2 ([tests/mutants.tsv](../tests/mutants.tsv)).

## Next

1. **Accuracy.** With `--iv-pair all --mapq pair`, CERTA is the most
   accurate aligner measured at every operating point, held out on HG005;
   with `--t2t` by 2.2–14× over minibwa
   ([above](#the-pair-model-refit-on-the-default-builds-placements)).
   - *Done (2026-11-19):* the fitted MAPQ cost DeepVariant 1.8–3.5×
     minibwa's errors because it was written as a Phred; on BWA-MEM's scale
     it gives fewer errors than the default MAPQ
     ([above](#the-fitted-mapq-on-bwas-scale)).
   - *Done (2026-11-20):* both models refitted on the default build's
     placements; XM is not an input
     ([above](#the-pair-model-refit-on-the-default-builds-placements)).
     Next: make `--mapq pair` the default MAPQ of paired interval mode,
     together with `--iv-pair all --iv-rescue local`, the records it was
     fitted on.
   - **Indels** (HG005, the refit `--mapq pair` run;
     [hg005_indel_sites.md](results/truth_2026-11_refit/hg005_indel_sites.md)).
     CERTA makes 18 indel errors that minibwa does not; minibwa makes 7 that
     CERTA does not.
     - 11 of the 18 are BWA-MEM2's errors too. BWA-MEM's single affine gap
       (6 + k, CERTA's objective) makes clipping a long insertion score more
       than aligning through it, where minibwa's two-piece gap costs do not.
       For example, a 32-bp insertion in a TG repeat scores 72 aligned
       through and 86 clipped. Aligning those reads needs another objective,
       and with it new certificates.
     - About 4 are search misses under CERTA's own objective. Tier R kept
       seed diagonals more than 16 apart in separate bands, so reads across
       a 24-bp insertion were clipped (68M83S, AS 63) where the alignment
       through it scores 97 (63M24I64M). *Done (2026-11-20):* tier R's wide
       band ([above](#tier-rs-wide-band)) removes three of the four.
     - Next for indels: minimap2-style two-piece gap costs, as an option with
       their own certificates, where clipping a long insertion now scores
       more than aligning through it.
   - **False SNV calls without `--t2t`** (2026-11-21;
     [snv_fp_diagnosis.md](results/truth_2026-11_fp/snv_fp_diagnosis.md)).
     With `--mapq pair` CERTA makes 364 against minibwa's 309 on HG002, and
     113 against 71 on HG005.
     - Of the calls only CERTA makes (127 and 73), BWA-MEM2 makes 36 and 10,
       so they do not come from the shared objective.
     - At those sites, most reads with the false allele are reads minibwa
       also counts with MAPQ ≥ 5 (662 of 1,123 on HG002). Many of the rest
       are PR mates in two-copy duplications.
     - *Tried, not adopted:* writing the model's Phred 3–4 below
       DeepVariant's MAPQ floor of 5, where minibwa leaves more of its least
       confident mates. Missed SNVs rose by 18–69, and false calls barely
       moved (HG002 364 → 369 and 358; HG005 113 → 110 and 114).
     - The calls differ in which reads each aligner trusts around paralogs.
       With `--t2t`, CERTA's errors equal or beat minibwa's.
   - On GRCh38 the placement limit is reached
     ([above](#how-far-from-the-limit)); what remains is the reference.
     The sample's own genome shows two orders of magnitude of room in
     confidence: CHM13's copies of GRCh38's collapsed duplications as sink
     contigs, then pangenome haplotypes, are the way there.
   - With `--iv-rescue local` and the seed votes, the default MAPQ is 118
     SNV and 115 indel errors behind minibwa on HG002
     ([above](#reads-certa-lost)). Before the votes, nearly all of the gap
     was the MAPQ of reads CERTA places where minibwa does: with minibwa's
     MAPQ there, 24 and 10 ([above](#a-local-mate-rescue)). The PR mates'
     MAPQ from the pairs found closes part of it, but costs the ranking
     unless T2T-CHM13 gates it ([above](#the-mapq-of-mates-placed-by-the-pair)).
   - 51,787 reads that minibwa places with MAPQ ≥ 20 are still elsewhere in
     CERTA: 24 k search misses (14 k with ≥ 10 differences) and 28 k
     clipped or chimeric reads that score as well where CERTA puts them.
     Score them against the truth, not against minibwa.
   - A second held-out genome (HG001 or another HPRC sample; a download).
   - The truth labels of tests B and H are certified unique on one
     haplotype, which under-samples segmental duplications (item 6); a
     model fitted on them may rank those mates less well.
   - `--t2t` holds two indexes (45.9 GiB). CHM13's copies of GRCh38's
     collapsed duplications as extra contigs of one index (FixItFelix's
     design; [the research report](../reports/Biologically%20correct%20read%20placement%20ideas.md))
     might give most of the gain at one index's memory.
   - Score accuracy at the locus and at matched sensitivity from now on.
2. **Cost.** The edit certificate that `--bio` runs first is now 48 % of
   its time; its band DP (`band_distance_vec8`) is also the largest cost of
   `--full-cert`.
3. **MAPQ of uncertified records.**
   - Test `--iv-mapq cert` in variant calling before making it the
     default: XQ for every record dropped the reads that carry variants
     (FINDINGS.md, Aim 3).
   - Place the reads `--bio` cannot certify with interval mode's tiers in
     the same process, so that `--bio` needs no second aligner either.
4. **MAPQ.** Replace the one-competitor-at-`YH` convention with a bound
   that does not make horizon-limited reads overconfident.
5. **Complex indels.** Allow an insertion next to a deletion in the model
   (`bio_band`, Lemma B, the oracle test).
6. **Test B labels.** Label PR-tier mates from the pair certificate, so the
   real-read truth covers segmental duplications.
7. **Phase 2 for `--bio`.** HG005's truth is built and interval mode is
   scored on it ([above](#held-out-hg005)); `--bio` has not been run on it.
   Run it once, with the build frozen after any change from 1–6.
