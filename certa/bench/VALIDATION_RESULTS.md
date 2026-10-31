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

## Next

1. **Fallback.** Lower the share of reads `--bio` sends to the fallback
   (`bio_horizon`, 4.7 % of reads) and its cost. Candidate enumeration
   dominates the time.
2. **MAPQ.** Replace the one-competitor-at-`YH` convention with a bound
   that does not make horizon-limited reads overconfident.
3. **Complex indels.** Allow an insertion next to a deletion in the model
   (`bio_band`, Lemma B, the oracle test).
4. **Test B labels.** Label PR-tier mates from the pair certificate, so the
   real-read truth covers segmental duplications.
5. **Phase 2.** HG005 (download 1.8 GB), the held-out verdict, with the
   build frozen after any change from 1–4.
