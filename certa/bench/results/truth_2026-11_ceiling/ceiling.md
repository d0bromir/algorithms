# A1 — the offline ceiling of a gap-graded MAPQ (HG002 test B)

**Date:** 2026-11-21. **Step:** research guideline Pillar A, A1 (the offline
ceiling test; go/no-go for the reference-content work in Pillar B).
**Code:** [`../../truth_ceiling.py`](../../truth_ceiling.py),
[`../../test_truth_ceiling.py`](../../test_truth_ceiling.py),
[`../../run_ceiling.sh`](../../run_ceiling.sh). Raw tables:
[`ceiling_holdout.md`](ceiling_holdout.md) (held out), [`ceiling_full.md`](ceiling_full.md) (in-sample).

## Question

How much of CERTA's residual confident misplacement can a *calibrated
competitor-gap signal* close on its own, before investing in reference
augmentation (Pillar B)? The placement report infers that ~85% of confident
errors sit where BWA-MEM2's score gap `AS−XS < 50`; A1 turns that inference
into a measured ceiling.

## Method

- Stream **BWA-MEM2** (`bwa-mem2 mem -t 64`, GRCh38) over HG002 test B
  (20 M pairs) directly into `truth_ceiling.py feat` — **no BAM written** —
  extracting per-mate `gap = AS−XS` (its `XS` is a *local* score, so it counts
  clipped partial paralogs that CERTA's end-to-end second-best misses).
- Join with CERTA **default mode**'s scored truth (`certa.npz`: placed /
  at-locus-150 bp / MAPQ) and `truth.npz`. 31,931,012 mates have a GRCh38 truth.
- `fit`: bin the gap, run weighted **PAVA** (isotonic, non-increasing error in
  gap), write `gap → MAPQ` (Phred of the pooled misplaced share, +0.5/+1, ≤60).
  **Fit on even pairs** (15,955,063 mates; 13,423 misplaced).
- `eval`: re-rank by `min(CERTA MAPQ, gap MAPQ)`; report misplaced per million
  at matched N (N most confident, boundary MAPQ pro rata) **on the odd pairs
  held out** (15,958,866 mates), and at each run's own MAPQ ≥ 40.
- Fitted gap MAPQ (even pairs): gap 0→3, 1→6, 2→7, 5→8, 10→18, 20→24, 50→32,
  100→38, no-competitor→58.

## Results (odd pairs held out)

**Matched N — misplaced per million among the N most confident mates:**

| run | 99.75% | 99.29% | 98.00% | 95.00% |
|---|---|---|---|---|
| baseline MAPQ | 190.4 | 141.0 | 108.1 | 108.1 |
| **min(MAPQ, gap MAPQ)** | 191.2 | 120.2 | **52.0** | **31.7** |
| gap MAPQ alone | 199.1 | 136.9 | 78.9 | 55.4 |

**At each run's own MAPQ ≥ 40 (the operating point):**

| run | confident % | misplaced per M |
|---|---|---|
| baseline MAPQ | 99.49 | 159.7 |
| **min(MAPQ, gap MAPQ)** | 77.89 | **11.2** |
| gap MAPQ alone | 77.89 | 15.2 |

**Diagnostic:** of 2,535 confident (MAPQ ≥ 40) misplaced mates, **1,971 (77.8%)
have competitor gap < 50** and 334 (13.2%) gap < 5. The full-data (in-sample)
tables are consistent (baseline 184.9/M at N=99.75%; min 31.6/M at N=95%).

## Interpretation

1. **The gap signal carries most of the error information.** 77.8% of confident
   errors have BWA gap < 50 (the report's marginal-stratification inference of
   ~85% is confirmed in spirit on held-out data). Re-ranking breaks the baseline
   MAPQ's saturation: the baseline is flat at 108.1/M for N=98% and 95% (its
   MAPQ-60 tie group cannot be ranked within — Limit 2), whereas the gap-graded
   score falls 108 → 52 → 32/M as N drops 98 → 95%.
2. **Calibrated abstention reaches ~11/M — near the ~3/M own-genome confidence
   ceiling** (Limit 4): a **14× reduction** in confident misplacement
   (159.7 → 11.2/M). The headroom the guideline claims is real and measurable.
3. **But gap-alone over-abstains.** The 14× gain costs **~22% sensitivity**
   (confident fraction 99.49% → 77.89%): the isotonic gap MAPQ caps at ~38–58,
   so it demotes *genuine-repeat correct* reads together with *missing-copy
   errors*. At matched N=99.75% (≈ the baseline's confident set) re-ranking
   does **not** help (191.2 vs 190.4) — keeping that sensitivity forces the
   errors back in. The gap cannot, by itself, tell "gap < 50 because a copy is
   missing" (error) from "gap < 50 because a real repeat" (correct).

## Go / no-go

- **Go for A2–A4.** The gap signal is validated; the next lever is to keep its
  error discrimination while recovering sensitivity — i.e. a **locus prior +
  phantom-copy posterior** (A3b/A4) that demotes missing-copy loci specifically,
  not all small-gap loci. A naive `min(MAPQ, gap)` is *not* the production rule.
- **Go for Pillar B.** The ~22% abstained mates are predominantly the
  missing-copy class; converting those abstentions into *correct placements*
  requires reference content (CHM13 v2.0 sinks / lift-and-commit, B1/B2). A1
  shows calibration alone approaches the ceiling but cannot place the missing
  copies — exactly the division of labour the guideline assumes.
- **Baseline caveat:** this used CERTA *default* mode (`certa.npz`), not the
  interval `--iv-mapq cert` production MAPQ (its scored npz was not present on
  the host). The ceiling conclusion is about the *gap signal*, which is
  baseline-agnostic; re-running on `--iv-mapq cert` / `--mapq pair` is a trivial
  follow-up (`feat` takes any scored npz).

## Disk safety (galaxy)

BWA-MEM2's SAM was streamed into the extractor (no BAM); every intermediate
(`joined.npz`, models, reports) lived in `/dev/shm` (RAM). Root disk stayed at
**44 G free (96%)** before and after — no spill. Index and FASTQ were read-only.

## Confirmation

- **Unit tests 14/14** on Windows x64 (Python 3.14.6 / numpy 2.5.0) **and galaxy
  ARM64** (numpy 2.4.4): PAVA monotonicity + weighted-total conservation, Phred,
  matched-N pro-rata ties, SAM gap extraction (missing XS / unmapped / secondary),
  and the error-demotion mechanics.
- **Data run confirmed on galaxy** (the truth data lives there).
- **a2 not run:** a2 has no truth data, no numpy/pip/ensurepip in any python,
  and its root disk is **100% full (6.9 G)**; obtaining numpy needs
  `sudo apt install python3.12-venv` (disk writes + password), which the
  disk-safety rule forbids. Code confirmation therefore rests on Windows +
  galaxy (ARM64); the code is pure numpy (arch-independent).

## Reproduce

```bash
# on galaxy, everything in RAM:
mkdir -p /dev/shm/ceiling && cp ~/certa-bio/certa/bench/truth_lift.py /dev/shm/ceiling/
cp bench/truth_ceiling.py bench/run_ceiling.sh /dev/shm/ceiling/   # + sed -i 's/\r$//'
bash /dev/shm/ceiling/run_ceiling.sh        # feat (stream bwa-mem2) + fit + eval
python3 -m unittest test_truth_ceiling      # from certa/bench
```
