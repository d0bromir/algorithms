# Paired production — the phantom prior in `bio_pair` / `bio_rescue` (HG002 test B)

**Date:** 2026-11-22. **Step:** paired production (extends A3–A5 from single-end
to the production paired MAPQ). **Code:** `include/certa/bio.h`
(`BioOut::phantom_w`, `bio_place` stores it, `bio_pair`/`bio_rescue` use it),
`tests/test_certa.cpp`. Raw table: [`paired_report.md`](paired_report.md).
Runner: [`../../paired_run.sh`](../../paired_run.sh).

## Change

A3–A5 added the phantom-copy term to the **single-end** posterior (`bio_place`).
Production `--bio` reports the **paired** MAPQ (`bio_pair`, with `bio_rescue`
fallback), so the gain did not yet reach production. Now:
- `BioOut` carries `phantom_w` — the weight `bio_place` adds to the single-end
  posterior (locus-specific via the A4 track).
- `bio_pair` adds each mate's `phantom_w` to that mate's pair posterior
  `(total−num)/num`, demoting only the divergent mate.
- `bio_rescue` adds the certified mate's stored weight and computes the rescued
  mate's from its traced alignment (`bio_align_locus` + `bio_mm_counts` + locus
  prior).

MAPQ only — the pair ML certificate is unchanged; phantom off ⇒ byte-identical.

## Correctness (galaxy ARM64, RAM build)

- `verify_lemmas` + `test_certa` **pass** (FAIL=0).
- New `bio_pair` unit test: with mate 1's `phantom_w=0.5`, mate 2's 0, the pair
  MAPQ is **5/60** (mate 1 demoted, mate 2 confident); with no phantom, **60/60**.
- Real data: **placements identical** (at-truth 20.061%, 5,004 misplaced both
  runs) — the pair certificate is untouched.

## Results (4 M pairs = 8 M mates, both scored; A5-frozen δ=0.006, π_∅=0.003 in segdups, 0 outside)

| run | MAPQ≥60 mispl. | MAPQ≥40 mispl. | MAPQ≥30 mispl. | MAPQ≥20 mispl. | MAPQ≥40 placed |
|---|---|---|---|---|---|
| paired `--bio` | 627 | 909 | 1,521 | 1,770 | 6,134,256 |
| paired phantom | **242** | **356** | **586** | 1,614 | 5,980,836 |

- **MAPQ≥40 errors −61% (909 → 356)**, MAPQ≥60 −61% (627 → 242), MAPQ≥30 −61%
  (1,521 → 586), MAPQ≥20 −9% (1,770 → 1,614), at a **~2.5% sensitivity cost**
  (153 k mates demoted below 40). Stronger than single-end (−49–59%) because
  **both mates** benefit.
- **Strata (MAPQ≥20):** segdup 1,275 → **1,141 (−11%)**, lowmap_segdup 1,487 →
  1,351 (−9%), difficult 1,606 → 1,470 (−8%); **easy 13 → 13 (unchanged)**.
- **Calibration:** the MAPQ 50–60 band's overconfidence falls from 763/6.6 ≈ 116×
  to 304/6.4 ≈ 47×.
- Missing-copy class ("no linear truth"): MAPQ≥20 429 → 384 (−10%).

## Interpretation

- The calibrated, locus-specific phantom prior now reaches the **production
  paired MAPQ**, cutting confident misplacement ~61% with the certificate intact
  and easy reads unaffected. This completes Pillar A (A1 ceiling → A2 local
  competitor signal → A3 phantom term → A4 locus prior → A5 calibration → paired
  production): a principled, certificate-safe, cross-sample-validated posterior.
- The remaining confident errors are predominantly the missing-copy class whose
  true copy GRCh38 lacks; placing them correctly (not just abstaining) is
  **Pillar B** (CHM13 v2.0 lift-and-commit / one-pass sinks). Speed is **Pillar C**.

## Go / no-go

- **Paired production validated** (keep; flag-gated, default byte-identical).
- Next: **Pillar B** (reference content to place the missing copies) and
  **Pillar C** (AVX-512 / GPU tier R / CERTA-X exact completion). Optional
  refinements: a sharper copy-deficit track (Aganezov/CHM13), conformal calibration.

## Disk safety & confirmation

- RAM build/run (`/dev/shm/a2`, `TMPDIR` in RAM, subsamples/BED in RAM, SAM
  streamed, no BAM); galaxy `/` ~96% from shared-host activity, this work wrote
  nothing to it. Confirmed galaxy ARM64; a2 not run (no truth data/numpy, disk
  full). C++ arch-independent; certificate tests are the cross-platform guarantee.

## Reproduce

```bash
bash bench/paired_run.sh   # fetch branch, build, tests, paired --bio +/- phantom (galaxy)
```
