# B1-impl — selective lift-and-commit (`--t2t-commit`): a negative result

**Date:** 2026-11-22. **Step:** Pillar B, B1 implementation. **Code:**
`sam_output.cpp` (`commit_lift`, `pair_interval`), `sam_output.h`
(`OutputOptions::t2t_commit/t2t_commit_mapq`), `map_config.cpp`, `cli.h`.
Measurement: [`b1impl_report.md`](b1impl_report.md); diagnostics:
[`b1impl_diag.md`](b1impl_diag.md). Runner: [`../../b1impl_run.sh`](../../b1impl_run.sh).

## What was implemented

`commit_lift` (levioSAM2-style commit/defer): a confident T2T-CHM13 lift
(`state==1`, MAPQ ≥ `--t2t-commit-mapq`, default 40) that lands *elsewhere* than
the GRCh38 placement proposes the missing copy's projected locus; the read is
re-aligned there on GRCh38 (`affine_align`, strand-encoded as the existing mate
rescue does) and committed **only if it scores better**. The committed mate
becomes `kTierMR` (moved, not certified — a lift is a proposal verified on the
reference, not a global certificate). Applied early in `pair_interval` so both
output paths use it. Flag-gated (`--t2t-commit`, needs `--t2t`); default
byte-identical. Certificate tests pass (FAIL=0).

## Result: ~0 benefit (negative)

Full HG002 test B subsample (4 M pairs), paired interval `--mapq pair --t2t`:

| run | at truth % | misplaced | MAPQ≥40 | MAPQ≥60 |
|---|---|---|---|---|
| `--t2t` | 20.085 | 6,115 | 1,232 | 229 |
| `--t2t --t2t-commit` | 20.085 | 6,115 | 1,232 | 229 |

**Byte-identical.** Diagnostics (500 k pairs, ~900 k `commit_lift` calls):

| outcome | share |
|---|---|
| lift not confident / unmapped / not certified ("early") | ~12% |
| **lift lands "here" (within 150 bp of the GRCh38 placement)** | **~89%** |
| elsewhere, re-aligned, but score ≤ GRCh38 (defer) | ~0.24% |
| **committed** | **2 / 900 k** |

## Diagnosis

The `hs1ToHg38` chain projects a CHM13 true-copy placement **back to the same
collapsed GRCh38 locus** CERTA already placed at — so ~89% of lifts are "here"
and propose nothing new. For the few that land elsewhere, re-aligning on GRCh38
cannot score better, because **GRCh38 lacks the true-copy sequence** (that is the
missing-copy problem itself). Lift-and-commit therefore cannot rescue
missing-copy reads: lifting back to a reference that lacks the copy collapses
the answer again.

**Reconciling the B1 ceiling:** the 38.9% "rescue" measured in B1 was
*BWA-MEM2-on-CHM13* (a different aligner with different placements) scored
against the assembly truth. For CERTA — whose GRCh38 placements are already
99.86% at-truth — the CHM13 lift lands "here", so there is no better locus to
commit to. The ceiling did not transfer to a lift-and-commit implementation.

## Conclusion / pivot

- **B1 lift-and-commit is ineffective for GRCh38 missing-copy rescue.** Keep
  `--t2t-commit` as a correct, harmless, flag-gated capability (it may help with
  a finer lift chain — CHM13 v2.0 + nf-LO — or other references), but it is not
  the effective path here.
- **The effective path is B2 — one-pass sink contigs**: append the CHM13-only
  (unliftable/collapsed) sequence as ordinary contigs so missing-copy reads map
  to the *true copy* with a genuinely better score (and abstain in GRCh38
  coordinates with a sink tag). Adding the sequence, not lifting back to a
  reference that lacks it, is what recovers the missing copies.
- This negative result is itself a thesis-worthy finding: it shows why
  lift-over-based rescue (levioSAM2) does not address the *absent-copy* class,
  motivating reference augmentation (sinks/pangenome) instead.

## Disk safety & confirmation

RAM build/run (`/dev/shm/a2`, streamed SAM, no BAM); galaxy `/` ~96% from
shared-host activity, this work wrote nothing to it. Confirmed galaxy ARM64
(tests + measurement); a2 not run (no data/numpy, disk full).

## Reproduce

```bash
bash bench/b1impl_run.sh   # fetch branch, build, tests, paired --t2t +/- --t2t-commit (galaxy)
```
