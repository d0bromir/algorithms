# B2 refinement — Cen/Sat filtering removes stealing *and* benefit: chain gaps are not missing copies

**Date:** 2026-11-28. **Step:** Pillar B, B2 refinement. **Code:**
[`../../b2_refine.py`](../../b2_refine.py), [`../../b2_refine_run.sh`](../../b2_refine_run.sh).
Results: [`b2_refine_report.md`](b2_refine_report.md),
[`b2_refine_sinks.md`](b2_refine_sinks.md), [`sinks_refined.bed`](sinks_refined.bed).

## What was done

Downloaded the **T2T CHM13 Cen/Sat v2.1** annotation (2,523 centromeric
satellite / higher-order-repeat / monomeric / rDNA / ct regions) and excluded
sink candidates overlapping it (plus the size/acrocentric/chrY heuristic).
Rebuilt the augmented index in RAM and re-mapped the 4 M-pair subsample.

- **Refined sinks: 174 contigs, 1.50 Mb** (from 269 / 11.2 Mb). Excluded 194
  heuristic-satellite + 95 Cen/Sat-overlapping gaps. Sizes: 5–20 kb ×167,
  20–100 kb ×7 (no ≥100 kb — those were all satellite). Euchromatic candidates.

## Result: stealing falls, but so does the benefit

| sink set | sinks | sink placements | stolen real (confident) | missing-copy abstain | MAPQ≥40 err Δ |
|---|---|---|---|---|---|
| heuristic (step 2) | 269 / 11.2 Mb | 82,698 | 916 (266) | 245 | −30 |
| **Cen/Sat-refined** | 174 / 1.5 Mb | 1,676 | **118 (36)** | **22** | **−1** |

The refined sinks steal far fewer real reads (916→118; confident 266→36) — the
FixItFelix failure mode is largely fixed — **but they also abstain far fewer
missing-copy reads (245→22) and reduce confident errors by only 1 (was 30).**
Placement is otherwise identical (at-truth 20.085→20.084).

## The clarifying finding

**Chain-uncovered gaps are mostly satellites, telomeres and assembly gaps — not
missing copies.** Removing the satellites (Cen/Sat) removes most of the sink
content *and* most of the benefit. The genuine missing copies (collapsed segmental
duplications) are **chain-*covered***: the `hs1ToHg38` chain maps a CHM13 copy onto
the GRCh38 *collapsed* locus (the B1 finding), so they are not "uncovered gaps" and
the gap-based extraction never sees them. What the gap extraction captures is the
sequence GRCh38 genuinely lacks *as sequence* (centromeres/satellites/telomeres),
which is not where the missing-copy reads come from.

So **B2 via chain-gaps is fundamentally limited**: it adds the wrong sequence.
The missing-copy reads originate from copies that exist in CHM13 *within
chain-covered, collapsed* regions.

## Correct path (next)

- **Aganezov collapsed-loci sinks**: take the 203 collapsed loci (~8.04 Mb,
  GRCh38 coords), lift them to CHM13 via the chain, and extract *those* CHM13
  regions as sinks (not the chain gaps). These are the genuine missing copies.
  Alt-mask near-identical orthologs; validate per-locus.
- Or **B4 pangenome** (HPRC2 haplotypes as the copy catalog), which represents the
  collapsed copies explicitly.
- This negative/clarifying result is thesis-worthy: it distinguishes
  "sequence absent from the reference" (satellites/gaps) from "copy collapsed in
  the reference" (the missing-copy class), and shows why naive decoy/sink addition
  (chain gaps) fails while collapsed-locus-targeted sinks are the right unit.

## Disk safety & confirmation

Annotation download (195 KB) + index/SAMs/npz in `/dev/shm` (RAM); root disk
unchanged (~40 Gb free). Confirmed galaxy ARM64; a2 not run (no data/numpy).

## Reproduce

```bash
curl -o censat.bed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed
python3 bench/b2_refine.py ~/_paper_work/hs1.fa chm13_gaps.bed censat.bed sinks_refined
bash bench/b2_refine_run.sh
```
