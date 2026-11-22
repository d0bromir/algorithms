# B2-extract — CHM13-only sink candidates: characterization

**Date:** 2026-11-22. **Step:** Pillar B, B2 (one-pass sink contigs), extraction
analysis. **Code:** [`../../b2_extract.py`](../../b2_extract.py). Gap BED:
[`chm13_gaps.bed`](chm13_gaps.bed) (1,117 gaps ≥ 5 kb).

## Method

Parse `hs1ToHg38.over.chain`, compute the CHM13 (hs1) intervals it covers, and
report the uncovered gaps — sequence present in CHM13 but absent from GRCh38
(the sink candidates). Read-only on the chain; no reference/index built.

## Result

- **Total CHM13-only gap bp: 278.1 Mb.** Gaps ≥ 5 kb: **1,117 regions, 273.4 Mb**.
- Size histogram (all gaps): `<1k: 351,797` (small indels), `1–5k: 914`,
  `5–50k: 681`, `50–500k: 397`, `≥500k: 39`.

**The raw sink content is dominated by centromeric satellites, not missing
copies.** The 39 gaps ≥ 500 kb (~200 Mb) are:

| gap | size | identity |
|---|---|---|
| chrY:26.6–56.8 Mb | 30.2 Mb | Y heterochromatin |
| chr15:0–19.8 Mb | 19.8 Mb | acrocentric short arm / satellite |
| chr14:0–18.2 Mb | 18.2 Mb | acrocentric short arm / satellite |
| chr13:0–18.2 Mb | 18.2 Mb | acrocentric short arm / satellite |
| chr1:125.2–143.2 Mb | 18.0 Mb | centromere |
| chr9:43.3–60.5 Mb | 17.2 Mb | centromere / heterochromatin |
| chr22:0–10.5 Mb | 10.5 Mb | acrocentric short arm |
| chr16:36.3–46.4 Mb | 10.1 Mb | centromere |
| chr18/21/17/X/11/5… | 3–5.5 Mb | centromeres / satellites |

These are unresolved centromeres/satellites in GRCh38, **not missing gene
copies**; adding them as sinks would attract reads spuriously (the FixItFelix
FANCD2 failure mode). They must be excluded.

## Sink-candidate estimate

Excluding the ≥500 kb satellite gaps leaves ~73 Mb of moderate gaps (5–500 kb),
but many of those are still satellite fragments or assembly gaps. The genuine
**missing-copy** content is the subset flanked by segmental duplications —
Aganezov's 203 collapsed loci (~8.04 Mb) and Vollger's CHM13-only SDs (~16–81 Mb).
So the useful sink content is **~8–80 Mb**, far below the raw 278 Mb.

## B2-build requirements (next step)

1. **Satellite exclusion**: use the T2T CHM13 Cen/Sat annotation (v2.1) to mask
   centromeric/satellite gaps (not on galaxy yet — a download); heuristic fallback
   = drop ≥500 kb gaps + acrocentric contig-start gaps (chr13/14/15/21/22 at 0) +
   chrY heterochromatin.
2. **Missing-copy filter**: keep moderate gaps flanked by segdups (GIAB segdup
   BED / Aganezov collapsed-loci BED); alt-mask near-identical orthologs so a sink
   does not steal reads from a real GRCh38 copy.
3. **Append as ordinary contigs** (not ALT — so they depress the competing GRCh38
   hit's MAPQ), build the augmented CERTA index, map; sink-winning mates → MAPQ 0
   + sink-coordinate tag; certificate over GRCh38 ∪ sinks.
4. **Per-locus validation** against HG002 truth (every sink must not steal real
   reads); freeze HG002 / hold out HG005.

**Disk cost:** sinks ~8–80 Mb add ~0.05–0.5 Gb to the index (GRCh38 s1v3 ≈ 20.4 Gb
for 3.1 Gb ⇒ ~6.6 bytes/base) plus the augmented fasta. Augmented index ≈ 21 Gb;
galaxy `/` has ~41 Gb free — feasible but tight; build in a scratch dir and remove
the un-augmented intermediate. The Cen/Sat annotation download is small.

## Go / no-go

- B2-extract confirms the sink approach is viable **but satellite exclusion is
  essential** (raw gaps are ~70% satellite). Proceed to B2-build with the Cen/Sat
  annotation + segdup-flanking filter; the useful sink content (~8–80 Mb) is small
  enough to index within the disk budget.

## Disk safety & confirmation

Read-only chain analysis on galaxy; the gap BED (37 KB) is the only output. No
reference/index built; no disk pressure. Confirmed galaxy ARM64 (numpy 2.4.4).

## Reproduce

```bash
python3 bench/b2_extract.py ~/headroom/hs1ToHg38.over.chain.gz chm13_gaps.bed --min-gap 5000
```
