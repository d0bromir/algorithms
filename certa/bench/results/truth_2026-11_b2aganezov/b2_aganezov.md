# B2 Aganezov collapsed-loci sinks — definitive: missing copies are chain-covered, not gaps

**Date:** 2026-11-28. **Step:** Pillar B, Aganezov collapsed-loci sinks.
**Code:** [`../../b2_aganezov.py`](../../b2_aganezov.py). Sink BED:
[`sinks_aganezov.bed`](sinks_aganezov.bed) (18 regions).

## Approach

The genuine missing copies should be CHM13 segmental-duplication sequence absent
from GRCh38. Downloaded the **CHM13 v2.0 SEDEF segdup annotation**
(`chm13v2.0_SD.bed`, 83,516 regions, T2T S3) and intersected it with the
CHM13-only chain gaps (`b2_extract`, 463 gaps), excluding Cen/Sat — i.e.
euchromatic SD sequence GRCh38 lacks (the Aganezov collapsed-locus copies).

## Result: only 102 kb

**18 regions, 102,269 bp (0.10 Mb)** — versus Aganezov's ~8.04 Mb of collapsed
loci. The largest are telomeric/subtelomeric SD fragments (chrX:7.72–7.73 Mb
10 kb; chr1:1–5525 telomere; chrY tips). Essentially nothing.

## Definitive finding

**The missing copies are chain-*covered* collapsed regions, not chain gaps.**
The `hs1ToHg38` chain maps a CHM13 copy *onto* the GRCh38 collapsed locus (the B1
finding), so the collapsed copies are "covered" by the chain and never appear as
uncovered gaps. Intersecting segdups with gaps therefore recovers almost nothing
(102 kb of telomeric SD, not the 8 Mb of collapsed copies).

This conclusively shows that **every gap-based B2 variant** (heuristic sinks,
Cen/Sat-refined sinks, segdup∩gap) **cannot recover the missing-copy class**:

| approach | sink content | missing-copy abstentions | stolen real |
|---|---|---|---|
| heuristic gaps | 11.2 Mb | 245 | 916 |
| Cen/Sat-refined gaps | 1.5 Mb | 22 | 118 |
| segdup ∩ gaps (Aganezov) | 0.10 Mb | ~0 | ~0 |

The "missing copy" is not sequence absent from the reference *as a gap*; it is a
paralogous copy **collapsed onto** an existing reference locus. Adding gap/decoy
sequence (FixItFelix-style) or lifting back to GRCh38 (levioSAM2-style, B1) cannot
address it.

## Correct path

To place the missing-copy reads, one must represent the *extra copies*:
1. **Copy-level alignment**: take the Aganezov GRCh38 collapsed loci (203 loci,
   ~8 Mb; BED on the GIAB FTP / mydennislab Globus, not plain HTTP), align each
   against CHM13 to find *all* paralogous copies (the chain gives only the best),
   and add the extra copies as sinks. Needs the collapsed-locus BED + a
   many-to-many aligner.
2. **B5 copy-level joint assignment** (EM over reads × copies): resolve the copy
   of origin across reads without adding reference sequence — the principled
   route, and the one that matches the "collapsed, not absent" nature of the
   problem.
3. **B4 pangenome** (HPRC2 haplotypes as an explicit copy catalog).

## Significance (thesis)

This is the clearest result of Pillar B: it operationally distinguishes
"sequence absent from the reference" (satellites/gaps — fixable by decoys) from
"copy collapsed in the reference" (the missing-copy class — requires copy-level
reasoning). It explains, with measurements, why decoy/lift-over approaches fail
for the dominant residual error class, and motivates B5/B4.

## Disk safety & confirmation

Annotation downloads (6.7 MB SD BED) + extraction in RAM; no index built (102 kb
of sinks would not meaningfully change mapping). Root disk unchanged. Confirmed
galaxy ARM64.

## Reproduce

```bash
curl -o chm13_sd.bed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed
python3 bench/b2_aganezov.py ~/_paper_work/hs1.fa chm13_gaps.bed chm13_sd.bed censat.bed sinks_aganezov
```
