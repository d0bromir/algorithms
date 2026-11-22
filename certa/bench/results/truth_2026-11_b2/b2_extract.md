# B2-extract — CHM13-only sink candidates: characterization

**Date:** 2026-11-22. **Step:** Pillar B, B2 (one-pass sink contigs), extraction.
**Code:** [`../../b2_extract.py`](../../b2_extract.py). Gap BED:
[`chm13_gaps.bed`](chm13_gaps.bed) (463 gaps ≥ 5 kb).

## Method and orientation correction

Parse `hs1ToHg38.over.chain` and report the uncovered **CHM13 (hs1)** intervals —
sequence present in CHM13 but absent from GRCh38 (the sink candidates).

**Orientation:** the chain header has `tName/tSize` = CHM13 (e.g. chr2 tSize
242,696,752) and `qName/qSize` = GRCh38 (chr2 qSize 242,193,529), i.e. the chain
is oriented **hg38(query) → hs1(target)**. CHM13-only sequence is the uncovered
**target** (hs1) intervals, so the parser tracks `tName/tStart` and `dt`. (An
earlier version tracked the query and reported hg38-only gaps — 278 Mb of hg38
patches/centromeres — which is corrected here.)

## Result

- **Total CHM13-only gap bp: 212.4 Mb.** Gaps ≥ 5 kb: **463 regions, 206.0 Mb**.
- Histogram (all gaps): `<1k: 425,799` (indels), `1–5k: 1,302`, `5–50k: 304`,
  `50–500k: 109`, `≥500k: 50`.

**Most raw gap content is centromeric satellite, not missing copies.** The
largest gaps are chrY heterochromatin (34.6 Mb), chr9 centromere/heterochromatin
(27.9 Mb), chr16 centromere (10.0 Mb), chr4/chr18/chr1 centromeres (5–5.6 Mb),
and acrocentric short arms (chr13:0–5 Mb). These must be excluded (FixItFelix
FANCD2 failure mode — a sink that steals real reads).

## Sink-candidate estimate

After satellite exclusion (B2-build step 1) the useful missing-copy content is
**~11 Mb** (269 contigs), consistent with Aganezov's ~8 Mb collapsed-loci
estimate. See [`b2_build.md`](b2_build.md).

## Disk safety & confirmation

Read-only chain analysis on galaxy; the gap BED (13 KB) is the only output. No
reference/index built. Confirmed galaxy ARM64 (numpy 2.4.4).

## Reproduce

```bash
python3 bench/b2_extract.py ~/headroom/hs1ToHg38.over.chain.gz chm13_gaps.bed --min-gap 5000
```
