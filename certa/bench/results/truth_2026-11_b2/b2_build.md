# B2-build step 1 — filter and extract CHM13 sink candidates

**Date:** 2026-11-22. **Step:** Pillar B, B2 (one-pass sink contigs), build step 1.
**Code:** [`../../b2_build.py`](../../b2_build.py). Sink BED:
[`sinks.bed`](sinks.bed) (269 contigs). Sink FASTA (`sinks.fa`, 11.3 Mb) lives on
galaxy `/dev/shm/a2/sinks.fa` and is regenerable (not committed — reference
sequence).

## Result

From the 463 CHM13-only gaps ≥ 5 kb (206 Mb, [`b2_extract.md`](b2_extract.md)),
satellite exclusion drops 194 gaps and keeps **269 sink candidates, 11.2 Mb** —
consistent with Aganezov's ~8 Mb collapsed-loci estimate.

- **Filter (heuristic):** drop gaps ≥ 500 kb (satellite arrays), acrocentric
  short arms (chr13/14/15/21/22 start < 20 Mb), chrY heterochromatin (> 26 Mb).
- **Sink size distribution:** 5–20 kb: 212, 20–100 kb: 25, 100–500 kb: 32.
- **Largest sinks** (pericentromeric / segdup regions): chr1:127.2–127.7 Mb
  (466 kb), chr5:50.6–51.1 Mb (456 kb), chr17:21.8–22.3 Mb (456 kb), chr19:25.3–25.7 Mb
  (440 kb), chr2:91.0–91.4 Mb (438 kb), chr10:42.1–42.5 Mb (405 kb).
- Some sinks are telomeric/subtelomeric repeats (e.g. `sink_0` = chr1:0–5525,
  TTAGGG-repeat); the Cen/Sat annotation would mask these.

## Refinement still needed (before the index build)

1. **T2T Cen/Sat v2.1 annotation** (download) for precise satellite masking —
   the heuristic leaves some pericentromeric/telomeric satellite fragments.
2. **Segdup-flanking filter** (Aganezov collapsed-loci / GIAB segdup BEDs) to keep
   genuine missing-copy candidates and drop assembly-gap sequence.
3. **Alt-mask near-identical orthologs** so a sink does not steal reads from a
   real GRCh38 copy (the FixItFelix FANCD2 lesson).

## Step 2 plan (index build + map + validate) — disk-heavy

1. Append `sinks.fa` as **ordinary contigs** to GRCh38 → augmented fasta (~3.1 Gb).
2. Build the CERTA index on the augmented reference (~20.5 Gb; sinks add only
   ~0.07 Gb). **Disk:** galaxy `/` ~41 Gb free ⇒ augmented fasta + index ≈ 23.6 Gb
   leaves ~17 Gb — feasible but tight; build in a scratch dir with `df` guards and
   remove intermediates. (A full second index is the cost; B2's payoff is one-pass
   mapping vs `--t2t`'s two indexes.)
3. Map reads; mates placed on `sink_*` contigs → **MAPQ 0 + sink-coordinate tag**
   (abstain in GRCh38 coordinates; the certificate is over GRCh38 ∪ sinks). This
   is identifiable by contig name in post-processing (no aligner logic change
   needed for the metric).
4. **Per-locus validation** vs HG002 truth: every sink must rescue missing-copy
   reads without stealing real GRCh38 reads; freeze HG002 / hold out HG005.

## Go / no-go

- Step 1 produces a compact, plausible sink set (11.2 Mb, 269 contigs). Proceed to
  annotation-based refinement, then the step-2 index build (disk-guarded) to test
  whether sinks rescue the missing-copy confident errors (the B1 ceiling's ~39%).

## Disk safety & confirmation

Step 1 reads hs1.fa into RAM; outputs `sinks.fa` (11.3 Mb) + `sinks.bed` (8.7 KB)
in `/dev/shm` (RAM). No index built; no disk pressure. Confirmed galaxy ARM64.

## Reproduce

```bash
python3 bench/b2_build.py ~/_paper_work/hs1.fa chm13_gaps.bed sinks   # galaxy
```
