# B2-build step 2 — augmented index mapping: sinks work, but heuristic filtering steals reads

**Date:** 2026-11-22. **Step:** Pillar B, B2 (one-pass sink contigs), build step 2.
**Code:** [`../../b2_index.sh`](../../b2_index.sh), [`../../b2_map.sh`](../../b2_map.sh),
[`../../b2_sink_analysis.py`](../../b2_sink_analysis.py). Results:
[`b2_report.md`](b2_report.md), [`b2_sinks.md`](b2_sinks.md).

## What was done

1. **Augmented index in RAM** (`b2_index.sh`): GRCh38 primary + 269 CHM13 sink
   contigs (11.2 Mb, heuristic-filtered) → `augmented.cidx`, q=22 s=1 (matches
   grch38.s1v3). **463 contigs, 2.96 G entries, 23 Gb**, built in `/dev/shm`
   (RAM) in ~2 min — **root disk unchanged (41 Gb free)**, disk-safe.
2. **Map + score** (`b2_map.sh`): 4 M-pair HG002 subsample, interval `--mapq pair`,
   GRCh38 baseline vs augmented. Sink contigs (`sink_*`) are absent from the
   GRCh38 truth reference, so `truth_eval` counts them as not-placed (abstain) —
   the B2 semantics.
3. **Sink analysis** (`b2_sink_analysis.py`): stolen real reads vs correct
   missing-copy abstentions.

## Result

| run | at truth % | misplaced | MAPQ≥40 err | MAPQ≥60 err |
|---|---|---|---|---|
| GRCh38 | 20.085 | 6,115 | 1,630 | 206 |
| augmented | 20.082 | 6,036 | **1,600** | 202 |

- **MAPQ≥40 errors −30 (−1.8%)**; MAPQ≥20 misplaced: segdup 1,495→1,482 (−13),
  difficult 2,158→2,129 (−29), tandem 1,108→1,079 (−29).
- **Missing-copy class** (no linear truth): placed 18.4%→16.5%, MAPQ≥20 901→877
  (−24) — sinks abstain on some missing-copy reads.
- **Sink placements: 82,698 mates** — but only **245** are correct missing-copy
  abstentions; **916 are stolen real reads** (266 at MAPQ≥40); **81,537** are
  unlabelled (mostly repetitive reads attracted to repeat-rich sinks).
- Calibration: sink competition demotes reads into lower MAPQ bands (0–9:
  3,599→4,822; 10–19: 3,660→5,493).

## Honest finding

The B2 **mechanics work** (RAM index, disk-safe; sinks abstain missing-copy reads
and depress competing GRCh38 MAPQ), but the **heuristic satellite filter leaves
too many promiscuous repeat sinks** (pericentromeric/telomeric): they steal 916
real reads (266 confident) to abstain only 245 missing-copy reads — a poor trade.
The confident-error reduction (−30) is small relative to the sensitivity cost.

**Precise sink filtering is essential** (the FixItFelix FANCD2 lesson, observed
directly): the T2T **Cen/Sat v2.1** annotation to mask satellites/telomeres, a
**segdup-flanking filter** (Aganezov collapsed-loci BED) to keep genuine
missing-copy candidates, and **alt-masking of near-identical orthologs** so a
sink cannot steal reads from a real GRCh38 copy. With specific sinks, the
missing-copy abstentions should rise and the stolen-read count fall toward zero.

## Go / no-go

- B2 is **viable but not yet effective** with heuristic sinks. Proceed to
  annotation-based refinement (Cen/Sat + segdup-flanking + alt-mask), rebuild the
  sink set, and re-test (the RAM index pipeline is ready and disk-safe).
- This quantifies *why* naive sink/decoy addition fails (read stealing) and
  motivates the precise filtering — a thesis-worthy negative/positive pair.

## Disk safety & confirmation

Index + SAMs + npz entirely in `/dev/shm` (RAM); root disk unchanged (41 Gb
free). Confirmed galaxy ARM64; a2 not run (no data/numpy, disk full).

## Reproduce

```bash
bash bench/b2_index.sh   # build augmented index in RAM (galaxy)
bash bench/b2_map.sh     # map GRCh38 vs augmented, score, sink analysis
```
