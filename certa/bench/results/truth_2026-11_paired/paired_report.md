Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| paired-bio | 20.077 | 20.061 | 5,004 | 4,986 | 1,770 | 909 | 627 |
| paired-phantom | 20.077 | 20.061 | 5,004 | 4,986 | 1,614 | 356 | 242 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| paired-bio | 0 | 0 | 0 | 0 | 0 |
| paired-phantom | 0 | 0 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| paired-bio | 5,643,748 / 627 | 5,922,837 / 763 | 6,134,256 / 909 | 6,258,217 / 1,521 | 6,326,352 / 1,770 | 6,378,409 / 2,084 | 6,410,819 / 4,986 | 6,410,839 / 5,004 |
| paired-phantom | 5,509,277 / 242 | 5,776,399 / 304 | 5,980,836 / 356 | 6,099,802 / 586 | 6,321,915 / 1,614 | 6,377,727 / 2,054 | 6,410,819 / 4,986 | 6,410,839 / 5,004 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | paired-bio | paired-phantom |
|---|---|---|---|
| bench | 31,272,338 | 20.080 / 482 | 20.080 / 455 |
| easy | 20,612,201 | 20.083 / 13 | 20.083 / 13 |
| difficult | 11,304,482 | 20.030 / 1,606 | 20.030 / 1,470 |
| lowmap_segdup | 1,558,316 | 19.648 / 1,487 | 19.648 / 1,351 |
| segdup | 857,838 | 19.441 / 1,275 | 19.441 / 1,141 |
| lowmap | 1,169,028 | 19.536 / 1,386 | 19.536 / 1,262 |
| tandem | 8,287,406 | 20.061 / 714 | 20.061 / 662 |
| mhc | 54,001 | 19.581 / 0 | 19.581 / 0 |
| cmrg_falsedup | 1,115 | 10.852 / 16 | 10.852 / 15 |
| par | 1,192 | 17.785 / 0 | 17.785 / 0 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| paired-bio | 32,430: 2,920 / 9,574.7 | 52,057: 314 / 2,545.9 | 68,135: 249 / 278.2 | 123,961: 612 / 45.6 | 211,419: 146 / 8.6 | 5,922,837: 763 / 6.6 |
| paired-phantom | 33,112: 2,950 / 9,699.6 | 55,812: 440 / 2,627.6 | 222,113: 1,028 / 600.2 | 118,966: 230 / 43.6 | 204,437: 52 / 8.3 | 5,776,399: 304 / 6.4 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| paired-bio | 6.698 % | 429 |
| paired-phantom | 6.698 % | 384 |
