Mates: 39,968,164; with a GRCh38 truth 30,127,192; no linear truth 8,776; other (excluded or not labelled) 9,832,196.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| hg005-bio | 10.054 | 10.045 | 2,539 | 2,524 | 745 | 351 | 264 |
| hg005-phantom | 10.054 | 10.045 | 2,539 | 2,524 | 634 | 180 | 138 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| hg005-bio | 0 | 0 | 0 | 0 | 0 |
| hg005-phantom | 0 | 0 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| hg005-bio | 2,871,252 / 264 | 2,911,757 / 288 | 2,967,966 / 351 | 2,998,756 / 642 | 3,007,909 / 745 | 3,019,216 / 861 | 3,028,862 / 2,524 | 3,028,878 / 2,539 |
| hg005-phantom | 2,801,040 / 138 | 2,839,441 / 155 | 2,893,730 / 180 | 2,922,489 / 265 | 3,005,214 / 634 | 3,018,810 / 845 | 3,028,862 / 2,524 | 3,028,878 / 2,539 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | hg005-bio | hg005-phantom |
|---|---|---|---|
| bench | 28,035,676 | 10.054 / 151 | 10.054 / 116 |
| easy | 19,575,635 | 10.055 / 5 | 10.055 / 5 |
| difficult | 10,510,710 | 10.031 / 686 | 10.031 / 585 |
| lowmap_segdup | 1,496,120 | 9.832 / 644 | 9.832 / 543 |
| segdup | 815,143 | 9.729 / 545 | 9.729 / 444 |
| lowmap | 1,135,643 | 9.789 / 624 | 9.789 / 535 |
| tandem | 7,626,213 | 10.048 / 303 | 10.048 / 278 |
| mhc | 55,975 | 9.963 / 0 | 9.963 / 0 |
| cmrg_falsedup | 1,252 | 4.792 / 0 | 4.792 / 0 |
| par | 25,039 | 9.709 / 1 | 9.709 / 1 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| hg005-bio | 9,662: 1,678 / 3,756.0 | 11,307: 116 / 487.1 | 9,153: 103 / 36.0 | 30,790: 291 / 17.0 | 56,209: 63 / 1.8 | 2,911,757: 288 / 3.0 |
| hg005-phantom | 10,068: 1,694 / 3,818.4 | 13,596: 211 / 542.7 | 82,725: 369 / 194.4 | 28,759: 85 / 15.9 | 54,289: 25 / 1.8 | 2,839,441: 155 / 2.9 |

**No linear truth** (8,776 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| hg005-bio | 3.111 % | 169 |
| hg005-phantom | 3.111 % | 150 |
