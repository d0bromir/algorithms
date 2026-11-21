Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| baseline | 99.947 | 99.861 | 27,210 | 13,775 | 9,590 | 4,923 | 3,345 |
| mapq-local | 99.947 | 99.861 | 27,210 | 13,763 | 9,562 | 4,891 | 3,325 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| baseline | 0 | 0 | 0 | 0 | 0 |
| mapq-local | 0 | 0 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| baseline | 31,596,726 / 3,345 | 31,674,465 / 3,999 | 31,751,509 / 4,923 | 31,852,703 / 8,782 | 31,858,645 / 9,590 | 31,880,808 / 13,349 | 31,883,289 / 13,775 | 31,913,929 / 27,210 |
| mapq-local | 31,596,078 / 3,325 | 31,674,274 / 3,979 | 31,751,306 / 4,891 | 31,852,637 / 8,756 | 31,858,604 / 9,562 | 31,880,778 / 13,325 | 31,883,272 / 13,763 | 31,913,929 / 27,210 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | baseline | mapq-local |
|---|---|---|---|
| bench | 31,272,338 | 99.932 / 3,071 | 99.932 / 3,065 |
| easy | 20,612,201 | 99.987 / 145 | 99.987 / 145 |
| difficult | 11,304,482 | 99.663 / 8,777 | 99.663 / 8,753 |
| lowmap_segdup | 1,558,316 | 98.104 / 7,696 | 98.104 / 7,672 |
| segdup | 857,838 | 97.228 / 6,662 | 97.228 / 6,639 |
| lowmap | 1,169,028 | 97.654 / 7,052 | 97.654 / 7,032 |
| tandem | 8,287,406 | 99.755 / 4,198 | 99.755 / 4,185 |
| mhc | 54,001 | 99.004 / 33 | 99.004 / 33 |
| cmrg_falsedup | 1,115 | 24.215 / 48 | 24.215 / 48 |
| par | 1,192 | 89.765 / 0 | 89.765 / 0 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| baseline | 33,121: 13,861 / 31,301.3 | 22,163: 3,759 / 789.6 | 5,942: 808 / 21.8 | 101,194: 3,859 / 27.5 | 77,044: 924 / 6.8 | 31,674,465: 3,999 / 31.9 |
| mapq-local | 33,151: 13,885 / 31,322.1 | 22,174: 3,763 / 790.2 | 5,967: 806 / 21.9 | 101,331: 3,865 / 27.6 | 77,032: 912 / 6.8 | 31,674,274: 3,979 / 31.9 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| baseline | 43.288 % | 3,215 |
| mapq-local | 43.288 % | 3,212 |
