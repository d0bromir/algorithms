Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| t2t | 20.104 | 20.085 | 6,115 | 5,410 | 1,745 | 1,232 | 229 |
| t2t-commit | 20.104 | 20.085 | 6,115 | 5,410 | 1,745 | 1,232 | 229 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| t2t | 0 | 0 | 0 | 0 | 0 |
| t2t-commit | 0 | 0 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| t2t | 6,157,820 / 229 | 6,405,858 / 953 | 6,409,773 / 1,232 | 6,412,113 / 1,542 | 6,413,100 / 1,745 | 6,414,302 / 2,254 | 6,418,559 / 5,410 | 6,419,371 / 6,115 |
| t2t-commit | 6,157,820 / 229 | 6,405,858 / 953 | 6,409,773 / 1,232 | 6,412,113 / 1,542 | 6,413,100 / 1,745 | 6,414,302 / 2,254 | 6,418,559 / 5,410 | 6,419,371 / 6,115 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | t2t | t2t-commit |
|---|---|---|---|
| bench | 31,272,338 | 20.098 / 537 | 20.098 / 537 |
| easy | 20,612,201 | 20.094 / 63 | 20.094 / 63 |
| difficult | 11,304,482 | 20.074 / 1,563 | 20.074 / 1,563 |
| lowmap_segdup | 1,558,316 | 19.761 / 1,167 | 19.761 / 1,167 |
| segdup | 857,838 | 19.559 / 969 | 19.559 / 969 |
| lowmap | 1,169,028 | 19.659 / 1,051 | 19.659 / 1,051 |
| tandem | 8,287,406 | 20.106 / 876 | 20.106 / 876 |
| mhc | 54,001 | 20.111 / 5 | 20.111 / 5 |
| cmrg_falsedup | 1,115 | 5.022 / 0 | 5.022 / 0 |
| par | 1,192 | 20.050 / 2 | 20.050 / 2 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| t2t | 5,069: 3,861 / 3,667.3 | 1,202: 509 / 85.6 | 987: 203 / 3.9 | 2,340: 310 / 0.6 | 3,915: 279 / 0.1 | 6,405,858: 953 / 6.6 |
| t2t-commit | 5,069: 3,861 / 3,667.3 | 1,202: 509 / 85.6 | 987: 203 / 3.9 | 2,340: 310 / 0.6 | 3,915: 279 / 0.1 | 6,405,858: 953 / 6.6 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| t2t | 18.415 % | 808 |
| t2t-commit | 18.415 % | 808 |
