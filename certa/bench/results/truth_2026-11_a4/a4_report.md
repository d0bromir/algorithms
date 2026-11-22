Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| bio | 10.026 | 10.017 | 2,654 | 2,640 | 819 | 461 | 358 |
| global | 10.026 | 10.017 | 2,654 | 2,640 | 819 | 460 | 280 |
| locus | 10.026 | 10.017 | 2,654 | 2,639 | 609 | 190 | 144 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| bio | 0 | 0 | 0 | 0 | 0 |
| global | 0 | 0 | 0 | 0 | 0 |
| locus | 0 | 0 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| bio | 3,061,069 / 358 | 3,100,682 / 393 | 3,147,250 / 461 | 3,171,346 / 712 | 3,184,107 / 819 | 3,194,030 / 925 | 3,201,288 / 2,640 | 3,201,304 / 2,654 |
| global | 2,635,545 / 280 | 3,091,678 / 390 | 3,146,508 / 460 | 3,171,252 / 712 | 3,184,107 / 819 | 3,194,021 / 925 | 3,201,288 / 2,640 | 3,201,304 / 2,654 |
| locus | 2,986,674 / 144 | 3,024,069 / 160 | 3,069,030 / 190 | 3,091,244 / 278 | 3,169,949 / 609 | 3,193,080 / 862 | 3,201,284 / 2,639 | 3,201,304 / 2,654 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | bio | global | locus |
|---|---|---|---|---|
| bench | 31,272,338 | 10.027 / 218 | 10.027 / 218 | 10.027 / 174 |
| easy | 20,612,201 | 10.026 / 7 | 10.026 / 7 | 10.026 / 7 |
| difficult | 11,304,482 | 10.006 / 738 | 10.006 / 738 | 10.006 / 545 |
| lowmap_segdup | 1,558,316 | 9.806 / 690 | 9.806 / 690 | 9.806 / 497 |
| segdup | 857,838 | 9.689 / 594 | 9.689 / 594 | 9.689 / 402 |
| lowmap | 1,169,028 | 9.754 / 643 | 9.754 / 643 | 9.754 / 467 |
| tandem | 8,287,406 | 10.025 / 320 | 10.025 / 320 | 10.025 / 247 |
| mhc | 54,001 | 9.744 / 0 | 9.744 / 0 | 9.744 / 0 |
| cmrg_falsedup | 1,115 | 5.740 / 3 | 5.740 / 3 | 5.740 / 3 |
| par | 1,192 | 8.557 / 0 | 8.557 / 0 | 8.557 / 0 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| bio | 7,274: 1,729 / 2,607.9 | 9,923: 106 / 482.5 | 12,761: 107 / 46.6 | 24,096: 251 / 7.5 | 46,568: 68 / 1.7 | 3,100,682: 393 / 3.2 |
| global | 7,283: 1,729 / 2,609.0 | 9,914: 106 / 481.6 | 12,855: 107 / 46.8 | 24,744: 252 / 8.0 | 54,830: 70 / 2.2 | 3,091,678: 390 / 3.5 |
| locus | 8,224: 1,792 / 2,977.1 | 23,131: 253 / 765.3 | 78,705: 331 / 255.8 | 22,214: 88 / 6.6 | 44,961: 30 / 1.6 | 3,024,069: 160 / 3.1 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| bio | 3.110 % | 221 |
| global | 3.110 % | 221 |
| locus | 3.110 % | 182 |
