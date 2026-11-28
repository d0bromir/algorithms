Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| grch38 | 20.104 | 20.085 | 6,115 | 5,838 | 2,412 | 1,630 | 206 |
| refined | 20.104 | 20.084 | 6,115 | 5,838 | 2,411 | 1,629 | 206 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| grch38 | 0 | 0 | 0 | 0 | 0 |
| refined | 1 | 1 | 0 | 0 | 0 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| grch38 | 5,975,668 / 206 | 6,403,487 / 1,211 | 6,409,241 / 1,630 | 6,410,506 / 1,837 | 6,412,112 / 2,412 | 6,415,772 / 4,069 | 6,419,052 / 5,838 | 6,419,371 / 6,115 |
| refined | 5,974,591 / 206 | 6,402,389 / 1,210 | 6,408,162 / 1,629 | 6,409,429 / 1,836 | 6,411,107 / 2,411 | 6,415,577 / 4,067 | 6,418,933 / 5,838 | 6,419,253 / 6,115 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | grch38 | refined |
|---|---|---|---|
| bench | 31,272,338 | 20.098 / 851 | 20.097 / 851 |
| easy | 20,612,201 | 20.094 / 59 | 20.094 / 59 |
| difficult | 11,304,482 | 20.074 / 2,158 | 20.074 / 2,157 |
| lowmap_segdup | 1,558,316 | 19.761 / 1,761 | 19.760 / 1,761 |
| segdup | 857,838 | 19.559 / 1,495 | 19.559 / 1,495 |
| lowmap | 1,169,028 | 19.659 / 1,608 | 19.659 / 1,608 |
| tandem | 8,287,406 | 20.106 / 1,108 | 20.106 / 1,107 |
| mhc | 54,001 | 20.111 / 6 | 20.092 / 6 |
| cmrg_falsedup | 1,115 | 5.022 / 19 | 5.022 / 19 |
| par | 1,192 | 20.050 / 1 | 20.050 / 1 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| grch38 | 3,599: 2,046 / 1,923.9 | 3,660: 1,657 / 135.7 | 1,606: 575 / 6.4 | 1,265: 207 / 0.5 | 5,754: 419 / 0.1 | 6,403,487: 1,211 / 6.8 |
| refined | 3,676: 2,048 / 1,958.8 | 4,470: 1,656 / 156.8 | 1,678: 575 / 6.6 | 1,267: 207 / 0.5 | 5,773: 419 / 0.1 | 6,402,389: 1,210 / 6.8 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| grch38 | 18.415 % | 901 |
| refined | 18.295 % | 897 |
