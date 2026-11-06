## HG005: BWA-MEM2 on T2T-CHM13 as the second opinion

Mates with a truth: 30,127,192. T2T-CHM13 records: 39,808,591 mapped, 37,612,665 on GRCh38; with a truth, 29,990,855 (99.547 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.29 %) | N = minibwa's (99.75 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.29 → 99.24 | 111 → 58 | 363 → 311 | 12,913 | 12.8 % |
| minibwa | 99.75 → 99.70 | 80 → 47 | 209 → 168 | 16,649 | 19.6 % |
| bwa-mem2 | 99.77 → 99.71 | 108 → 60 | 203 → 162 | 17,473 | 22.6 % |
| strobealign | 99.57 → 99.52 | 108 → 56 | 324 → 207 | 14,618 | 14.0 % |
| minimap2 | 99.67 → 99.62 | 152 → 75 | 253 → 176 | 16,012 | 18.4 % |
| certa interval mode, --iv-pair all --mapq cal | 72.78 → 72.77 | 86 → 50 | 320 → 185 | 3,612 | 2.4 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 111 (99.29 %) | 55 (99.24 %) | 113 (99.76 %) | 69 (98.93 %) | 105 (99.45 %) | 58 (99.24 %) | 115 (99.76 %) | 76 (98.94 %) | 107 (99.45 %) | 67 (99.25 %) | 119 (99.68 %) | 98 (98.95 %) | 112 (99.37 %) |
| minibwa | 80 (99.75 %) | 46 (99.69 %) | 111 (99.77 %) | 50 (99.38 %) | 103 (99.46 %) | 47 (99.70 %) | 112 (99.78 %) | 53 (99.38 %) | 103 (99.46 %) | 52 (99.70 %) | 113 (99.77 %) | 62 (99.39 %) | 106 (99.46 %) |
| bwa-mem2 | 108 (99.77 %) | 58 (99.71 %) | 124 (99.77 %) | 50 (99.39 %) | 110 (99.46 %) | 60 (99.71 %) | 126 (99.78 %) | 53 (99.40 %) | 112 (99.46 %) | 66 (99.72 %) | 132 (99.77 %) | 63 (99.41 %) | 119 (99.46 %) |
| strobealign | 108 (99.57 %) | 54 (99.52 %) | 130 (99.77 %) | 61 (99.21 %) | 113 (99.45 %) | 56 (99.52 %) | 132 (99.77 %) | 64 (99.21 %) | 115 (99.46 %) | 63 (99.53 %) | 139 (99.73 %) | 77 (99.22 %) | 124 (99.43 %) |

## HG005: CERTA on T2T-CHM13 as the second opinion, certified placements only

Mates with a truth: 30,127,192. T2T-CHM13 records: 39,556,477 mapped, 37,378,846 on GRCh38; with a truth, 29,990,031 (99.545 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.29 %) | N = minibwa's (99.75 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.29 → 99.26 | 111 → 72 | 363 → 319 | 9,414 | 13.2 % |
| minibwa | 99.75 → 99.71 | 80 → 54 | 209 → 175 | 11,650 | 20.7 % |
| bwa-mem2 | 99.77 → 99.73 | 108 → 72 | 203 → 169 | 12,295 | 24.0 % |
| strobealign | 99.57 → 99.54 | 108 → 68 | 324 → 251 | 10,222 | 15.0 % |
| minimap2 | 99.67 → 99.63 | 152 → 93 | 253 → 199 | 11,293 | 19.9 % |
| certa interval mode, --iv-pair all --mapq cal | 72.78 → 72.78 | 86 → 59 | 320 → 229 | 2,125 | 3.8 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 111 (99.29 %) | 60 (99.25 %) | 99 (99.45 %) | 80 (98.95 %) | 115 (99.14 %) | 72 (99.26 %) | 105 (99.45 %) | 111 (98.95 %) | 128 (99.15 %) | 82 (99.26 %) | 111 (99.41 %) | 117 (98.96 %) | 135 (99.11 %) |
| minibwa | 80 (99.75 %) | 49 (99.70 %) | 92 (99.77 %) | 55 (99.39 %) | 105 (99.46 %) | 54 (99.71 %) | 94 (99.78 %) | 68 (99.40 %) | 114 (99.47 %) | 60 (99.72 %) | 99 (99.78 %) | 77 (99.42 %) | 120 (99.47 %) |
| bwa-mem2 | 108 (99.77 %) | 62 (99.72 %) | 117 (99.78 %) | 55 (99.41 %) | 103 (99.46 %) | 72 (99.73 %) | 127 (99.78 %) | 69 (99.42 %) | 115 (99.48 %) | 80 (99.74 %) | 135 (99.78 %) | 78 (99.43 %) | 125 (99.48 %) |
| strobealign | 108 (99.57 %) | 59 (99.53 %) | 116 (99.72 %) | 67 (99.22 %) | 99 (99.41 %) | 68 (99.54 %) | 125 (99.72 %) | 84 (99.23 %) | 110 (99.42 %) | 75 (99.54 %) | 132 (99.70 %) | 96 (99.24 %) | 120 (99.40 %) |

## HG005, control: BWA-MEM2 on GRCh38 as the second opinion

Mates with a truth: 30,127,192. bwa-mem2 records: 39,799,556 mapped, 39,799,556 on GRCh38; with a truth, 30,098,782 (99.906 %) at the locus.

**The rule** (bwa-mem2 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.29 %) | N = minibwa's (99.75 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.29 → 99.29 | 111 → 109 | 363 → 347 | 49 | 91.8 % |
| minibwa | 99.75 → 99.75 | 80 → 79 | 209 → 208 | 37 | 51.4 % |
| bwa-mem2 | 99.77 → 99.77 | 108 → 108 | 203 → 203 | 0 | 0.0 % |
| strobealign | 99.57 → 99.57 | 108 → 105 | 324 → 319 | 165 | 71.5 % |
| minimap2 | 99.67 → 99.67 | 152 → 141 | 253 → 241 | 411 | 90.3 % |
| certa interval mode, --iv-pair all --mapq cal | 72.78 → 72.78 | 86 → 86 | 320 → 316 | 0 | 0.0 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 111 (99.29 %) | 107 (99.29 %) | 119 (99.78 %) | 107 (99.29 %) | 119 (99.78 %) | 109 (99.29 %) | 120 (99.78 %) | 109 (99.29 %) | 120 (99.78 %) | 110 (99.29 %) | 121 (99.65 %) | 110 (99.29 %) | 121 (99.65 %) |
| minibwa | 80 (99.75 %) | 79 (99.75 %) | 110 (99.78 %) | 79 (99.75 %) | 110 (99.78 %) | 79 (99.75 %) | 111 (99.78 %) | 79 (99.75 %) | 111 (99.78 %) | 80 (99.75 %) | 111 (99.77 %) | 80 (99.75 %) | 111 (99.77 %) |
| bwa-mem2 | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) | 108 (99.77 %) |
| strobealign | 108 (99.57 %) | 104 (99.57 %) | 142 (99.78 %) | 104 (99.57 %) | 142 (99.78 %) | 105 (99.57 %) | 143 (99.78 %) | 105 (99.57 %) | 143 (99.78 %) | 106 (99.57 %) | 145 (99.72 %) | 106 (99.57 %) | 145 (99.72 %) |
