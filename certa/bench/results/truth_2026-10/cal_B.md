## Test B: BWA-MEM2 on T2T-CHM13 as the second opinion

Mates with a truth: 31,931,012. T2T-CHM13 records: 39,815,804 mapped, 37,593,189 on GRCh38; with a truth, 31,791,464 (99.563 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.44 %) | N = --bio + minibwa's (95.58 %) | N = minibwa's (99.84 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.44 → 99.40 | 122 → 65 | 77 → 41 | 346 → 319 | 13,650 | 14.1 % |
| certa + minibwa (MAPQ / 4) | 99.44 → 99.40 | 123 → 65 | 77 → 41 | 295 → 227 | 13,629 | 14.0 % |
| certa --bio + minibwa (MAPQ / 4) | 95.58 → 95.53 | 332 → 137 | 138 → 72 | 402 → 175 | 13,984 | 14.4 % |
| minibwa | 99.84 → 99.79 | 89 → 52 | 46 → 27 | 201 → 172 | 17,066 | 19.6 % |
| bwa-mem2 | 99.86 → 99.80 | 114 → 63 | 114 → 63 | 195 → 166 | 17,944 | 22.5 % |
| strobealign | 99.69 → 99.64 | 111 → 60 | 111 → 60 | 282 → 195 | 14,855 | 13.7 % |
| minimap2 | 99.81 → 99.75 | 130 → 71 | 121 → 64 | 219 → 175 | 16,590 | 18.7 % |
| certa interval mode, --iv-pair all --mapq cal | 74.51 → 74.50 | 100 → 58 | 34 → 21 | 281 → 195 | 3,721 | 4.6 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 122 (99.44 %) | 63 (99.39 %) | 124 (99.84 %) | 87 (99.09 %) | 118 (99.54 %) | 65 (99.40 %) | 125 (99.84 %) | 93 (99.10 %) | 119 (99.54 %) | 72 (99.40 %) | 129 (99.77 %) | 112 (99.10 %) | 123 (99.47 %) |
| minibwa | 89 (99.84 %) | 51 (99.78 %) | 139 (99.85 %) | 62 (99.48 %) | 136 (99.55 %) | 52 (99.79 %) | 139 (99.85 %) | 64 (99.48 %) | 136 (99.55 %) | 56 (99.79 %) | 140 (99.85 %) | 71 (99.49 %) | 140 (99.55 %) |
| bwa-mem2 | 114 (99.86 %) | 62 (99.80 %) | 146 (99.85 %) | 64 (99.50 %) | 137 (99.55 %) | 63 (99.80 %) | 148 (99.85 %) | 65 (99.50 %) | 139 (99.55 %) | 68 (99.81 %) | 152 (99.85 %) | 74 (99.51 %) | 144 (99.55 %) |
| strobealign | 111 (99.69 %) | 58 (99.64 %) | 152 (99.85 %) | 80 (99.34 %) | 140 (99.54 %) | 60 (99.64 %) | 153 (99.85 %) | 84 (99.34 %) | 142 (99.55 %) | 65 (99.65 %) | 159 (99.82 %) | 93 (99.35 %) | 148 (99.53 %) |

## Test B: CERTA on T2T-CHM13 as the second opinion

Mates with a truth: 31,931,012. T2T-CHM13 records: 39,704,967 mapped, 37,496,581 on GRCh38; with a truth, 31,791,057 (99.562 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.44 %) | N = --bio + minibwa's (95.58 %) | N = minibwa's (99.84 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.44 → 99.40 | 122 → 78 | 77 → 48 | 346 → 323 | 11,364 | 13.3 % |
| certa + minibwa (MAPQ / 4) | 99.44 → 99.40 | 123 → 79 | 77 → 48 | 295 → 224 | 11,343 | 13.2 % |
| certa --bio + minibwa (MAPQ / 4) | 95.58 → 95.54 | 332 → 187 | 138 → 89 | 402 → 226 | 10,956 | 13.6 % |
| minibwa | 99.84 → 99.80 | 89 → 59 | 46 → 29 | 201 → 185 | 13,368 | 18.9 % |
| bwa-mem2 | 99.86 → 99.82 | 114 → 76 | 114 → 76 | 195 → 177 | 14,081 | 21.9 % |
| strobealign | 99.69 → 99.65 | 111 → 70 | 111 → 70 | 282 → 242 | 11,783 | 13.4 % |
| minimap2 | 99.81 → 99.76 | 130 → 86 | 121 → 78 | 219 → 184 | 13,058 | 18.1 % |
| certa interval mode, --iv-pair all --mapq cal | 74.51 → 74.50 | 100 → 70 | 34 → 23 | 281 → 231 | 2,986 | 5.6 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 122 (99.44 %) | 67 (99.40 %) | 112 (99.54 %) | 98 (99.10 %) | 130 (99.25 %) | 78 (99.40 %) | 119 (99.55 %) | 129 (99.11 %) | 142 (99.25 %) | 88 (99.41 %) | 125 (99.51 %) | 145 (99.12 %) | 147 (99.22 %) |
| minibwa | 89 (99.84 %) | 53 (99.79 %) | 116 (99.85 %) | 66 (99.49 %) | 132 (99.55 %) | 59 (99.80 %) | 118 (99.86 %) | 79 (99.50 %) | 143 (99.56 %) | 65 (99.80 %) | 123 (99.86 %) | 88 (99.51 %) | 145 (99.56 %) |
| bwa-mem2 | 114 (99.86 %) | 65 (99.81 %) | 135 (99.85 %) | 68 (99.50 %) | 126 (99.55 %) | 76 (99.82 %) | 145 (99.86 %) | 83 (99.52 %) | 138 (99.56 %) | 83 (99.82 %) | 153 (99.86 %) | 91 (99.53 %) | 146 (99.57 %) |
| strobealign | 111 (99.69 %) | 61 (99.65 %) | 133 (99.80 %) | 85 (99.35 %) | 130 (99.50 %) | 70 (99.65 %) | 142 (99.80 %) | 102 (99.35 %) | 139 (99.51 %) | 77 (99.66 %) | 149 (99.78 %) | 110 (99.36 %) | 144 (99.49 %) |
