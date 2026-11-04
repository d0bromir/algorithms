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

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 122 (99.44 %) | 63 (99.39 %) | 124 (99.84 %) | 87 (99.09 %) | 118 (99.54 %) | 65 (99.40 %) | 125 (99.84 %) | 93 (99.10 %) | 119 (99.54 %) | 72 (99.40 %) | 129 (99.77 %) | 112 (99.10 %) | 123 (99.47 %) |
| minibwa | 89 (99.84 %) | 51 (99.78 %) | 139 (99.85 %) | 62 (99.48 %) | 136 (99.55 %) | 52 (99.79 %) | 139 (99.85 %) | 64 (99.48 %) | 136 (99.55 %) | 56 (99.79 %) | 140 (99.85 %) | 71 (99.49 %) | 140 (99.55 %) |
| bwa-mem2 | 114 (99.86 %) | 62 (99.80 %) | 146 (99.85 %) | 64 (99.50 %) | 137 (99.55 %) | 63 (99.80 %) | 148 (99.85 %) | 65 (99.50 %) | 139 (99.55 %) | 68 (99.81 %) | 152 (99.85 %) | 74 (99.51 %) | 144 (99.55 %) |
| strobealign | 111 (99.69 %) | 58 (99.64 %) | 152 (99.85 %) | 80 (99.34 %) | 140 (99.54 %) | 60 (99.64 %) | 153 (99.85 %) | 84 (99.34 %) | 142 (99.55 %) | 65 (99.65 %) | 159 (99.82 %) | 93 (99.35 %) | 148 (99.53 %) |

## Test B, control: BWA-MEM2 on GRCh38 as the second opinion

Mates with a truth: 31,931,012. bwa-mem2 records: 39,807,452 mapped, 39,807,452 on GRCh38; with a truth, 31,901,790 (99.908 %) at the locus.

**The rule** (bwa-mem2 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (99.44 %) | N = --bio + minibwa's (95.58 %) | N = minibwa's (99.84 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.44 → 99.44 | 122 → 121 | 77 → 76 | 346 → 333 | 61 | 90.2 % |
| certa + minibwa (MAPQ / 4) | 99.44 → 99.44 | 123 → 121 | 77 → 75 | 295 → 288 | 60 | 91.7 % |
| certa --bio + minibwa (MAPQ / 4) | 95.58 → 95.58 | 332 → 325 | 138 → 138 | 402 → 391 | 11 | 45.5 % |
| minibwa | 99.84 → 99.84 | 89 → 89 | 46 → 46 | 201 → 201 | 10 | 50.0 % |
| bwa-mem2 | 99.86 → 99.86 | 114 → 114 | 114 → 114 | 195 → 195 | 0 | 0.0 % |
| strobealign | 99.69 → 99.69 | 111 → 110 | 111 → 110 | 282 → 280 | 74 | 62.2 % |
| minimap2 | 99.81 → 99.81 | 130 → 128 | 121 → 119 | 219 → 215 | 169 | 83.4 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 122 (99.44 %) | 119 (99.44 %) | 126 (99.86 %) | 119 (99.44 %) | 126 (99.86 %) | 121 (99.44 %) | 127 (99.87 %) | 121 (99.44 %) | 127 (99.87 %) | 122 (99.44 %) | 129 (99.75 %) | 122 (99.44 %) | 129 (99.75 %) |
| minibwa | 89 (99.84 %) | 88 (99.84 %) | 116 (99.86 %) | 88 (99.84 %) | 116 (99.86 %) | 89 (99.84 %) | 116 (99.86 %) | 89 (99.84 %) | 116 (99.86 %) | 89 (99.84 %) | 116 (99.85 %) | 89 (99.84 %) | 116 (99.85 %) |
| bwa-mem2 | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) | 114 (99.86 %) |
| strobealign | 111 (99.69 %) | 109 (99.69 %) | 143 (99.87 %) | 109 (99.69 %) | 143 (99.87 %) | 110 (99.69 %) | 144 (99.87 %) | 110 (99.69 %) | 144 (99.87 %) | 110 (99.69 %) | 145 (99.81 %) | 110 (99.69 %) | 145 (99.81 %) |
