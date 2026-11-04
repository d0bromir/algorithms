## Test B, window 20 bp

Mates with a truth: 31,931,012. Placed at the truth by no GRCh38 run: 18,536 (0.058 %).

**The mates no GRCh38 run places at the truth**: placed there by each new run, and with MAPQ ≥ 20 / 40.

| run | at the truth | MAPQ ≥ 20 | MAPQ ≥ 40 | misplaced with MAPQ ≥ 40 |
|---|---|---|---|---|
| BWA-MEM2 on T2T-CHM13, lifted | 7,367 (39.7 %) | 7,181 | 6,984 | 8,082 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 6,055 (32.7 %) | 3,098 | 2,689 | 2,795 |

**Matched confident mates**: misplaced / avoidable per million among each run's N most confident mates with a truth, N = the confident mates (MAPQ ≥ 40, or 60) of the run named.

| run | at truth % | N = certa interval mode, --iv-mapq cert ≥ 40 (99.44 %) | N = certa --bio + minibwa (MAPQ / 4) ≥ 40 (95.58 %) | N = minibwa ≥ 40 (99.84 %) | N = minibwa ≥ 60 (98.65 %) | N = BWA-MEM2 on T2T-CHM13, lifted ≥ 40 (99.61 %) | N = vg giraffe, HPRC v2.1 (HG002 held out) ≥ 40 (99.44 %) |
|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.899 | 155 / 15.4 | 106 / 10.8 | 399 / 60.0 | 106 / 10.8 | 175 / 16.2 | 156 / 15.6 |
| certa + minibwa (MAPQ / 4) | 99.902 | 156 / 15.9 | 106 / 10.9 | 356 / 44.1 | 106 / 10.9 | 175 / 16.6 | 156 / 15.9 |
| certa --bio + minibwa (MAPQ / 4) | 99.900 | 394 / 35.7 | 151 / 3.5 | 470 / 57.8 | 291 / 15.6 | 423 / 43.0 | 394 / 35.8 |
| minibwa | 99.901 | 154 / 17.0 | 107 / 16.1 | 269 / 22.5 | 107 / 16.1 | 173 / 17.2 | 154 / 17.0 |
| bwa-mem2 | 99.900 | 189 / 26.5 | 189 / 26.5 | 273 / 30.6 | 189 / 26.5 | 189 / 26.5 | 189 / 26.5 |
| strobealign | 99.901 | 173 / 16.9 | 173 / 16.9 | 351 / 35.7 | 173 / 16.9 | 181 / 17.2 | 173 / 16.9 |
| minimap2 | 99.899 | 198 / 20.6 | 188 / 19.9 | 290 / 31.4 | 188 / 19.9 | 219 / 22.6 | 198 / 20.6 |
| BWA-MEM2 on T2T-CHM13, lifted | 99.545 | 941 / 715.2 | 941 / 715.2 | 1,273 / 987.2 | 941 / 715.2 | 1,025 / 770.7 | 941 / 715.2 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 99.849 | 125 / 37.0 | 102 / 30.9 | 673 / 415.9 | 102 / 30.9 | 194 / 59.0 | 125 / 37.2 |

**At each run's own MAPQ ≥ 40**: confident %, misplaced per million confident, and misplaced in easy regions at MAPQ ≥ 20.

| run | confident % | misplaced per M | easy, MAPQ ≥ 20 |
|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.44 | 155 | 145 |
| certa + minibwa (MAPQ / 4) | 99.44 | 155 | 145 |
| certa --bio + minibwa (MAPQ / 4) | 95.58 | 151 | 47 |
| minibwa | 99.84 | 269 | 285 |
| bwa-mem2 | 99.86 | 309 | 368 |
| strobealign | 99.69 | 199 | 280 |
| minimap2 | 99.81 | 258 | 301 |
| BWA-MEM2 on T2T-CHM13, lifted | 99.61 | 1,025 | 1,139 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 99.44 | 125 | 198 |

## Test B, window 150 bp (the locus)

Mates with a truth: 31,931,012. Placed at the truth by no GRCh38 run: 16,639 (0.052 %).

**The mates no GRCh38 run places at the truth**: placed there by each new run, and with MAPQ ≥ 20 / 40.

| run | at the truth | MAPQ ≥ 20 | MAPQ ≥ 40 | misplaced with MAPQ ≥ 40 |
|---|---|---|---|---|
| BWA-MEM2 on T2T-CHM13, lifted | 6,372 (38.3 %) | 6,177 | 5,984 | 7,329 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 4,636 (27.9 %) | 1,766 | 1,420 | 2,529 |

**Matched confident mates**: misplaced / avoidable per million among each run's N most confident mates with a truth, N = the confident mates (MAPQ ≥ 40, or 60) of the run named.

| run | at truth % | N = certa interval mode, --iv-mapq cert ≥ 40 (99.44 %) | N = certa --bio + minibwa (MAPQ / 4) ≥ 40 (95.58 %) | N = minibwa ≥ 40 (99.84 %) | N = minibwa ≥ 60 (98.65 %) | N = BWA-MEM2 on T2T-CHM13, lifted ≥ 40 (99.61 %) | N = vg giraffe, HPRC v2.1 (HG002 held out) ≥ 40 (99.44 %) |
|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.908 | 122 / 8.9 | 77 / 5.5 | 346 / 45.9 | 77 / 5.5 | 142 / 9.7 | 123 / 9.0 |
| certa + minibwa (MAPQ / 4) | 99.910 | 123 / 9.5 | 77 / 5.7 | 295 / 28.1 | 77 / 5.7 | 142 / 10.2 | 123 / 9.5 |
| certa --bio + minibwa (MAPQ / 4) | 99.907 | 332 / 24.5 | 138 / 1.3 | 402 / 45.7 | 266 / 10.9 | 358 / 31.4 | 332 / 24.6 |
| minibwa | 99.909 | 89 / 0.6 | 46 / 0.3 | 201 / 5.1 | 46 / 0.3 | 108 / 0.7 | 89 / 0.6 |
| bwa-mem2 | 99.908 | 114 / 0.7 | 114 / 0.7 | 195 / 3.9 | 114 / 0.7 | 114 / 0.7 | 114 / 0.7 |
| strobealign | 99.908 | 111 / 3.6 | 111 / 3.6 | 282 / 21.6 | 111 / 3.6 | 118 / 3.9 | 111 / 3.6 |
| minimap2 | 99.907 | 130 / 4.5 | 121 / 4.0 | 219 / 13.9 | 121 / 4.0 | 150 / 5.9 | 130 / 4.5 |
| BWA-MEM2 on T2T-CHM13, lifted | 99.563 | 782 / 578.7 | 782 / 578.7 | 1,092 / 834.9 | 782 / 578.7 | 859 / 628.4 | 782 / 578.7 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 99.856 | 94 / 14.7 | 73 / 9.5 | 609 / 364.2 | 73 / 9.5 | 159 / 33.2 | 95 / 14.9 |

**At each run's own MAPQ ≥ 40**: confident %, misplaced per million confident, and misplaced in easy regions at MAPQ ≥ 20.

| run | confident % | misplaced per M | easy, MAPQ ≥ 20 |
|---|---|---|---|
| certa interval mode, --iv-mapq cert | 99.44 | 122 | 28 |
| certa + minibwa (MAPQ / 4) | 99.44 | 122 | 28 |
| certa --bio + minibwa (MAPQ / 4) | 95.58 | 138 | 8 |
| minibwa | 99.84 | 201 | 29 |
| bwa-mem2 | 99.86 | 231 | 32 |
| strobealign | 99.69 | 135 | 35 |
| minimap2 | 99.81 | 188 | 29 |
| BWA-MEM2 on T2T-CHM13, lifted | 99.61 | 859 | 698 |
| vg giraffe, HPRC v2.1 (HG002 held out) | 99.44 | 95 | 169 |

## Test B, near misses (truth_nearmiss.py)

**ivcert**: 1,055 confident near misses:
- soft clip ≥ 10: 71; offset less the clip, quartiles 11, 16, 24.5
- indel ≥ 10: 7; offset, quartiles 28.5, 30, 32
- neither: 977; offset, quartiles 24, 30, 39

**bio**: 399 confident near misses:
- soft clip ≥ 10: 0; offset less the clip, quartiles –
- indel ≥ 10: 9; offset, quartiles 26, 30, 34
- neither: 390; offset, quartiles 25, 31, 40

**minibwa**: 2,183 confident near misses:
- soft clip ≥ 10: 992; offset less the clip, quartiles -9, 5, 13
- indel ≥ 10: 78; offset, quartiles 27.25, 31.5, 48.5
- neither: 1,113; offset, quartiles 24, 30, 39

**bwa-mem2**: 2,504 confident near misses:
- soft clip ≥ 10: 1,262; offset less the clip, quartiles -12, 0, 11
- indel ≥ 10: 49; offset, quartiles 28, 32, 41
- neither: 1,193; offset, quartiles 24, 30, 38

**strobealign**: 2,022 confident near misses:
- soft clip ≥ 10: 927; offset less the clip, quartiles -9, 6, 15
- indel ≥ 10: 57; offset, quartiles 26, 30, 39
- neither: 1,038; offset, quartiles 24, 30, 38

**minimap2**: 2,256 confident near misses:
- soft clip ≥ 10: 988; offset less the clip, quartiles -9, 5, 14
- indel ≥ 10: 76; offset, quartiles 28, 36, 60
- neither: 1,192; offset, quartiles 25, 30, 40
