## Test A, window 20 bp

Mates with a truth: 18,818,766. Placed at the truth by no GRCh38 run: 301,659 (1.603 %).

**Matched confident mates**: misplaced / avoidable per million among each run's N most confident mates with a truth, N = the confident mates (MAPQ ≥ 40, or 60) of the run named.

| run | at truth % | N = certa interval mode, --iv-mapq cert ≥ 40 (94.13 %) | N = certa --bio + minibwa (MAPQ / 4) ≥ 40 (83.10 %) | N = minibwa ≥ 40 (96.63 %) | N = minibwa ≥ 60 (92.95 %) |
|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 97.690 | 1,316 / 112.7 | 818 / 63.8 | 4,930 / 668.0 | 1,201 / 104.4 |
| certa + minibwa (MAPQ / 4) | 97.699 | 1,349 / 127.7 | 841 / 71.8 | 4,729 / 654.5 | 1,232 / 119.1 |
| certa --bio + minibwa (MAPQ / 4) | 97.630 | 4,349 / 234.5 | 1,363 / 30.4 | 5,347 / 448.8 | 4,078 / 191.3 |
| minibwa | 97.777 | 607 / 42.7 | 421 / 38.5 | 2,081 / 80.0 | 421 / 38.5 |
| bwa-mem2 | 97.791 | 981 / 74.8 | 981 / 74.8 | 2,332 / 122.8 | 981 / 74.8 |
| strobealign | 97.727 | 1,182 / 90.7 | 1,182 / 90.7 | 3,387 / 400.1 | 1,182 / 90.7 |
| minimap2 | 97.747 | 1,043 / 59.0 | 1,043 / 59.0 | 2,430 / 128.2 | 1,043 / 59.0 |

**At each run's own MAPQ ≥ 40**: confident %, misplaced per million confident, and misplaced in easy regions at MAPQ ≥ 20.

| run | confident % | misplaced per M | easy, MAPQ ≥ 20 |
|---|---|---|---|
| certa interval mode, --iv-mapq cert | 94.13 | 1,316 | 148 |
| certa + minibwa (MAPQ / 4) | 94.13 | 1,344 | 148 |
| certa --bio + minibwa (MAPQ / 4) | 83.10 | 1,363 | 55 |
| minibwa | 96.63 | 2,081 | 345 |
| bwa-mem2 | 96.87 | 2,566 | 440 |
| strobealign | 95.90 | 1,504 | 348 |
| minimap2 | 96.20 | 1,805 | 334 |

## Test A, window 150 bp (the locus)

Mates with a truth: 18,818,766. Placed at the truth by no GRCh38 run: 296,586 (1.576 %).

**Matched confident mates**: misplaced / avoidable per million among each run's N most confident mates with a truth, N = the confident mates (MAPQ ≥ 40, or 60) of the run named.

| run | at truth % | N = certa interval mode, --iv-mapq cert ≥ 40 (94.13 %) | N = certa --bio + minibwa (MAPQ / 4) ≥ 40 (83.10 %) | N = minibwa ≥ 40 (96.63 %) | N = minibwa ≥ 60 (92.95 %) |
|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 97.728 | 1,232 / 95.5 | 762 / 57.4 | 4,749 / 622.1 | 1,133 / 94.3 |
| certa + minibwa (MAPQ / 4) | 97.732 | 1,266 / 111.0 | 785 / 65.6 | 4,519 / 606.4 | 1,164 / 109.3 |
| certa --bio + minibwa (MAPQ / 4) | 97.662 | 4,245 / 207.4 | 1,309 / 13.2 | 5,086 / 392.5 | 3,983 / 165.8 |
| minibwa | 97.811 | 430 / 3.3 | 262 / 2.2 | 1,877 / 34.5 | 262 / 2.2 |
| bwa-mem2 | 97.828 | 767 / 3.4 | 767 / 3.4 | 2,095 / 46.4 | 767 / 3.4 |
| strobealign | 97.759 | 995 / 52.5 | 995 / 52.5 | 3,174 / 360.9 | 995 / 52.5 |
| minimap2 | 97.782 | 831 / 13.1 | 831 / 13.1 | 2,186 / 74.6 | 831 / 13.1 |

**At each run's own MAPQ ≥ 40**: confident %, misplaced per million confident, and misplaced in easy regions at MAPQ ≥ 20.

| run | confident % | misplaced per M | easy, MAPQ ≥ 20 |
|---|---|---|---|
| certa interval mode, --iv-mapq cert | 94.13 | 1,232 | 49 |
| certa + minibwa (MAPQ / 4) | 94.13 | 1,261 | 49 |
| certa --bio + minibwa (MAPQ / 4) | 83.10 | 1,309 | 35 |
| minibwa | 96.63 | 1,877 | 73 |
| bwa-mem2 | 96.87 | 2,326 | 77 |
| strobealign | 95.90 | 1,311 | 93 |
| minimap2 | 96.20 | 1,571 | 63 |

## Test A, near misses (truth_nearmiss.py)

**ivcert**: 1,490 confident near misses:
- soft clip ≥ 10: 49; offset less the clip, quartiles 13, 18, 29
- indel ≥ 10: 25; offset, quartiles 28, 31, 44
- neither: 1,416; offset, quartiles 27, 34, 49

**bio**: 840 confident near misses:
- soft clip ≥ 10: 0; offset less the clip, quartiles –
- indel ≥ 10: 2; offset, quartiles 25, 26, 27
- neither: 838; offset, quartiles 31, 43.5, 78.75

**minibwa**: 3,715 confident near misses:
- soft clip ≥ 10: 1,803; offset less the clip, quartiles -10, 2, 13
- indel ≥ 10: 226; offset, quartiles 28, 33, 48
- neither: 1,686; offset, quartiles 25, 31, 40

**bwa-mem2**: 4,364 confident near misses:
- soft clip ≥ 10: 2,340; offset less the clip, quartiles -13, -1, 11
- indel ≥ 10: 156; offset, quartiles 26, 30, 42.5
- neither: 1,868; offset, quartiles 25, 31, 41

**strobealign**: 3,477 confident near misses:
- soft clip ≥ 10: 1,803; offset less the clip, quartiles -9, 3, 14
- indel ≥ 10: 206; offset, quartiles 26, 32, 45
- neither: 1,468; offset, quartiles 25, 30, 39

**minimap2**: 4,247 confident near misses:
- soft clip ≥ 10: 1,892; offset less the clip, quartiles -11, 2, 13
- indel ≥ 10: 265; offset, quartiles 28, 34, 49
- neither: 2,090; offset, quartiles 26, 33, 48
