## Test A: BWA-MEM2 on T2T-CHM13 as the second opinion

Mates with a truth: 18,818,766. T2T-CHM13 records: 19,996,538 mapped, 18,783,376 on GRCh38; with a truth, 18,308,587 (97.289 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (94.13 %) | N = --bio + minibwa's (83.10 %) | N = minibwa's (96.63 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 94.13 → 94.04 | 1,232 → 1,023 | 762 → 615 | 4,749 → 4,775 | 16,854 | 27.6 % |
| certa + minibwa (MAPQ / 4) | 94.13 → 94.04 | 1,266 → 1,050 | 785 → 634 | 4,519 → 3,721 | 16,949 | 28.1 % |
| certa --bio + minibwa (MAPQ / 4) | 83.10 → 83.01 | 4,245 → 3,536 | 1,309 → 1,056 | 5,086 → 4,252 | 16,973 | 24.0 % |
| minibwa | 96.63 → 96.49 | 430 → 354 | 262 → 209 | 1,877 → 1,650 | 25,773 | 30.0 % |
| bwa-mem2 | 96.87 → 96.72 | 767 → 597 | 767 → 597 | 2,095 → 1,734 | 28,483 | 33.6 % |
| strobealign | 95.90 → 95.80 | 995 → 812 | 995 → 812 | 3,174 → 2,768 | 18,186 | 25.4 % |
| minimap2 | 96.20 → 96.08 | 831 → 648 | 831 → 648 | 2,186 → 1,977 | 22,977 | 27.9 % |
| certa interval mode, --iv-pair all --mapq cal | 36.49 → 36.48 | 1,063 → 894 | 267 → 224 | 4,141 → 3,770 | 1,527 | 6.8 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 1,232 (94.13 %) | 927 (94.02 %) | 778 (96.63 %) | 1,033 (93.67 %) | 725 (96.28 %) | 1,023 (94.04 %) | 831 (96.65 %) | 1,156 (93.72 %) | 795 (96.32 %) | 1,158 (94.07 %) | 902 (95.73 %) | 1,341 (93.77 %) | 883 (95.43 %) |
| minibwa | 430 (96.63 %) | 324 (96.45 %) | 487 (96.84 %) | 344 (96.08 %) | 457 (96.46 %) | 354 (96.49 %) | 503 (96.88 %) | 391 (96.14 %) | 480 (96.53 %) | 384 (96.55 %) | 521 (96.67 %) | 436 (96.24 %) | 505 (96.35 %) |
| bwa-mem2 | 767 (96.87 %) | 535 (96.67 %) | 738 (96.89 %) | 471 (96.29 %) | 675 (96.51 %) | 597 (96.72 %) | 800 (96.94 %) | 550 (96.37 %) | 754 (96.59 %) | 660 (96.79 %) | 863 (96.88 %) | 631 (96.47 %) | 834 (96.56 %) |
| strobealign | 995 (95.90 %) | 744 (95.77 %) | 985 (96.71 %) | 652 (95.41 %) | 895 (96.35 %) | 812 (95.80 %) | 1,053 (96.74 %) | 748 (95.46 %) | 990 (96.40 %) | 899 (95.83 %) | 1,139 (96.19 %) | 865 (95.52 %) | 1,106 (95.88 %) |

## Test A: CERTA on T2T-CHM13 as the second opinion

Mates with a truth: 18,818,766. T2T-CHM13 records: 19,992,735 mapped, 18,787,825 on GRCh38; with a truth, 18,298,596 (97.236 %) at the locus.

**The rule** (T2T-CHM13 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (94.13 %) | N = --bio + minibwa's (83.10 %) | N = minibwa's (96.63 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 94.13 → 94.07 | 1,232 → 1,124 | 762 → 671 | 4,749 → 4,781 | 11,360 | 23.8 % |
| certa + minibwa (MAPQ / 4) | 94.13 → 94.07 | 1,266 → 1,156 | 785 → 693 | 4,519 → 4,226 | 11,360 | 23.9 % |
| certa --bio + minibwa (MAPQ / 4) | 83.10 → 83.05 | 4,245 → 3,922 | 1,309 → 1,186 | 5,086 → 4,717 | 9,702 | 20.6 % |
| minibwa | 96.63 → 96.56 | 430 → 383 | 262 → 227 | 1,877 → 1,783 | 14,380 | 26.8 % |
| bwa-mem2 | 96.87 → 96.79 | 767 → 669 | 767 → 669 | 2,095 → 1,925 | 15,798 | 30.7 % |
| strobealign | 95.90 → 95.84 | 995 → 898 | 995 → 898 | 3,174 → 3,015 | 10,744 | 21.7 % |
| minimap2 | 96.20 → 96.13 | 831 → 730 | 831 → 730 | 2,186 → 2,099 | 13,372 | 24.8 % |
| certa interval mode, --iv-pair all --mapq cal | 36.49 → 36.49 | 1,063 → 966 | 267 → 239 | 4,141 → 3,924 | 1,030 | 7.8 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 1,232 (94.13 %) | 1,040 (94.05 %) | 1,165 (94.61 %) | 1,171 (93.72 %) | 1,130 (94.28 %) | 1,124 (94.07 %) | 1,250 (94.63 %) | 1,288 (93.76 %) | 1,237 (94.32 %) | 1,175 (94.09 %) | 1,254 (94.45 %) | 1,351 (93.79 %) | 1,248 (94.15 %) |
| minibwa | 430 (96.63 %) | 354 (96.51 %) | 440 (96.70 %) | 391 (96.17 %) | 419 (96.36 %) | 383 (96.56 %) | 457 (96.74 %) | 433 (96.24 %) | 441 (96.43 %) | 397 (96.58 %) | 465 (96.69 %) | 452 (96.28 %) | 452 (96.40 %) |
| bwa-mem2 | 767 (96.87 %) | 600 (96.74 %) | 786 (96.88 %) | 554 (96.39 %) | 741 (96.54 %) | 669 (96.79 %) | 854 (96.94 %) | 640 (96.47 %) | 826 (96.62 %) | 699 (96.81 %) | 884 (96.90 %) | 677 (96.51 %) | 863 (96.60 %) |
| strobealign | 995 (95.90 %) | 833 (95.81 %) | 1,035 (96.31 %) | 768 (95.48 %) | 971 (95.98 %) | 898 (95.84 %) | 1,100 (96.34 %) | 857 (95.53 %) | 1,059 (96.02 %) | 926 (95.85 %) | 1,128 (96.16 %) | 895 (95.56 %) | 1,097 (95.86 %) |
