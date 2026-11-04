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

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 1,232 (94.13 %) | 927 (94.02 %) | 778 (96.63 %) | 1,033 (93.67 %) | 725 (96.28 %) | 1,023 (94.04 %) | 831 (96.65 %) | 1,156 (93.72 %) | 795 (96.32 %) | 1,158 (94.07 %) | 902 (95.73 %) | 1,341 (93.77 %) | 883 (95.43 %) |
| minibwa | 430 (96.63 %) | 324 (96.45 %) | 487 (96.84 %) | 344 (96.08 %) | 457 (96.46 %) | 354 (96.49 %) | 503 (96.88 %) | 391 (96.14 %) | 480 (96.53 %) | 384 (96.55 %) | 521 (96.67 %) | 436 (96.24 %) | 505 (96.35 %) |
| bwa-mem2 | 767 (96.87 %) | 535 (96.67 %) | 738 (96.89 %) | 471 (96.29 %) | 675 (96.51 %) | 597 (96.72 %) | 800 (96.94 %) | 550 (96.37 %) | 754 (96.59 %) | 660 (96.79 %) | 863 (96.88 %) | 631 (96.47 %) | 834 (96.56 %) |
| strobealign | 995 (95.90 %) | 744 (95.77 %) | 985 (96.71 %) | 652 (95.41 %) | 895 (96.35 %) | 812 (95.80 %) | 1,053 (96.74 %) | 748 (95.46 %) | 990 (96.40 %) | 899 (95.83 %) | 1,139 (96.19 %) | 865 (95.52 %) | 1,106 (95.88 %) |

## Test A, control: BWA-MEM2 on GRCh38 as the second opinion

Mates with a truth: 18,818,766. bwa-mem2 records: 19,990,054 mapped, 19,990,054 on GRCh38; with a truth, 18,410,041 (97.828 %) at the locus.

**The rule** (bwa-mem2 MAPQ ≥ 40, elsewhere on GRCh38: MAPQ 0). Misplaced per million among each run's N most confident mates with a truth, as is → with the rule; confident % (MAPQ ≥ 40); and the confident mates the rule sets to 0: how many, and the share misplaced.

| run | confident % | N = interval mode's (94.13 %) | N = --bio + minibwa's (83.10 %) | N = minibwa's (96.63 %) | confident mates set to 0 | of them misplaced |
|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 94.13 → 94.13 | 1,232 → 1,200 | 762 → 739 | 4,749 → 4,567 | 699 | 94.8 % |
| certa + minibwa (MAPQ / 4) | 94.13 → 94.13 | 1,266 → 1,226 | 785 → 759 | 4,519 → 4,341 | 802 | 94.9 % |
| certa --bio + minibwa (MAPQ / 4) | 83.10 → 83.10 | 4,245 → 4,177 | 1,309 → 1,308 | 5,086 → 4,978 | 41 | 61.0 % |
| minibwa | 96.63 → 96.63 | 430 → 429 | 262 → 261 | 1,877 → 1,876 | 48 | 47.9 % |
| bwa-mem2 | 96.87 → 96.87 | 767 → 767 | 767 → 767 | 2,095 → 2,095 | 0 | 0.0 % |
| strobealign | 95.90 → 95.89 | 995 → 976 | 995 → 976 | 3,174 → 3,054 | 684 | 91.2 % |
| minimap2 | 96.20 → 96.20 | 831 → 828 | 831 → 828 | 2,186 → 2,180 | 170 | 61.8 % |

**Variants**, misplaced per million at N = interval mode's (and confident %): the T2T MAPQ threshold; T2T records that do not lift counted as disagreeing; and a placement T2T confirms raised to T2T's MAPQ.

| run | as is | ≥ 20 | ≥ 20, confirm | ≥ 20, unlifted | ≥ 20, unlifted, confirm | ≥ 40 | ≥ 40, confirm | ≥ 40, unlifted | ≥ 40, unlifted, confirm | ≥ 60 | ≥ 60, confirm | ≥ 60, unlifted | ≥ 60, unlifted, confirm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| certa interval mode, --iv-mapq cert | 1,232 (94.13 %) | 1,180 (94.13 %) | 1,060 (96.89 %) | 1,180 (94.13 %) | 1,060 (96.89 %) | 1,200 (94.13 %) | 1,073 (96.89 %) | 1,200 (94.13 %) | 1,073 (96.89 %) | 1,231 (94.13 %) | 1,095 (95.72 %) | 1,231 (94.13 %) | 1,095 (95.72 %) |
| minibwa | 430 (96.63 %) | 429 (96.63 %) | 781 (96.90 %) | 429 (96.63 %) | 781 (96.90 %) | 429 (96.63 %) | 781 (96.90 %) | 429 (96.63 %) | 781 (96.90 %) | 430 (96.63 %) | 781 (96.65 %) | 430 (96.63 %) | 781 (96.65 %) |
| bwa-mem2 | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) | 767 (96.87 %) |
| strobealign | 995 (95.90 %) | 965 (95.89 %) | 1,262 (96.93 %) | 965 (95.89 %) | 1,262 (96.93 %) | 976 (95.89 %) | 1,272 (96.94 %) | 976 (95.89 %) | 1,272 (96.94 %) | 991 (95.90 %) | 1,287 (96.20 %) | 991 (95.90 %) | 1,287 (96.20 %) |
