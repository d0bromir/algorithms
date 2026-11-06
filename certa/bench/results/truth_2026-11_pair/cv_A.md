## Test A (simulated): fitted on the even pairs of tests A and B, scored on the odd ones

At the locus (150 bp). cv_pair, cv_t2t: depth 4, 300 trees; _d6: depth 6, 400 trees. "the run's MAPQ": interval mode, --iv-mapq cert --iv-pair all.

| run | 94.13 % | 83.10 % | MAPQ >= 40 |
|---|---|---|---|
| the run's MAPQ | 1,166 | 708 | 94.13 % |
| --mapq pair (cv_pair.h) | 405 | 70 | 68.22 % |
| --mapq pair (cv_t2t.h) | 133 | 23 | 84.02 % |
| certa --mapq cal | 1,069 | 267 | 36.49 % |
| minibwa | 433 | 262 | 96.63 % |
| bwa-mem2 | 778 | 778 | 96.88 % |
| strobealign | 1,005 | 1,005 | 95.90 % |
| minimap2 | 844 | 844 | 96.22 % |
| --mapq pair (cv_pair_d6.h) | 370 | 67 | 74.44 % |
| --mapq pair (cv_t2t_d6.h) | 111 | 23 | 85.49 % |

9,409,729 mates with a truth; misplaced per million among each run's N most confident.
