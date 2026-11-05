## Test B (real NovaSeq X): fitted on the even pairs of tests A and B, scored on the odd ones

At the locus (150 bp). cv_pair, cv_t2t: depth 4, 300 trees; _d6: depth 6, 400 trees. "the run's MAPQ": interval mode, --iv-mapq cert --iv-pair all.

| run | 99.75 % | 99.44 % | 95.58 % | MAPQ >= 40 |
|---|---|---|---|---|
| the run's MAPQ | 220 | 121 | 75 | 99.44 % |
| --mapq pair (cv_pair.h) | 118 | 71 | 18 | 92.92 % |
| --mapq pair (cv_t2t.h) | 72 | 33 | 10 | 94.74 % |
| certa --mapq cal | 172 | 103 | 34 | 74.51 % |
| minibwa | 137 | 93 | 48 | 99.84 % |
| bwa-mem2 | 143 | 119 | 119 | 99.86 % |
| strobealign | 169 | 115 | 115 | 99.69 % |
| minimap2 | 168 | 134 | 125 | 99.81 % |
| --mapq pair (cv_pair_d6.h) | 115 | 70 | 18 | 94.62 % |
| --mapq pair (cv_t2t_d6.h) | 63 | 32 | 9 | 96.24 % |

15,967,388 mates with a truth; misplaced per million among each run's N most confident.
