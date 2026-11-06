## HG005 (held out): the frozen build, run once

Build: frozen_source.sha (the models fitted on HG002 tests A and B; mapq_model_t2t.h of depth 6 by the rule fixed on HG002). 20 M NovaSeq 6000 pairs of SRR14724528, truth from HG005's HPRC release 2 haplotypes (truth_hg005.sh). At the locus (150 bp); each run's N most confident mates with a truth, the boundary MAPQ pro rata. "the run's MAPQ": interval mode, --iv-mapq cert --iv-pair all. The other runs are truth_hg005.sh's placements; the veto is truth_veto.py with BWA-MEM2's placements on hs1.

CERTA's --mapq pair MAPQ equals truth_mapq_pair.py's replay on all 39,375,943 mapped mates, with and without --t2t. Mapping (galaxy, cores 0-63): 43.8 s at 24.6 GiB; with --t2t 76.2 s at 45.9 GiB.

| run | 99.75 % | 99.29 % | 98.00 % | 95.00 % | MAPQ >= 40 |
|---|---|---|---|---|---|
| the run's MAPQ | 328 | 106 | 61 | 61 | 99.29 % |
| certa --mapq cal | 319 | 86 | 33 | 22 | 72.79 % |
| certa --mapq pair | 186 | 60 | 20 | 9 | 89.28 % |
| certa --mapq pair --t2t | 93 | 25 | 7 | 3 | 95.61 % |
| minibwa | 208 | 80 | 37 | 37 | 99.75 % |
| bwa-mem2 | 202 | 108 | 108 | 108 | 99.77 % |
| strobealign | 323 | 108 | 108 | 108 | 99.57 % |
| minimap2 | 252 | 152 | 110 | 110 | 99.67 % |
| minibwa, BWA-MEM2's T2T veto | 167 | 48 | 22 | 22 | 99.70 % |
| bwa-mem2, BWA-MEM2's T2T veto | 160 | 60 | 60 | 60 | 99.71 % |

30,127,192 mates with a truth; misplaced per million among each run's N most confident.
