## Replay of the frozen build (frozen_source.sha): CERTA's --mapq pair MAPQ against bench/truth_mapq_pair.py replay

truth_pair.sh check (tests A and B; per test: --mapq pair, then --mapq pair --t2t):
```
== 03:08:16 A: replay
[replay] 19,887,522 mapped mates, 0 differ
[replay] 19,887,522 mapped mates, 0 differ
== 03:37:26 B: replay
[replay] 39,546,158 mapped mates, 0 differ
[replay] 39,546,158 mapped mates, 0 differ
== 04:33:49 check finished
```

truth_pair.sh holdout (HG005):
```
== 03:08:16 H: --mapq pair 
[read_sam] 39968164 primary records
[replay] 39,375,943 mapped mates, 0 differ
== 03:27:15 H: --mapq pair --t2t /home/mpiuser/certa-t2t/hs1.cidx --t2t-chain /home/mpiuser/headroom/hs1ToHg38.over.chain.gz
[read_sam] 39968164 primary records
[replay] 39,375,943 mapped mates, 0 differ
[truth_veto] 214,005 placements set to MAPQ 0
[truth_veto] 202,979 placements set to MAPQ 0
== 04:09:56 holdout finished
```
