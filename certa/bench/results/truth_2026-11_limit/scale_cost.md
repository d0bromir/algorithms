# The ranking cost of --mapq-scale bwa (bench/results/truth_2026-11_limit)

The default build's --mapq pair runs (with and without --t2t), placed by XP (truth_t2t.py place --xp); the MAPQ written is kBwaScale[XP]. Misplaced per million among the N most confident mates with a truth, at the locus; MAPQ 60: the share of those mates written with MAPQ 60.

## A
| run | 94.13 % | 83.10 % | MAPQ 60 |
|---|---|---|---|
| fpair, XP | 388.8 | 73.2 | 0.00 % |
| fpair, MAPQ (BWA scale) | 391.7 | 75.8 | 67.17 % |
| fpairt2t, XP | 93.1 | 20.5 | 0.00 % |
| fpairt2t, MAPQ (BWA scale) | 94.8 | 22.3 | 85.64 % |

## B
| run | 99.75 % | 99.44 % | 95.58 % | MAPQ 60 |
|---|---|---|---|---|
| fpair, XP | 110.2 | 66.6 | 17.4 | 0.00 % |
| fpair, MAPQ (BWA scale) | 110.2 | 67.3 | 17.5 | 92.62 % |
| fpairt2t, XP | 56.5 | 30.1 | 8.4 | 0.00 % |
| fpairt2t, MAPQ (BWA scale) | 56.5 | 31.4 | 8.9 | 96.29 % |

## H
| run | 99.75 % | 99.29 % | 98.00 % | 95.00 % | MAPQ 60 |
|---|---|---|---|---|---|
| fpair, XP | 207.3 | 61.7 | 20.7 | 9.4 | 0.00 % |
| fpair, MAPQ (BWA scale) | 207.3 | 61.7 | 21.2 | 9.4 | 89.29 % |
| fpairt2t, XP | 99.2 | 25.2 | 7.2 | 2.6 | 0.00 % |
| fpairt2t, MAPQ (BWA scale) | 99.2 | 26.1 | 7.4 | 3.0 | 95.53 % |

