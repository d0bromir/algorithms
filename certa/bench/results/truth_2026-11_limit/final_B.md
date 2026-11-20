# Test B: the default build (joint seeds, occ 4096, budget 8192) against the build with them off

## Matched sensitivity and moved mates

| run | 99.75 % | 99.44 % | 95.58 % |
|---|---|---|---|
| votes off | 170.9 (5,444) | 134.9 (4,282) | 81.0 (2,473) |
| final | 170.8 (5,440) | 134.7 (4,277) | 80.9 (2,469) |

31,931,012 mates with a truth; misplaced per million among each run's N most confident (count).

mates placed differently: 318,375; with a truth 1,110: at the locus before 101, after 256 (wrong to right 255, right to wrong 100); unplaced before 8, after 0
MAPQ >= 20, with a truth: 31,872,738 -> 31,872,715 mates; misplaced 9,053 -> 9,026
MAPQ >= 40, with a truth: 31,764,759 -> 31,764,762 mates; misplaced 4,405 -> 4,401
MAPQ >= 40 status changes (mates with a truth):
  9 becomes confident, moved, at the locus
  6 stops being confident, moved, elsewhere
  6 becomes confident, moved, elsewhere
  4 stops being confident, moved, at the locus
  1 stops being confident, same place, elsewhere
  1 stops being confident, same place, at the locus

## The oracle's classes

Mates with a truth: 31,931,012. Pairs scored: 24,535 (a mate away from its locus in some run).

### minibwa

At the locus: 31,902,040 (99.909 %); misplaced 28,941, unplaced 31.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 382 | 12 | 10 | 26 | 214 | 118 | 14 |
| tie | 9,435 | 295 | 45 | 44 | 273 | 9,070 | 3 |
| reference (the truth scores less) | 19,155 | 600 | 6,340 | 3,922 | 4,353 | 4,526 | 14 |

Reference class, how much less the truth scores: 1–4 540, 5–9 3,481, 10–16 2,797, 17–29 7,231, ≥ 30 5,106.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 99.75 % | 131.7 | 0.1 | 0.7 | 130.9 |
| 99.44 % | 88.7 | 0.1 | 0.7 | 88.0 |
| 95.58 % | 45.8 | 0.0 | 0.1 | 45.7 |

### bwa-mem2

At the locus: 31,901,790 (99.908 %); misplaced 29,209, unplaced 13.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 287 | 9 | 9 | 37 | 115 | 119 | 7 |
| tie | 9,623 | 301 | 53 | 72 | 228 | 9,269 | 1 |
| reference (the truth scores less) | 19,312 | 605 | 7,292 | 4,958 | 2,286 | 4,771 | 5 |

Reference class, how much less the truth scores: 1–4 552, 5–9 3,510, 10–16 2,781, 17–29 7,284, ≥ 30 5,185.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 99.75 % | 137.6 | 0.1 | 0.9 | 136.6 |
| 99.44 % | 114.3 | 0.0 | 0.8 | 113.5 |
| 95.58 % | 114.3 | 0.0 | 0.8 | 113.5 |

### certa-off

At the locus: 31,902,932 (99.912 %); misplaced 27,974, unplaced 106.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 312 | 10 | 28 | 16 | 81 | 142 | 45 |
| tie | 8,617 | 270 | 60 | 52 | 167 | 8,324 | 14 |
| reference (the truth scores less) | 19,151 | 600 | 4,317 | 4,580 | 4,010 | 6,197 | 47 |

Reference class, how much less the truth scores: 1–4 678, 5–9 3,607, 10–16 2,736, 17–29 7,075, ≥ 30 5,055.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 99.75 % | 170.9 | 0.9 | 2.1 | 167.9 |
| 99.44 % | 134.9 | 0.9 | 1.8 | 132.2 |
| 95.58 % | 81.0 | 0.6 | 0.9 | 79.5 |

### certa-fin

At the locus: 31,903,087 (99.913 %); misplaced 27,827, unplaced 98.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 147 | 5 | 26 | 10 | 24 | 45 | 42 |
| tie | 8,557 | 268 | 60 | 50 | 152 | 8,281 | 14 |
| reference (the truth scores less) | 19,221 | 602 | 4,315 | 4,565 | 3,766 | 6,533 | 42 |

Reference class, how much less the truth scores: 1–4 623, 5–9 3,569, 10–16 2,729, 17–29 7,157, ≥ 30 5,143.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 99.75 % | 170.8 | 0.8 | 2.1 | 167.8 |
| 99.44 % | 134.7 | 0.8 | 1.8 | 132.1 |
| 95.58 % | 80.9 | 0.6 | 0.9 | 79.4 |

