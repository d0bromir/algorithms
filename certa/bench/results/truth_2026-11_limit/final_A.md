# Test A: the default build (joint seeds, occ 4096, budget 8192) against the build with them off

## Matched sensitivity and moved mates

| run | 94.13 % | 83.10 % |
|---|---|---|
| votes off | 1,265.4 (22,416) | 753.9 (11,790) |
| final | 1,262.3 (22,361) | 752.4 (11,766) |

18,818,766 mates with a truth; misplaced per million among each run's N most confident (count).

mates placed differently: 184,597; with a truth 42,539: at the locus before 792, after 8,569 (wrong to right 8,556, right to wrong 779); unplaced before 17, after 0
MAPQ >= 20, with a truth: 18,081,983 -> 18,081,707 mates; misplaced 57,503 -> 57,233
MAPQ >= 40, with a truth: 17,825,160 -> 17,825,120 mates; misplaced 23,954 -> 23,894
MAPQ >= 40 status changes (mates with a truth):
  47 stops being confident, moved, elsewhere
  44 becomes confident, moved, at the locus
  23 stops being confident, moved, at the locus
  20 stops being confident, same place, elsewhere
  17 becomes confident, moved, elsewhere
  13 stops being confident, same place, at the locus
  2 becomes confident, same place, at the locus

## The oracle's classes

Mates with a truth: 18,818,766. Pairs scored: 250,820 (a mate away from its locus in some run).

### minibwa

At the locus: 18,406,849 (97.811 %); misplaced 411,891, unplaced 26.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 9,668 | 514 | 53 | 135 | 3,183 | 6,288 | 9 |
| tie | 139,864 | 7,432 | 155 | 185 | 4,808 | 134,714 | 2 |
| reference (the truth scores less) | 262,385 | 13,943 | 33,916 | 29,344 | 55,474 | 143,636 | 15 |

Reference class, how much less the truth scores: 1–4 4,858, 5–9 77,254, 10–16 73,874, 17–29 60,717, ≥ 30 45,682.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 94.13 % | 429.1 | 0.9 | 3.0 | 425.1 |
| 83.10 % | 261.7 | 0.6 | 1.4 | 259.6 |

### bwa-mem2

At the locus: 18,410,041 (97.828 %); misplaced 408,714, unplaced 11.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 2,797 | 149 | 47 | 207 | 806 | 1,735 | 2 |
| tie | 138,979 | 7,385 | 186 | 540 | 2,996 | 135,255 | 2 |
| reference (the truth scores less) | 266,949 | 14,185 | 42,179 | 44,599 | 28,216 | 151,948 | 7 |

Reference class, how much less the truth scores: 1–4 5,184, 5–9 77,659, 10–16 74,903, 17–29 62,030, ≥ 30 47,173.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 94.13 % | 766.6 | 0.2 | 3.2 | 763.1 |
| 83.10 % | 766.6 | 0.2 | 3.2 | 763.1 |

### certa-off

At the locus: 18,405,625 (97.805 %); misplaced 412,981, unplaced 160.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 21,699 | 1,153 | 361 | 354 | 2,085 | 18,838 | 61 |
| tie | 134,586 | 7,152 | 669 | 499 | 3,138 | 130,255 | 25 |
| reference (the truth scores less) | 256,856 | 13,649 | 22,924 | 32,696 | 33,070 | 168,092 | 74 |

Reference class, how much less the truth scores: 1–4 6,615, 5–9 76,746, 10–16 73,180, 17–29 57,623, ≥ 30 42,692.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 94.13 % | 1,265.4 | 19.6 | 35.4 | 1,210.5 |
| 83.10 % | 753.9 | 11.3 | 23.6 | 719.0 |

### certa-fin

At the locus: 18,413,402 (97.846 %); misplaced 405,221, unplaced 143.

| class | mates | per million | MAPQ ≥ 40 | 20–39 | 1–19 | 0 | unplaced |
|---|---|---|---|---|---|---|---|
| missed (the truth scores more) | 6,165 | 328 | 322 | 249 | 429 | 5,118 | 47 |
| tie | 135,378 | 7,194 | 665 | 484 | 2,776 | 131,428 | 25 |
| reference (the truth scores less) | 263,821 | 14,019 | 22,907 | 32,606 | 31,261 | 176,976 | 71 |

Reference class, how much less the truth scores: 1–4 5,997, 5–9 77,770, 10–16 74,759, 17–29 60,222, ≥ 30 45,073.

| share of the mates with a truth | misplaced per million | missed | tie | reference |
|---|---|---|---|---|
| 94.13 % | 1,262.3 | 17.5 | 35.1 | 1,209.7 |
| 83.10 % | 752.4 | 10.2 | 23.5 | 718.7 |

