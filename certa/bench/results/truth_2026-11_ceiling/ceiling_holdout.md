| run | 99.75 % | 99.29 % | 98.00 % | 95.00 % | MAPQ >= 40 |
|---|---|---|---|---|---|
| baseline MAPQ | 190.4 | 141.0 | 108.1 | 108.1 | 99.49 % |
| min(MAPQ, gap MAPQ) | 191.2 | 120.2 | 52.0 | 31.7 | 77.89 % |
| gap MAPQ alone | 199.1 | 136.9 | 78.9 | 55.4 | 77.89 % |

15,958,866 mates with a truth (odd pairs held out); misplaced per million among each run's N most confident.

**At each run's own MAPQ >= 40**: confident %, misplaced per million confident.

| run | confident % | misplaced per M |
|---|---|---|
| baseline MAPQ | 99.49 | 159.7 |
| min(MAPQ, gap MAPQ) | 77.89 | 11.2 |
| gap MAPQ alone | 77.89 | 15.2 |

Confident (MAPQ >= 40) misplaced: 2,535; of them competitor gap < 5: 334; gap < 50: 1,971; no competitor score: 0.
  N=99.75%: boundary MAPQ 15, tie group 8,239
  N=99.29%: boundary MAPQ 20, tie group 16,468
  N=98.00%: boundary MAPQ 27, tie group 58,447
  N=95.00%: boundary MAPQ 33, tie group 266,518
