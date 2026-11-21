| run | 99.75 % | 99.29 % | 98.00 % | 95.00 % | MAPQ >= 40 |
|---|---|---|---|---|---|
| baseline MAPQ | 184.9 | 137.7 | 105.9 | 105.9 | 99.49 % |
| min(MAPQ, gap MAPQ) | 187.0 | 117.5 | 51.2 | 31.6 | 77.89 % |
| gap MAPQ alone | 195.3 | 133.5 | 77.1 | 54.4 | 77.89 % |

31,913,929 mates with a truth (all pairs); misplaced per million among each run's N most confident.

**At each run's own MAPQ >= 40**: confident %, misplaced per million confident.

| run | confident % | misplaced per M |
|---|---|---|
| baseline MAPQ | 99.49 | 155.0 |
| min(MAPQ, gap MAPQ) | 77.89 | 11.3 |
| gap MAPQ alone | 77.89 | 15.1 |

Confident (MAPQ >= 40) misplaced: 4,923; of them competitor gap < 5: 614; gap < 50: 3,811; no competitor score: 0.
  N=99.75%: boundary MAPQ 14, tie group 16,511
  N=99.29%: boundary MAPQ 19, tie group 40,325
  N=98.00%: boundary MAPQ 28, tie group 116,998
  N=95.00%: boundary MAPQ 33, tie group 770,316
