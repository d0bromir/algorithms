**minibwa's MAPQ on the fallback mates** (with a truth): mates, misplaced, expected from MAPQ, and the Phred of the misplaced fraction.

| MAPQ | `certa` + minibwa | `certa --bio` + minibwa |
|---|---|---|
| 0–0 | 60,522: 58,141 / 60,522.0, Phred 0.2 | 18,974: 17,096 / 18,974.0, Phred 0.5 |
| 1–9 | 19,785: 14,032 / 11,818.3, Phred 1.5 | 10,353: 8,038 / 6,283.8, Phred 1.1 |
| 10–19 | 11,524: 3,496 / 540.4, Phred 5.2 | 6,495: 4,004 / 272.8, Phred 2.1 |
| 20–29 | 5,773: 1,113 / 25.9, Phred 7.1 | 3,237: 1,592 / 15.5, Phred 3.1 |
| 30–39 | 5,998: 828 / 2.7, Phred 8.6 | 2,802: 1,176 / 1.2, Phred 3.8 |
| 40–49 | 21,442: 1,364 / 0.6, Phred 12.0 | 9,332: 2,121 / 0.2, Phred 6.4 |
| 50–59 | 23,030: 1,148 / 0.1, Phred 13.0 | 14,839: 1,780 / 0.1, Phred 9.2 |
| 60–60 | 62,088: 1,340 / 0.1, Phred 16.7 | 68,225: 2,178 / 0.1, Phred 15.0 |

**The same, divided by 4** (MAPQ / 4, rounding down).

| MAPQ / 4 | `certa` + minibwa | `certa --bio` + minibwa |
|---|---|---|
| 0–0 | 74,167: 68,837 / 74,167.0, 0.9× | 26,242: 22,984 / 26,242.0, 0.9× |
| 1–4 | 17,664: 6,832 / 10,274.2, 0.7× | 9,580: 6,154 / 5,428.2, 1.1× |
| 5–9 | 11,771: 1,941 / 2,477.6, 0.8× | 6,039: 2,768 / 1,311.9, 2.1× |
| 10–14 | 44,472: 2,512 / 2,631.7, 1.0× | 24,171: 3,901 / 1,360.5, 2.9× |
| 15–19 | 62,088: 1,340 / 1,963.4, 0.7× | 68,225: 2,178 / 2,157.5, 1.0× |

**Interval mode by tier** (mates with a truth): placed, misplaced; with MAPQ ≥ 40; with XQ (`--iv-mapq cons`) ≥ 20.

| tier | placed: misplaced | MAPQ ≥ 40 | XQ ≥ 20 |
|---|---|---|---|
| S0 | 7,196,395: 16,204 | 7,112,282: 4,482 | 6,864,268: 3,768 |
| S1 | 9,831,457: 52,162 | 9,701,271: 9,021 | 213: 8 |
| SR | 36,266: 34,834 | 33: 2 | 0: 0 |
| PR | 1,521,393: 225,985 | 900,737: 9,800 | 960,953: 30,192 |
| R | 75,581: 24,541 | 39,566: 461 | 5,706: 53 |
| IV | 59,698: 16,980 | 44,157: 4,773 | 0: 0 |
| MR | 97,168: 63,222 | 20,107: 1,240 | 0: 0 |

`--iv-mapq cert` against the runs with the default MAPQ and `--iv-mapq cons`: the same placements, and XQ for tiers R, IV and MR, the point estimate otherwise: yes.


**MAPQ rules** (mates with a truth). Confident: MAPQ ≥ 40. Avoidable: misplaced while another run of the report placed the mate at the truth.

| run | at truth % | confident % | misplaced / avoidable per M, MAPQ ≥ 40 | ≥ 60 | easy regions, misplaced at MAPQ ≥ 20 |
|---|---|---|---|---|---|
| minibwa | 97.777 | 96.63 | 2,081 / 80.0 | 421 / 38.5 | 345 |
| bwa-mem2 | 97.791 | 96.87 | 2,566 / 138.9 | 981 / 74.8 | 440 |
| strobealign | 97.727 | 95.90 | 1,504 / 121.0 | 1,182 / 90.7 | 348 |
| minimap2 | 97.747 | 96.20 | 1,805 / 96.3 | 1,043 / 59.0 | 334 |
| `certa + minibwa` | 97.699 | 94.70 | 1,552 / 151.8 | 915 / 84.6 | 297 |
| `certa + minibwa`, fallback MAPQ ≤ 39 | 97.699 | 94.13 | 1,344 / 126.5 | 841 / 71.8 | 297 |
| `certa + minibwa`, fallback MAPQ / 4 | 97.699 | 94.13 | 1,344 / 126.5 | 841 / 71.8 | 148 |
| `certa --bio + minibwa` | 97.630 | 83.59 | 1,741 / 69.2 | 1,020 / 48.2 | 307 |
| `certa --bio + minibwa`, fallback MAPQ ≤ 39 | 97.630 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 307 |
| `certa --bio + minibwa`, fallback MAPQ / 4 | 97.630 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 55 |
| `certa interval mode` | 97.690 | 94.68 | 1,671 / 233.8 | 1,055 / 159.9 | 497 |
| `certa interval mode`, `--iv-mapq cons` | 97.690 | 3.70 | 8,671 / 378.9 | 5,710 / 47.2 | 88 |
| `certa interval mode`, XQ for IV, MR | 97.690 | 94.34 | 1,339 / 115.1 | 834 / 65.3 | 153 |
| `certa interval mode`, `--iv-mapq cert` | 97.690 | 94.13 | 1,316 / 112.7 | 818 / 63.8 | 148 |
| `certa --bio` + interval mode (composed) | 97.611 | 83.49 | 1,779 / 164.7 | 1,222 / 165.3 | 462 |
| `certa --bio` + interval mode, XQ | 97.611 | 83.20 | 1,373 / 32.3 | 837 / 21.3 | 110 |
| `certa --bio` + interval mode, MAPQ / 4 | 97.611 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 55 |
| *`certa --bio` alone* | 97.118 | 83.10 | 1,363 / 30.4 | 837 / 21.3 | 55 |
| *`certa` alone* | 97.015 | 94.13 | 1,344 / 126.5 | 841 / 71.8 | 148 |
