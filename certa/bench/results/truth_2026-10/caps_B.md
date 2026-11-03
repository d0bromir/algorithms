**minibwa's MAPQ on the fallback mates** (with a truth): mates, misplaced, expected from MAPQ, and the Phred of the misplaced fraction.

| MAPQ | `certa` + minibwa | `certa --bio` + minibwa |
|---|---|---|
| 0–0 | 1,252: 1,015 / 1,252.0, Phred 0.9 | 2,105: 1,645 / 2,105.0, Phred 1.1 |
| 1–9 | 1,512: 1,153 / 899.4, Phred 1.2 | 2,142: 1,629 / 1,252.7, Phred 1.2 |
| 10–19 | 881: 553 / 35.2, Phred 2.0 | 1,489: 1,010 / 57.7, Phred 1.7 |
| 20–29 | 447: 238 / 2.0, Phred 2.7 | 748: 427 / 3.7, Phred 2.4 |
| 30–39 | 476: 229 / 0.2, Phred 3.2 | 718: 366 / 0.3, Phred 2.9 |
| 40–49 | 1,164: 318 / 0.0, Phred 5.6 | 2,077: 632 / 0.1, Phred 5.2 |
| 50–59 | 1,773: 244 / 0.0, Phred 8.6 | 3,233: 456 / 0.0, Phred 8.5 |
| 60–60 | 9,547: 449 / 0.0, Phred 13.3 | 29,909: 980 / 0.0, Phred 14.8 |

**The same, divided by 4** (MAPQ / 4, rounding down).

| MAPQ / 4 | `certa` + minibwa | `certa --bio` + minibwa |
|---|---|---|
| 0–0 | 2,273: 1,830 / 2,273.0, 0.8× | 3,513: 2,754 / 3,513.0, 0.8× |
| 1–4 | 1,372: 891 / 788.4, 1.1× | 2,223: 1,530 / 1,254.6, 1.2× |
| 5–9 | 923: 467 / 192.7, 2.4× | 1,466: 793 / 315.0, 2.5× |
| 10–14 | 2,937: 562 / 168.2, 3.3× | 5,310: 1,088 / 302.7, 3.6× |
| 15–19 | 9,547: 449 / 301.9, 1.5× | 29,909: 980 / 945.8, 1.0× |

**Interval mode by tier** (mates with a truth): placed, misplaced; with MAPQ ≥ 40; with XQ (`--iv-mapq cons`) ≥ 20.

| tier | placed: misplaced | MAPQ ≥ 40 | XQ ≥ 20 |
|---|---|---|---|
| S0 | 26,044,561: 3,114 | 25,962,942: 1,617 | 26,030,842: 2,374 |
| S1 | 5,614,298: 5,797 | 5,584,799: 2,204 | 61: 0 |
| SR | 7: 4 | 0: 0 | 0: 0 |
| PR | 254,275: 17,753 | 204,123: 1,115 | 208,241: 3,787 |
| R | 4,360: 401 | 3,752: 73 | 1,089: 4 |
| IV | 9,734: 3,470 | 7,337: 1,535 | 0: 0 |
| MR | 3,485: 1,531 | 1,506: 166 | 0: 0 |

`--iv-mapq cert` against the runs with the default MAPQ and `--iv-mapq cons`: the same placements, and XQ for tiers R, IV and MR, the point estimate otherwise: yes.


**MAPQ rules** (mates with a truth). Confident: MAPQ ≥ 40. Avoidable: misplaced while another run of the report placed the mate at the truth.

| run | at truth % | confident % | misplaced / avoidable per M, MAPQ ≥ 40 | ≥ 60 | easy regions, misplaced at MAPQ ≥ 20 |
|---|---|---|---|---|---|
| minibwa | 99.901 | 99.84 | 269 / 22.5 | 107 / 16.1 | 285 |
| bwa-mem2 | 99.900 | 99.86 | 309 / 35.2 | 189 / 26.5 | 368 |
| strobealign | 99.901 | 99.69 | 199 / 18.1 | 173 / 16.9 | 280 |
| minimap2 | 99.899 | 99.81 | 258 / 29.0 | 188 / 19.9 | 301 |
| `certa + minibwa` | 99.902 | 99.48 | 187 / 20.2 | 120 / 13.4 | 253 |
| `certa + minibwa`, fallback MAPQ ≤ 39 | 99.902 | 99.44 | 155 / 15.7 | 106 / 10.9 | 253 |
| `certa + minibwa`, fallback MAPQ / 4 | 99.902 | 99.44 | 155 / 15.7 | 106 / 10.9 | 145 |
| `certa --bio + minibwa` | 99.900 | 95.69 | 219 / 11.6 | 147 / 6.2 | 244 |
| `certa --bio + minibwa`, fallback MAPQ ≤ 39 | 99.900 | 95.58 | 151 / 3.5 | 112 / 1.1 | 244 |
| `certa --bio + minibwa`, fallback MAPQ / 4 | 99.900 | 95.58 | 151 / 3.5 | 112 / 1.1 | 47 |
| `certa interval mode` | 99.899 | 99.48 | 211 / 40.4 | 149 / 32.0 | 464 |
| `certa interval mode`, `--iv-mapq cons` | 99.899 | 0.56 | 4,706 / 330.2 | 3,978 / 86.5 | 71 |
| `certa interval mode`, XQ for IV, MR | 99.899 | 99.45 | 158 / 15.9 | 107 / 11.1 | 151 |
| `certa interval mode`, `--iv-mapq cert` | 99.899 | 99.44 | 155 / 15.4 | 106 / 10.8 | 145 |
| `certa --bio` + interval mode (composed) | 99.897 | 95.67 | 220 / 30.7 | 167 / 25.8 | 423 |
| `certa --bio` + interval mode, XQ | 99.897 | 95.61 | 154 / 4.1 | 112 / 1.1 | 97 |
| `certa --bio` + interval mode, MAPQ / 4 | 99.897 | 95.58 | 151 / 3.5 | 112 / 1.1 | 47 |
| *`certa --bio` alone* | 99.790 | 95.58 | 151 / 3.5 | 112 / 1.1 | 47 |
| *`certa` alone* | 99.861 | 99.44 | 155 / 15.7 | 106 / 10.9 | 145 |
