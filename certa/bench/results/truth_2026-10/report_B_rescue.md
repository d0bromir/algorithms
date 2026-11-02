Mates: 39,887,166; with a GRCh38 truth 31,931,012; no linear truth 10,899; other (excluded or not labelled) 7,945,255.

**Placement (mates with a truth).** Misplaced: in the SAM file, not at the truth.

| Aligner | in SAM % | at truth % | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|---|---|
| certa --bio (rescue) | 99.867 | 99.790 | 24,737 | 24,628 | 9,000 | 4,609 | 3,163 |
| certa --bio | 99.689 | 99.614 | 24,146 | 24,040 | 8,845 | 4,621 | 3,173 |
| certa | 99.947 | 99.861 | 27,210 | 13,775 | 9,590 | 4,923 | 3,345 |
| minibwa | 100.000 | 99.901 | 31,498 | 17,723 | 12,663 | 8,578 | 3,374 |
| bwa-mem2 | 100.000 | 99.900 | 32,061 | 17,812 | 15,046 | 9,858 | 6,017 |
| strobealign | 100.000 | 99.901 | 31,711 | 18,505 | 8,646 | 6,331 | 5,503 |
| minimap2 | 99.999 | 99.899 | 31,966 | 18,294 | 11,369 | 8,236 | 5,979 |

**Avoidable errors**: misplaced mates that some other aligner placed at the truth (the rest no aligner placed there), by MAPQ.

| Aligner | misplaced | MAPQ ≥ 1 | ≥ 20 | ≥ 40 | ≥ 60 |
|---|---|---|---|---|---|
| certa --bio (rescue) | 12,143 | 12,100 | 439 | 104 | 31 |
| certa --bio | 12,023 | 11,980 | 398 | 103 | 31 |
| certa | 12,137 | 2,317 | 1,010 | 495 | 345 |
| minibwa | 12,851 | 2,033 | 961 | 699 | 500 |
| bwa-mem2 | 13,398 | 2,362 | 1,570 | 1,100 | 835 |
| strobealign | 13,038 | 2,284 | 723 | 552 | 515 |
| minimap2 | 13,335 | 2,421 | 1,142 | 899 | 615 |

**Placed against misplaced as the MAPQ threshold falls** (mates with a truth; placed / misplaced).

| Aligner | MAPQ ≥ 60 | MAPQ ≥ 50 | MAPQ ≥ 40 | MAPQ ≥ 30 | MAPQ ≥ 20 | MAPQ ≥ 10 | MAPQ ≥ 1 | MAPQ ≥ 0 |
|---|---|---|---|---|---|---|---|---|
| certa --bio (rescue) | 28,132,233 / 3,163 | 29,476,098 / 3,825 | 30,518,882 / 4,609 | 31,129,440 / 7,767 | 31,466,225 / 9,000 | 31,726,621 / 10,547 | 31,888,441 / 24,628 | 31,888,560 / 24,737 |
| certa --bio | 28,150,348 / 3,173 | 29,487,936 / 3,836 | 30,526,106 / 4,621 | 31,116,175 / 7,721 | 31,434,925 / 8,845 | 31,679,287 / 10,224 | 31,831,639 / 24,040 | 31,831,753 / 24,146 |
| certa | 31,596,726 / 3,345 | 31,674,465 / 3,999 | 31,751,509 / 4,923 | 31,852,703 / 8,782 | 31,858,645 / 9,590 | 31,880,808 / 13,349 | 31,883,289 / 13,775 | 31,913,929 / 27,210 |
| minibwa | 31,501,529 / 3,374 | 31,858,820 / 6,524 | 31,879,502 / 8,578 | 31,885,054 / 9,556 | 31,894,420 / 12,663 | 31,899,492 / 14,745 | 31,904,989 / 17,723 | 31,930,981 / 31,498 |
| bwa-mem2 | 31,828,486 / 6,017 | 31,855,567 / 6,899 | 31,886,149 / 9,858 | 31,890,239 / 10,837 | 31,900,277 / 15,046 | 31,902,216 / 16,123 | 31,905,008 / 17,812 | 31,930,999 / 32,061 |
| strobealign | 31,791,993 / 5,503 | 31,809,588 / 5,781 | 31,831,409 / 6,331 | 31,850,700 / 7,191 | 31,861,800 / 8,646 | 31,887,808 / 14,343 | 31,904,652 / 18,505 | 31,931,008 / 31,711 |
| minimap2 | 31,722,644 / 5,979 | 31,808,458 / 6,966 | 31,868,976 / 8,236 | 31,880,413 / 9,380 | 31,889,614 / 11,369 | 31,896,142 / 13,349 | 31,905,559 / 18,294 | 31,930,824 / 31,966 |

**Strata** (mates with a truth in each; at truth % / misplaced with MAPQ ≥ 20).

| Stratum | mates | certa --bio (rescue) | certa --bio | certa | minibwa | bwa-mem2 | strobealign | minimap2 |
|---|---|---|---|---|---|---|---|---|
| bench | 31,272,338 | 99.875 / 2,526 | 99.704 / 2,431 | 99.932 / 3,071 | 99.969 / 3,877 | 99.968 / 4,448 | 99.968 / 2,637 | 99.968 / 3,517 |
| easy | 20,612,201 | 99.944 / 47 | 99.841 / 38 | 99.987 / 145 | 99.999 / 285 | 99.998 / 368 | 99.999 / 280 | 99.998 / 301 |
| difficult | 11,304,482 | 99.546 / 8,309 | 99.237 / 8,166 | 99.663 / 8,777 | 99.755 / 11,231 | 99.751 / 13,080 | 99.753 / 7,425 | 99.750 / 10,096 |
| lowmap_segdup | 1,558,316 | 97.877 / 7,772 | 97.470 / 7,706 | 98.104 / 7,696 | 98.340 / 9,529 | 98.325 / 11,118 | 98.320 / 5,789 | 98.312 / 8,382 |
| segdup | 857,838 | 96.944 / 6,743 | 96.655 / 6,713 | 97.228 / 6,662 | 97.426 / 8,392 | 97.393 / 9,802 | 97.405 / 4,984 | 97.393 / 7,290 |
| lowmap | 1,169,028 | 97.417 / 7,279 | 96.968 / 7,218 | 97.654 / 7,052 | 97.901 / 8,677 | 97.883 / 10,146 | 97.877 / 5,023 | 97.867 / 7,571 |
| tandem | 8,287,406 | 99.633 / 3,691 | 99.278 / 3,567 | 99.755 / 4,198 | 99.854 / 5,539 | 99.850 / 6,390 | 99.854 / 3,951 | 99.850 / 5,069 |
| mhc | 54,001 | 97.689 / 16 | 96.939 / 15 | 99.004 / 33 | 99.813 / 72 | 99.806 / 86 | 99.793 / 57 | 99.769 / 62 |
| cmrg_falsedup | 1,115 | 56.682 / 62 | 56.682 / 62 | 24.215 / 48 | 61.704 / 62 | 61.256 / 64 | 57.668 / 20 | 59.193 / 32 |
| par | 1,192 | 86.829 / 0 | 82.299 / 0 | 89.765 / 0 | 97.651 / 0 | 97.315 / 0 | 97.819 / 1 | 97.819 / 0 |

**Calibration** (mates with a truth): misplaced / expected from MAPQ (Σ 10^(−MAPQ/10)).

| Aligner | MAPQ 0–9 | MAPQ 10–19 | MAPQ 20–29 | MAPQ 30–39 | MAPQ 40–49 | MAPQ 50–60 |
|---|---|---|---|---|---|---|
| certa --bio (rescue) | 161,939: 14,190 / 48,100.3 | 260,396: 1,547 / 12,683.8 | 336,785: 1,233 / 1,365.2 | 610,558: 3,158 / 223.9 | 1,042,784: 784 / 42.0 | 29,476,098: 3,825 / 32.6 |
| certa --bio | 152,466: 13,922 / 45,663.5 | 244,362: 1,379 / 11,983.8 | 318,750: 1,124 / 1,283.3 | 590,069: 3,100 / 214.8 | 1,038,170: 785 / 41.7 | 29,487,936: 3,836 / 32.6 |
| certa | 33,121: 13,861 / 31,301.3 | 22,163: 3,759 / 789.6 | 5,942: 808 / 21.8 | 101,194: 3,859 / 27.5 | 77,044: 924 / 6.8 | 31,674,465: 3,999 / 31.9 |
| minibwa | 31,489: 16,753 / 28,894.5 | 5,072: 2,082 / 202.1 | 9,366: 3,107 / 62.3 | 5,552: 978 / 2.7 | 20,682: 2,054 / 0.5 | 31,858,820: 6,524 / 33.0 |
| bwa-mem2 | 28,783: 15,938 / 26,906.7 | 1,939: 1,077 / 73.1 | 10,038: 4,209 / 31.9 | 4,090: 979 / 1.6 | 30,582: 2,959 / 1.4 | 31,855,567: 6,899 / 31.9 |
| strobealign | 43,200: 17,368 / 32,689.8 | 26,008: 5,697 / 1,608.0 | 11,100: 1,455 / 66.7 | 19,291: 860 / 13.8 | 21,821: 550 / 1.6 | 31,809,588: 5,781 / 31.9 |
| minimap2 | 34,682: 18,617 / 31,143.4 | 6,528: 1,980 / 259.5 | 9,201: 1,989 / 40.6 | 11,437: 1,144 / 4.6 | 60,518: 1,270 / 1.4 | 31,808,458: 6,966 / 32.1 |

**certa --bio (rescue) against certa** (mates with a truth).

| | mates |
|---|---|
| at truth in both | 31,854,789 |
| moved to the truth (certa wrong or absent) | 9,034 |
| moved away (certa --bio (rescue) wrong or absent) | 31,930 |
| (certa --bio (rescue) wrong, certa at truth) | 6,407 |
| wrong in both | 18,059 |
| certa --bio (rescue) wrong with MAPQ ≥ 20: of them YA > 0 (tandem offset) | 9,000: 25 |
| certa --bio (rescue) wrong with MAPQ ≥ 40: of them YA > 0 (tandem offset) | 4,609: 17 |
| placed only by the pair (YR:i:1) | 56,807 |
| of them at truth; misplaced with MAPQ ≥ 20 / 40 | 56,188; 187 / 3 |
| of them certa at truth / wrong / absent | 55,196 / 807 / 804 |

**No linear truth** (10,899 mates: origin not 1:1 in GRCh38): placed, and with MAPQ ≥ 20.

| Aligner | placed | MAPQ ≥ 20 |
|---|---|---|
| certa --bio (rescue) | 32.535 % | 1,996 |
| certa --bio | 28.874 % | 1,829 |
| certa | 43.288 % | 3,215 |
| minibwa | 96.935 % | 5,190 |
| bwa-mem2 | 98.587 % | 5,945 |
| strobealign | 99.624 % | 5,689 |
| minimap2 | 94.899 % | 4,584 |

**Soundness**: other aligners' end-to-end alignments scored under the certa --bio (rescue) model against its certified optimum p1 (38,222,668 certified mates; for a mate placed only by the pair, its YH, which bounds every alignment of it). Alignments with an insertion next to a deletion are outside the model (THEORY §12) and counted apart.

| Aligner | scored | on certified mates | cheaper than p1 (violations) | equal to p1 | above p1 | I next to D (outside the model) |
|---|---|---|---|---|---|---|
| certa | 37,979,843 | 37,760,004 | 0 | 37,339,842 | 420,162 | 72 |
| minibwa | 38,835,199 | 38,068,970 | 0 | 37,398,876 | 670,094 | 381 |
| bwa-mem2 | 38,795,068 | 38,048,172 | 0 | 37,410,757 | 637,415 | 0 |
| strobealign | 38,757,573 | 37,993,538 | 0 | 37,303,454 | 690,084 | 4,237 |
| minimap2 | 38,771,299 | 38,047,973 | 0 | 37,337,690 | 710,283 | 299 |

**Why certa --bio (rescue) misplaces**: for its misplaced, certified mates, the cheapest alignment at the truth that another aligner reported, against certa --bio (rescue)'s optimum p1 (model penalty, Phred).

| MAPQ of the misplaced mate | no other aligner at the truth | equal to p1 (a tie) | 1–9 above p1 | 10–29 above | 30–59 above | 60 or more above |
|---|---|---|---|---|---|---|
| 0–19 | 3,835 | 10,269 | 241 | 169 | 510 | 281 |
| 20–39 | 3,960 | 14 | 7 | 26 | 121 | 79 |
| 40–60 | 4,533 | 10 | 5 | 8 | 14 | 36 |
