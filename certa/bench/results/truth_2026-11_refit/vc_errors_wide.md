# DeepVariant errors with tier R's wide band, by region

DeepVariant 1.10.0 on the variant-calling subsets (HG002 against v5.0q and CMRG, HG005 against v4.2.1); hap.py PASS
errors per run from `vc/errors_table.py`. CERTA's runs with the refit models, without and with the wide band
(`--iv-wide`); costs in VALIDATION_RESULTS.md, "Tier R's wide band".


v5q, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,522 | 348 | 1,377 | 517 | 3,764 |
| CERTA, default MAPQ, wide band | 1,511 | 344 | 1,357 | 512 | 3,724 |
| CERTA --mapq pair | 1,465 | 365 | 1,333 | 498 | 3,661 |
| CERTA --mapq pair, wide band | 1,458 | 364 | 1,321 | 497 | 3,640 |
| CERTA --mapq pair --t2t | 1,309 | 362 | 1,304 | 502 | 3,477 |
| CERTA --mapq pair --t2t, wide band | 1,301 | 355 | 1,287 | 503 | 3,446 |
| minibwa | 1,436 | 309 | 1,281 | 494 | 3,520 |
| BWA-MEM2 | 1,419 | 342 | 1,329 | 519 | 3,609 |

v5q, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 584 | 83 | 139 | 32 | 838 |
| CERTA, default MAPQ, wide band | 581 | 81 | 139 | 33 | 834 |
| CERTA --mapq pair | 537 | 88 | 135 | 36 | 796 |
| CERTA --mapq pair, wide band | 534 | 86 | 135 | 37 | 792 |
| CERTA --mapq pair --t2t | 515 | 84 | 124 | 30 | 753 |
| CERTA --mapq pair --t2t, wide band | 512 | 81 | 124 | 31 | 748 |
| minibwa | 556 | 71 | 127 | 30 | 784 |
| BWA-MEM2 | 548 | 82 | 130 | 29 | 789 |

v5q, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,448 | 331 | 1,363 | 515 | 3,657 |
| CERTA, default MAPQ, wide band | 1,438 | 329 | 1,344 | 510 | 3,621 |
| CERTA --mapq pair | 1,377 | 352 | 1,318 | 496 | 3,543 |
| CERTA --mapq pair, wide band | 1,371 | 352 | 1,307 | 495 | 3,525 |
| CERTA --mapq pair --t2t | 1,237 | 350 | 1,290 | 500 | 3,377 |
| CERTA --mapq pair --t2t, wide band | 1,230 | 343 | 1,274 | 501 | 3,348 |
| minibwa | 1,388 | 290 | 1,271 | 490 | 3,439 |
| BWA-MEM2 | 1,367 | 322 | 1,318 | 514 | 3,521 |

v5q, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 74 | 17 | 17 | 2 | 110 |
| CERTA, default MAPQ, wide band | 73 | 15 | 16 | 2 | 106 |
| CERTA --mapq pair | 88 | 13 | 17 | 2 | 120 |
| CERTA --mapq pair, wide band | 87 | 12 | 17 | 2 | 118 |
| CERTA --mapq pair --t2t | 72 | 12 | 17 | 2 | 103 |
| CERTA --mapq pair --t2t, wide band | 71 | 12 | 16 | 2 | 101 |
| minibwa | 48 | 19 | 13 | 3 | 83 |
| BWA-MEM2 | 52 | 20 | 13 | 4 | 89 |

cmrg, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 714 | 129 | 239 | 60 | 1,142 |
| CERTA, default MAPQ, wide band | 711 | 125 | 234 | 57 | 1,127 |
| CERTA --mapq pair | 687 | 152 | 232 | 62 | 1,133 |
| CERTA --mapq pair, wide band | 687 | 151 | 230 | 62 | 1,130 |
| CERTA --mapq pair --t2t | 645 | 138 | 225 | 60 | 1,068 |
| CERTA --mapq pair --t2t, wide band | 643 | 133 | 219 | 60 | 1,055 |
| minibwa | 687 | 113 | 223 | 55 | 1,078 |
| BWA-MEM2 | 676 | 129 | 232 | 59 | 1,096 |

cmrg, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 503 | 58 | 69 | 10 | 640 |
| CERTA, default MAPQ, wide band | 503 | 58 | 69 | 10 | 640 |
| CERTA --mapq pair | 464 | 79 | 70 | 14 | 627 |
| CERTA --mapq pair, wide band | 464 | 79 | 70 | 14 | 627 |
| CERTA --mapq pair --t2t | 454 | 74 | 65 | 14 | 607 |
| CERTA --mapq pair --t2t, wide band | 454 | 74 | 65 | 14 | 607 |
| minibwa | 483 | 60 | 67 | 12 | 622 |
| BWA-MEM2 | 478 | 67 | 68 | 12 | 625 |

cmrg, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 704 | 123 | 235 | 60 | 1,122 |
| CERTA, default MAPQ, wide band | 701 | 121 | 230 | 57 | 1,109 |
| CERTA --mapq pair | 675 | 148 | 230 | 62 | 1,115 |
| CERTA --mapq pair, wide band | 675 | 148 | 228 | 62 | 1,113 |
| CERTA --mapq pair --t2t | 638 | 137 | 222 | 60 | 1,057 |
| CERTA --mapq pair --t2t, wide band | 636 | 132 | 216 | 60 | 1,044 |
| minibwa | 683 | 106 | 221 | 55 | 1,065 |
| BWA-MEM2 | 673 | 125 | 230 | 58 | 1,086 |

cmrg, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 10 | 6 | 4 | 0 | 20 |
| CERTA, default MAPQ, wide band | 10 | 4 | 4 | 0 | 18 |
| CERTA --mapq pair | 12 | 4 | 2 | 0 | 18 |
| CERTA --mapq pair, wide band | 12 | 3 | 2 | 0 | 17 |
| CERTA --mapq pair --t2t | 7 | 1 | 3 | 0 | 11 |
| CERTA --mapq pair --t2t, wide band | 7 | 1 | 3 | 0 | 11 |
| minibwa | 4 | 7 | 2 | 0 | 13 |
| BWA-MEM2 | 3 | 4 | 2 | 1 | 10 |

v421, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 515 | 122 | 67 | 12 | 716 |
| CERTA, default MAPQ, wide band | 515 | 122 | 66 | 12 | 715 |
| CERTA --mapq pair | 484 | 113 | 66 | 9 | 672 |
| CERTA --mapq pair, wide band | 484 | 113 | 63 | 9 | 669 |
| CERTA --mapq pair --t2t | 386 | 148 | 61 | 13 | 608 |
| CERTA --mapq pair --t2t, wide band | 386 | 148 | 59 | 13 | 606 |
| minibwa | 472 | 71 | 55 | 8 | 606 |
| BWA-MEM2 | 457 | 96 | 64 | 9 | 626 |

v421, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 267 | 80 | 15 | 3 | 365 |
| CERTA, default MAPQ, wide band | 267 | 80 | 14 | 3 | 364 |
| CERTA --mapq pair | 243 | 57 | 15 | 1 | 316 |
| CERTA --mapq pair, wide band | 243 | 57 | 14 | 1 | 315 |
| CERTA --mapq pair --t2t | 208 | 101 | 15 | 5 | 329 |
| CERTA --mapq pair --t2t, wide band | 208 | 101 | 14 | 5 | 328 |
| minibwa | 243 | 48 | 15 | 2 | 308 |
| BWA-MEM2 | 243 | 69 | 14 | 0 | 326 |

v421, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 496 | 117 | 62 | 12 | 687 |
| CERTA, default MAPQ, wide band | 496 | 117 | 62 | 12 | 687 |
| CERTA --mapq pair | 469 | 108 | 62 | 9 | 648 |
| CERTA --mapq pair, wide band | 469 | 108 | 61 | 9 | 647 |
| CERTA --mapq pair --t2t | 371 | 144 | 58 | 13 | 586 |
| CERTA --mapq pair --t2t, wide band | 371 | 144 | 57 | 13 | 585 |
| minibwa | 460 | 68 | 53 | 8 | 589 |
| BWA-MEM2 | 443 | 92 | 62 | 9 | 606 |

v421, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 19 | 5 | 5 | 0 | 29 |
| CERTA, default MAPQ, wide band | 19 | 5 | 4 | 0 | 28 |
| CERTA --mapq pair | 15 | 5 | 4 | 0 | 24 |
| CERTA --mapq pair, wide band | 15 | 5 | 2 | 0 | 22 |
| CERTA --mapq pair --t2t | 15 | 4 | 3 | 0 | 22 |
| CERTA --mapq pair --t2t, wide band | 15 | 4 | 2 | 0 | 21 |
| minibwa | 12 | 3 | 2 | 0 | 17 |
| BWA-MEM2 | 14 | 4 | 2 | 0 | 20 |
