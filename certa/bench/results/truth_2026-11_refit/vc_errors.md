# DeepVariant errors of the refit, by region

DeepVariant 1.10.0 on the variant-calling subsets (HG002 against v5.0q and CMRG, HG005 against v4.2.1), hap.py
PASS errors per run from `vc/errors_table.py`. CERTA's runs: the default build's placements, which the refit build
reproduces, so only MAPQ differs between them. The cost of each run is in VALIDATION_RESULTS.md, "The pair model
refit on the default build's placements".


v5q, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,522 | 348 | 1,377 | 517 | 3,764 |
| CERTA --mapq pair, 2026-11-05 trees | 1,474 | 367 | 1,332 | 494 | 3,667 |
| CERTA --mapq pair, refit | 1,465 | 365 | 1,333 | 498 | 3,661 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 1,339 | 357 | 1,307 | 500 | 3,503 |
| CERTA --mapq pair --t2t, refit | 1,309 | 362 | 1,304 | 502 | 3,477 |
| minibwa | 1,436 | 309 | 1,281 | 494 | 3,520 |
| BWA-MEM2 | 1,419 | 342 | 1,329 | 519 | 3,609 |

v5q, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 584 | 83 | 139 | 32 | 838 |
| CERTA --mapq pair, 2026-11-05 trees | 543 | 95 | 137 | 37 | 812 |
| CERTA --mapq pair, refit | 537 | 88 | 135 | 36 | 796 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 516 | 80 | 127 | 30 | 753 |
| CERTA --mapq pair --t2t, refit | 515 | 84 | 124 | 30 | 753 |
| minibwa | 556 | 71 | 127 | 30 | 784 |
| BWA-MEM2 | 548 | 82 | 130 | 29 | 789 |

v5q, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 1,448 | 331 | 1,363 | 515 | 3,657 |
| CERTA --mapq pair, 2026-11-05 trees | 1,386 | 355 | 1,316 | 492 | 3,549 |
| CERTA --mapq pair, refit | 1,377 | 352 | 1,318 | 496 | 3,543 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 1,264 | 346 | 1,293 | 498 | 3,401 |
| CERTA --mapq pair --t2t, refit | 1,237 | 350 | 1,290 | 500 | 3,377 |
| minibwa | 1,388 | 290 | 1,271 | 490 | 3,439 |
| BWA-MEM2 | 1,367 | 322 | 1,318 | 514 | 3,521 |

v5q, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 74 | 17 | 17 | 2 | 110 |
| CERTA --mapq pair, 2026-11-05 trees | 88 | 12 | 19 | 2 | 121 |
| CERTA --mapq pair, refit | 88 | 13 | 17 | 2 | 120 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 75 | 11 | 17 | 2 | 105 |
| CERTA --mapq pair --t2t, refit | 72 | 12 | 17 | 2 | 103 |
| minibwa | 48 | 19 | 13 | 3 | 83 |
| BWA-MEM2 | 52 | 20 | 13 | 4 | 89 |

cmrg, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 714 | 129 | 239 | 60 | 1,142 |
| CERTA --mapq pair, 2026-11-05 trees | 695 | 156 | 234 | 61 | 1,146 |
| CERTA --mapq pair, refit | 687 | 152 | 232 | 62 | 1,133 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 649 | 140 | 222 | 58 | 1,069 |
| CERTA --mapq pair --t2t, refit | 645 | 138 | 225 | 60 | 1,068 |
| minibwa | 687 | 113 | 223 | 55 | 1,078 |
| BWA-MEM2 | 676 | 129 | 232 | 59 | 1,096 |

cmrg, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 503 | 58 | 69 | 10 | 640 |
| CERTA --mapq pair, 2026-11-05 trees | 470 | 84 | 73 | 16 | 643 |
| CERTA --mapq pair, refit | 464 | 79 | 70 | 14 | 627 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 453 | 73 | 66 | 12 | 604 |
| CERTA --mapq pair --t2t, refit | 454 | 74 | 65 | 14 | 607 |
| minibwa | 483 | 60 | 67 | 12 | 622 |
| BWA-MEM2 | 478 | 67 | 68 | 12 | 625 |

cmrg, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 704 | 123 | 235 | 60 | 1,122 |
| CERTA --mapq pair, 2026-11-05 trees | 683 | 153 | 232 | 61 | 1,129 |
| CERTA --mapq pair, refit | 675 | 148 | 230 | 62 | 1,115 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 643 | 139 | 220 | 58 | 1,060 |
| CERTA --mapq pair --t2t, refit | 638 | 137 | 222 | 60 | 1,057 |
| minibwa | 683 | 106 | 221 | 55 | 1,065 |
| BWA-MEM2 | 673 | 125 | 230 | 58 | 1,086 |

cmrg, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 10 | 6 | 4 | 0 | 20 |
| CERTA --mapq pair, 2026-11-05 trees | 12 | 3 | 2 | 0 | 17 |
| CERTA --mapq pair, refit | 12 | 4 | 2 | 0 | 18 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 6 | 1 | 2 | 0 | 9 |
| CERTA --mapq pair --t2t, refit | 7 | 1 | 3 | 0 | 11 |
| minibwa | 4 | 7 | 2 | 0 | 13 |
| BWA-MEM2 | 3 | 4 | 2 | 1 | 10 |

v421, all:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 515 | 122 | 67 | 12 | 716 |
| CERTA --mapq pair, 2026-11-05 trees | 481 | 111 | 67 | 10 | 669 |
| CERTA --mapq pair, refit | 484 | 113 | 66 | 9 | 672 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 404 | 151 | 63 | 14 | 632 |
| CERTA --mapq pair --t2t, refit | 386 | 148 | 61 | 13 | 608 |
| minibwa | 472 | 71 | 55 | 8 | 606 |
| BWA-MEM2 | 457 | 96 | 64 | 9 | 626 |

v421, segdups:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 267 | 80 | 15 | 3 | 365 |
| CERTA --mapq pair, 2026-11-05 trees | 236 | 60 | 15 | 2 | 313 |
| CERTA --mapq pair, refit | 243 | 57 | 15 | 1 | 316 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 211 | 98 | 15 | 5 | 329 |
| CERTA --mapq pair --t2t, refit | 208 | 101 | 15 | 5 | 329 |
| minibwa | 243 | 48 | 15 | 2 | 308 |
| BWA-MEM2 | 243 | 69 | 14 | 0 | 326 |

v421, all difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 496 | 117 | 62 | 12 | 687 |
| CERTA --mapq pair, 2026-11-05 trees | 465 | 106 | 63 | 10 | 644 |
| CERTA --mapq pair, refit | 469 | 108 | 62 | 9 | 648 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 389 | 147 | 60 | 14 | 610 |
| CERTA --mapq pair --t2t, refit | 371 | 144 | 58 | 13 | 586 |
| minibwa | 460 | 68 | 53 | 8 | 589 |
| BWA-MEM2 | 443 | 92 | 62 | 9 | 606 |

v421, not difficult:

| run | SNV FN | SNV FP | indel FN | indel FP | errors |
|---|---|---|---|---|---|
| CERTA, default MAPQ | 19 | 5 | 5 | 0 | 29 |
| CERTA --mapq pair, 2026-11-05 trees | 16 | 5 | 4 | 0 | 25 |
| CERTA --mapq pair, refit | 15 | 5 | 4 | 0 | 24 |
| CERTA --mapq pair --t2t, 2026-11-05 trees | 15 | 4 | 3 | 0 | 22 |
| CERTA --mapq pair --t2t, refit | 15 | 4 | 3 | 0 | 22 |
| minibwa | 12 | 3 | 2 | 0 | 17 |
| BWA-MEM2 | 14 | 4 | 2 | 0 | 20 |
