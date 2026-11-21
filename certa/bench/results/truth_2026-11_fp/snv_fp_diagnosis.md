# CERTA's extra false SNV calls without `--t2t` (2026-11-21)

DeepVariant 1.10.0 on the variant-calling subsets, hap.py against v5.0q (HG002) and v4.2.1 (HG005).
CERTA is the default build with the refit models and the wide band (`--intervals --iv-pair all --iv-rescue local
--mapq pair`, runs `S_wpair`). Scripts: [vc/snv_tiers.py](../../vc/snv_tiers.py),
[vc/snv_missing.py](../../vc/snv_missing.py), [vc/snv_ties.py](../../vc/snv_ties.py) and
[vc/lowend.sh](../../vc/lowend.sh) on a2, [low_levels.py](../../low_levels.py) on galaxy. The cost of the
mappings is in the experiment's table.

## The calls

| genome | CERTA false SNV calls | minibwa's | CERTA's only | wrong genotype at a true variant | no variant there | BWA-MEM2 wrong too |
|---|---|---|---|---|---|---|
| HG002 | 364 | 309 | 127 | 60 | 67 | 36 |
| HG005 | 113 | 71 | 73 | 33 | 40 | 10 |

Hap.py's false SNV calls, by the call's type as in its summary: on HG002, 8 of CERTA's and 11 of minibwa's indel
calls at true SNVs count as indel errors (minibwa makes CERTA's 8 too). A call within 10 bp of one of minibwa's
false calls counts as shared. Most of CERTA's extra false calls are not BWA-MEM2's: they do not come from the
alignment objective the two share.

## The reads behind the false alleles

The reads with the false allele (the called alt base; the ref base where the truth is homalt) and MAPQ >= 5 in
CERTA (DeepVariant's floor), by CERTA's tier and what minibwa does with the same read:

| genome | tier | reads | CERTA MAPQ 5-19, minibwa below 5 or not in the targets | minibwa has it at the site with MAPQ >= 5 | other |
|---|---|---|---|---|---|
| HG005 | PR | 320 | 140 | 170 | 10 |
| HG005 | IV | 95 | 38 | 6 | 51 |
| HG005 | MR | 97 | 34 | 35 | 28 |
| HG005 | S1 | 68 | 6 | 52 | 10 |
| HG005 | S0 | 53 | 7 | 33 | 13 |
| HG005 | R | 45 | 29 | 2 | 14 |
| HG005 | all | 678 | 254 | 298 | 126 |
| HG002 | PR | 438 | 241 | 168 | 29 |
| HG002 | IV | 307 | 14 | 255 | 38 |
| HG002 | S1 | 159 | 9 | 133 | 17 |
| HG002 | MR | 129 | 29 | 67 | 33 |
| HG002 | R | 50 | 0 | 28 | 22 |
| HG002 | S0 | 40 | 9 | 11 | 20 |
| HG002 | all | 1,123 | 302 | 662 | 159 |

"Other": CERTA MAPQ >= 20 where minibwa gives below 5, places the read elsewhere in the targets, or not in them.
Of the PR mates, 212 of 320 (HG005) and 314 of 438 (HG002) have two equally best single-end loci (XB 2): mates
in two-copy duplications that the pair placed. So do 29 of HG005's 45 tier R reads (HG002: 1 of 50). Of the
other tiers' reads, all but 3 (HG005) and 5 (HG002) have one best locus (XB 1).

At HG005's 73 sites CERTA has more reads with MAPQ >= 5 than minibwa (1,605 against 1,236). Of minibwa's, 174
are missing from CERTA's pileup: 133 of them CERTA places within 1 kb with MAPQ below 5 (MR 67, PR 51, IV 6, S1
6, R 3), 36 outside the targets or unmapped, 5 elsewhere with MAPQ >= 5.

## A lower floor for the least confident: not adopted

The model's Phred 3 and 4 are written as MAPQ 10 and 16 (`kBwaScale`), above DeepVariant's floor, while minibwa
gives MAPQ 0 to more of its least confident mates (test B: 0.081 % of mates, against CERTA's 0.056 % below
MAPQ 5). Calibration cannot place these levels: on real reads minibwa misplaces 19-62 % of the mates it gives
each MAPQ from 0 to 22 (test B, below). The test, by the MAPQ rewritten from the same records (`vc/bwscale.sh`), with a rule fixed
beforehand (adopt it if HG002's errors fall by >= 1 % with no rise on CMRG, then HG005 must not rise):

| run | HG002 SNV FN | HG002 SNV FP | HG002 indel FN | HG002 indel FP | HG005 SNV FN | HG005 SNV FP | HG005 indel FN | HG005 indel FP | CPU (HG002 subset) | memory |
|---|---|---|---|---|---|---|---|---|---|---|
| CERTA `--mapq pair` | **1,458** | 364 | **1,321** | 497 | **484** | 113 | **63** | 9 | 3,311 s | 21.6 GiB |
| Phred 3 below MAPQ 5 | 1,483 | 369 | 1,322 | **488** | 502 | **110** | 64 | **8** | 3,311 s | 21.6 GiB |
| Phred 3 and 4 below MAPQ 5 | 1,527 | **358** | 1,334 | 494 | 553 | 114 | 69 | 9 | 3,311 s | 21.6 GiB |
| *minibwa* | *1,436* | *309* | *1,281* | *494* | *472* | *71* | *55* | *8* | *75,512 s (30×)* | *9.5 GiB* |

HG002 CMRG, SNV and indel errors: 838 and 292, 848 and 289, 861 and 288 for the three CERTA rows. The false calls
barely move (364 -> 369 and 358 on HG002, 113 -> 110 and 114 on HG005), while the missed ones rise by 18-69: the
reads at Phred 3-4 carry more true alleles than false ones. The rule fails on HG002.

## Misplaced share by MAPQ (galaxy `low_levels.py`)

The 150-bp locus rule; per MAPQ value (minibwa) or per XP, the model's Phred (CERTA `--mapq pair`):

```
B (real)  minibwa  0:53%/25992 1:62%/2502 2:46%/277 3:42%/361 4:61%/242 5:60%/280 6:34%/850 7:38%/392 8:41%/249
                   9:43%/344 10:48%/518 12:19%/561 15:37%/484 18:41%/840 20:41%/3032 22:36%/2880 24:15%/649 27:18%/418
          CERTA XP 0:82%/2424 1:65%/3263 2:49%/12288 3:46%/3328 4:45%/15941 5:32%/2511 6:33%/5128 7:19%/3965
                   8:14%/3926 9:14%/3521 10:10%/3334 11:7%/4630 12:5%/6333
H (HG005) minibwa  0:27%/49341 1:61%/2492 2:51%/275 3:39%/371 4:54%/201 5:55%/265 6:33%/890 10:47%/511 20:41%/2855
          CERTA XP 0:78%/1453 1:55%/4754 2:40%/15870 3:34%/5051 4:22%/33837 5:30%/4277 6:29%/5271 7:14%/4954
```

## What it says

- CERTA's extra false SNV calls are not its least confident reads: excluding them loses true variants and keeps
  the false calls.
- On HG002, most of the false-allele reads at those sites are reads minibwa also counts (662 of 1,123; on HG005,
  298 of 678). The calls differ in the pileups around paralogous sites: CERTA keeps more reads at MAPQ >= 5 there
  and drops others that minibwa keeps. That is a property of the pair model's ranking in duplications, where the
  truth sets' labels (certified unique on one haplotype) are thin.
- With `--t2t` the T2T placement separates such reads, and CERTA's errors equal or beat minibwa's on both genomes.
