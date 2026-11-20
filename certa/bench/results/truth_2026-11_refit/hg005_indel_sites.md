# HG005: indel errors of CERTA's refit `--mapq pair` run that minibwa does not make

DeepVariant 1.10.0 on the variant-calling subsets, hap.py against v4.2.1 (error sites within 10 bp of each other
count as one). Per site and run, the reads that cover it (`samtools view -F 0x904`): all, with MAPQ >= 20, with
an indel in their CIGAR within 5 bp of the site, those with MAPQ >= 20, and those soft-clipped within 10 bp of it.

| site | error | ref | alt | BWA-MEM2 wrong too | CERTA | minibwa | BWA-MEM2 |
|---|---|---|---|---|---|---|---|
| 16:71132958 | FP | TGAG | T | no | 36 / 20 / 7 / 3 / 0 | 44 / 20 / 6 / 3 / 0 | 45 / 19 / 6 / 2 / 0 |
| 18:663360 | FN | GTTGT | G | yes | 11 / 1 / 6 / 1 / 0 | 7 / 1 / 5 / 1 / 0 | 7 / 6 / 4 / 3 / 0 |
| 20:5471754 | FN | CT | C | no | 31 / 31 / 3 / 3 / 0 | 39 / 31 / 11 / 3 / 0 | 32 / 31 / 4 / 3 / 0 |
| 20:12382164 | FN | A | AAAAAAC | yes | 24 / 24 / 15 / 15 / 0 | 24 / 24 / 15 / 15 / 0 | 24 / 24 / 15 / 15 / 0 |
| 20:16663209 | FP | A | T,ATTTTTT | no | 19 / 19 / 11 / 11 / 0 | 19 / 19 / 11 / 11 / 0 | 19 / 19 / 11 / 11 / 0 |
| 20:18898847 | FN | A | ATGTATATATAT | yes | 21 / 21 / 2 / 2 / 0 | 21 / 21 / 3 / 3 / 0 | 21 / 21 / 2 / 2 / 0 |
| 20:19850990 | FN | AAAGGAAGGA | AAAGGAAGG,A | yes | 20 / 20 / 3 / 3 / 0 | 20 / 20 / 4 / 4 / 0 | 21 / 21 / 4 / 4 / 0 |
| 20:38081181 | FN | A | ATTTTTCCTATC | no | 28 / 28 / 0 / 0 / 0 | 30 / 30 / 4 / 4 / 0 | 30 / 30 / 3 / 3 / 0 |
| 20:38166305 | FN | A | ATCCCTCCCTCG | yes | 44 / 44 / 0 / 0 / 2 | 52 / 51 / 5 / 4 / 2 | 48 / 48 / 1 / 1 / 2 |
| 20:38790115 | FN | G | GCCATGCTGTGA | yes | 34 / 33 / 0 / 0 / 0 | 38 / 37 / 3 / 2 / 0 | 35 / 35 / 0 / 0 / 0 |
| 20:42618200 | FN | CA | C | no | 28 / 25 / 6 / 3 / 0 | 27 / 25 / 5 / 4 / 0 | 26 / 26 / 4 / 4 / 0 |
| 20:50366199 | FN | G | GCTCACTGCAGC | no | 28 / 28 / 1 / 1 / 8 | 35 / 35 / 10 / 10 / 6 | 30 / 30 / 4 / 4 / 7 |
| 20:51753454 | FN | A | AAAC | yes | 32 / 32 / 3 / 3 / 0 | 33 / 33 / 3 / 3 / 0 | 34 / 33 / 3 / 3 / 1 |
| 20:52354715 | FN | G | GTATGTTATGTT | yes | 38 / 38 / 0 / 0 / 0 | 43 / 41 / 3 / 3 / 1 | 40 / 39 / 1 / 1 / 1 |
| 20:55245013 | FN | CA | C | yes | 32 / 31 / 5 / 5 / 2 | 33 / 33 / 6 / 6 / 2 | 33 / 33 / 6 / 6 / 2 |
| 20:63732368 | FN | G | GTGTGTGTGCGC | yes | 53 / 53 / 0 / 0 / 2 | 57 / 56 / 8 / 7 / 2 | 54 / 54 / 0 / 0 / 2 |
| 21:45989998 | FN | AGGGGTCCCC | A | no | 29 / 29 / 1 / 1 / 0 | 30 / 30 / 3 / 3 / 0 | 29 / 29 / 2 / 2 / 0 |
| 5:1430797 | FN | C | CTCACCTCCCGC | no | 30 / 30 / 0 / 0 / 28 | 39 / 39 / 24 / 24 / 13 | 38 / 38 / 21 / 21 / 15 |

Cells: reads / MAPQ >= 20 / with the indel / with it and MAPQ >= 20 / clipped near it.
CERTA's indel error sites 70, minibwa's 59, BWA-MEM2's 66;
CERTA only (not within 10 bp of a minibwa error): 18.
