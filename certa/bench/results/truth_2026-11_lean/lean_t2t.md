# A lean `--t2t`: CHM13 lookups only where they can matter (2026-11-21)

## Why

`--t2t` maps every read pair a second time, on T2T-CHM13, and the pair model reads where CHM13 puts each mate.
It is CERTA's most accurate mode, at about 2× minibwa's CPU and 4.6× its memory (HG005: 7,544 s and 46.1 GiB,
against 3,981 s and 24.8 GiB without it, and minibwa's 3,741 s and 10.1 GiB).

The first idea for a cheaper mode was CHM13 sink contigs in one index: the CHM13 sequence GRCh38 lacks, added
to the GRCh38 index to attract reads that belong there. A check on the existing `--t2t --feature-tags` runs of
the variant-calling subsets ([vc/snv_t2t.py](../../vc/snv_t2t.py), 2026-11-07 build) ruled that out. At the
false SNV calls only CERTA makes ([snv_fp_diagnosis.md](../truth_2026-11_fp/snv_fp_diagnosis.md)), CHM13 puts
the reads with the false allele (MAPQ >= 5) as follows:

| genome | reads | T2T here | elsewhere in GRCh38 | CHM13 only (a sink's) | of those, T2T MAPQ >= 20 | other |
|---|---|---|---|---|---|---|
| HG002 | 1,123 | 634 | 325 | 143 (12.7 %) | 57 (5.1 %) | 21 |
| HG005 | 678 | 249 | 345 | 80 (11.8 %) | 25 (3.7 %) | 4 |

A sink can only take the "CHM13 only" reads. The larger group, reads that CHM13 places at another GRCh38 locus,
needs the CHM13 placement itself.

## Where T2T disagrees

Across all records of the subsets ([vc/t2t_tiers.py](../../vc/t2t_tiers.py); disagreement: class 2 or 3 with T2T
MAPQ >= 20), the tiers S0 and S1 outside segdups and low-mappability regions hold 77 % of the records and
17-18 % of the confident disagreements. On the truth sets ([t2t_skip.py](../../t2t_skip.py), the refit's
features of 2026-11-20), the gap CERTA records for S0/S1 mates (XD, else XO, minus XE) is absent for 96 % of
them: no second locus within the radius. Skipping the CHM13 lookup for a pair whose two mates are both S0/S1
with gap >= G (the mapped mates in skipped pairs; the misplaced mates in them, of all misplaced mates; and those
of them T2T flags, of all misplaced mates T2T flags):

| G | mates skipped, A / B / H | misplaced mates skipped, A / B / H | flagged misplaced mates lost, A / B / H |
|---|---|---|---|
| 0 | 85.1 / 85.9 / 83.6 % | 5.3 / 15.8 / 15.0 % | 9.1 / 10.9 / 11.8 % |
| 3 | 83.6 / 84.3 / 82.0 % | 1.1 / 5.2 / 4.0 % | 2.2 / 3.5 / 3.2 % |
| 6 | 82.3 / 82.9 / 80.7 % | 0.5 / 2.8 / 1.6 % | 0.9 / 1.8 / 1.4 % |
| 11 | 81.0 / 81.4 / 79.3 % | 0.3 / 2.1 / 1.1 % | 0.6 / 1.3 / 1.0 % |

On real reads T2T flags most misplaced mates: 20,981 of test B's 27,828 and 20,708 of HG005's 28,206 (on
simulated test A, 51,672 of 405,221). At G = 11 a lean run would keep 99 % of those flags while mapping a fifth
of the mates on CHM13.

## The rule, fixed before the fits

1. Cross-validation on HG002 (fit on the even pairs of tests A and B, report the odd pairs) of the T2T model
   (depth 6, 400 rounds) with the skip emulated for G = 3, 6 and 11 (`truth_mapq_pair.py fit --t2t-skip G`:
   skipped mates get class 4, no T2T placement). A G passes if, on test B's odd pairs, its misplaced mates per
   million at matched N (99.75 / 99.44 / 95.58 %) are at most 5 % above the full T2T model's (`cv_t2t.h` of
   2026-11-20) at every N. The smallest passing G (the most pairs skipped) is chosen.
2. The skip for that G is built into CERTA: a pair whose mates' results on GRCh38 (before the pair placement)
   are both certified S0/S1 with gap >= G is not mapped on CHM13, and its mates get class 4. The lean model is
   fitted on the lean build's own features of A and B (as the refit was), and evaluated once on HG005. It is
   adopted unless it is more than 5 % worse than the refit T2T model at two or more of HG005's four N.
3. The lean build's CPU and memory are measured on the HG005 run and on test B back to back, and the
   variant-calling subsets are called with DeepVariant.

(Step 2 was amended before any result: the first draft evaluated the emulated model on HG005, which would have
spent the held-out genome on the emulation rather than on the build.)

## Cross-validation on HG002 (step 1)

Test B, the odd pairs: misplaced mates per million among each model's N most confident. The cost is test B
mapped on galaxy's 64 cores with the T2T lookups each model reads (`--t2t`, and `--t2t-skip G`), once each
([time_B_tags.log](time_B_tags.log); `--feature-tags` without `--mapq pair`, so MAPQ limit 8, about 4 % less
CPU than with it; records to /dev/null; build ~/certa-lim/b/lean1; the fits emulate the skip):

| model | 99.75 % | 99.44 % | 95.58 % | MAPQ >= 40 | CPU (B) | memory (B) |
|---|---|---|---|---|---|---|
| full T2T (`cv_t2t.h` of 2026-11-20) | **60.8** | **31.3** | **8.7** | 96.17 % | 6,099 s | 42.9 GiB |
| G = 3 | 77.6 (+28 %) | 52.0 (+66 %) | 17.2 (+98 %) | 94.50 % | 5,604 s (-8.1 %) | 42.9 GiB |
| G = 6 | 67.5 (+11 %) | 40.8 (+30 %) | 17.3 (+99 %) | 94.50 % | 5,644 s (-7.5 %) | 42.9 GiB |
| G = 11 | 64.6 (+6 %) | 37.9 (+21 %) | 17.3 (+99 %) | 94.11 % | 5,666 s (-7.1 %) | 42.9 GiB |

Test A (simulated), at 94.13 % and 83.10 %: full 110.2 and 23.4; G = 3 196.4 and 66.4; G = 6 161.5 and 66.6;
G = 11 146.1 and 66.2. No G passes: each is more than 5 % worse at every N of test B. The rule stops here: the
skip was not built into CERTA, and HG005 was not used.

## Why, and what the skip would have saved

- **The flags it loses are the confident errors.** At 95.58 % the full model leaves about 133 misplaced mates
  among the odd pairs' top N, every lean model about 264. The difference, about 130, is about the number of
  flagged misplaced mates in skipped pairs among the odd pairs (264 on all of test B at G = 11, half of them
  odd). These mates look unique on
  GRCh38: certified S0 or S1, with no second locus within the radius. Only their CHM13 placement shows they
  are misplaced; without it they rank with the confident mates, and at loose N they are most of the errors.
- **It would have saved little CPU.** On test B the full `--t2t` takes 6,099 s, against 3,453 s and 21.6 GiB
  without it at the same MAPQ limit ([wide_cost.md](../truth_2026-11_refit/wide_cost.md)). Skipping 81-84 % of
  the mates saves 433-495 s of the 2,646 s it adds (16-19 %), and no memory. With `--mapq pair`, back to back
  twice (the lean model a stand-in: CPU only): full 6,349 s and 6,348 s, `--t2t-skip 11` 5,934 s and 5,934 s
  (-6.5 %), all four at 42.8-42.9 GiB ([time_B.log](time_B.log)). Pairs that are easy on GRCh38 are cheap on
  CHM13 too: the cost of `--t2t` is in the hard fifth of the mates, where most of its signal is.

## What it says

- Neither sinks nor skips make `--t2t` cheap. Sinks can take the reads from copies GRCh38 lacks, a small part
  of what the CHM13 placement shows. Skipping the pairs that are easy on GRCh38 loses the reference-caused errors
  that look unique there, and saves little.
- The cost of `--t2t` is the CHM13 mapping of the reads that are hard on GRCh38 (its CPU), and a second
  full index (its memory: 42.9 against 21.6 GiB on test B). Those are what a cheaper mode would have to cut.
- The emulation is `truth_mapq_pair.py fit --t2t-skip G` (and `eval`, which reads G from the header); the tallies
  are [t2t_skip.py](../../t2t_skip.py). The C++ of the skip (an option `--t2t-skip`, class 4 for skipped pairs)
  was written and tested before the rule decided, and then removed.
