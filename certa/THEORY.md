# CERTA: the theory behind every bound

This document holds the written proofs for Aim 5 of
[EXACT_ALIGNMENT_PROPOSAL.md](../EXACT_ALIGNMENT_PROPOSAL.md). Every bound
the implementation uses is stated here and proven. Each statement also
names:
- the code that implements it;
- the oracle checks in `tests/test_certa.cpp` that test it against brute
  force;
- the mutants in `tests/mutants.tsv` that break it;
- for finite statements, the exhaustive check in `tests/verify_lemmas.cpp`.

**Review status.** Written, then cross-checked by:
- whole-reference oracles (4 seeds, ARM and x86);
- 111 mutants;
- exhaustive verification of the finite lemmas.

It has **not yet been reviewed by a second person**. §11 lists what a
reviewer should check first.

**Two bugs found while writing these proofs.** Both were in tier PR (§5) and
are now fixed and tested:
- S2 did not bound a second partner inside the winning anchor's window;
- the window search's band pad went negative for high thresholds.

## Contents

1. [Setting](#1-setting)
2. [The index](#2-the-index)
3. [End-to-end certificates (S0, S1, SR)](#3-end-to-end-certificates-s0-s1-sr)
4. [The affine local objective (Lemma L, SL, Lemma E)](#4-the-affine-local-objective)
5. [Pairs (Theorem P, S2)](#5-pairs-theorem-p-and-s2)
6. [Repeats (Theorems G and R, tier R)](#6-repeats)
7. [Interval output and XQ](#7-interval-output-and-xq)
8. [The checker (Theorem V)](#8-the-checker-theorem-v)
9. [Bounds proposed but not built (Lemma H, Proposition A, Theorem C)](#9-bounds-proposed-but-not-built)
   9a. [Mutants that survive by design](#9a-mutants-that-survive-by-design)
10. [Traceability](#10-traceability)
11. [For a reviewer](#11-for-a-reviewer)

---

## 1. Setting

**Reference.** G is the contigs concatenated, each surrounded by runs of
`N` (512 bases). `N` never matches anything, itself included.

**Read.** A read r has length L. Both orientations are aligned: strand 0 is
r, strand 1 is its reverse complement.

**Alignment.** An alignment A aligns the read span [x, y) to a reference
interval. Its transcript is a sequence of aligned bases (match or mismatch),
insertions (read bases with no reference base) and deletions (reference
bases with no read base). A starts and ends with an aligned base.

**Scores** (BWA-MEM's defaults):
- σ(A) = #matches − 4·#mismatches − Σ over gap runs of (6 + run length)
  − 5·[x > 0] − 5·[y < L].
- **Loss** λ(A) = (y − x) − σ(A) − 5·[x > 0] − 5·[y < L]. Event by event, λ
  is the sum of:
  - mismatches, 5 each;
  - insertion runs of k bases, 6 + 2k each;
  - deletion runs of k bases, 6 + k each.

  Matches cost nothing.
- **Clipping:** a clipped end of c bases costs c + 5 relative to an
  alignment of the whole read.

**Optimum and loci.**
- **OPT(r)** = max σ(A) over all alignments, both strands, all of G.
- **Same locus:** two alignments are at the same locus when they are on the
  same strand and their reference intervals overlap.
- **SUB** = the best score of an alignment not at the reported one's locus.

**Diagonals.** The *diagonal* of an aligned base is its reference position
minus its read position. A gap run of length k moves the diagonal by k.
An insertion's cells and a deletion's cells lie between the diagonals of
the aligned bases around them.

**Parts and occurrences.**
- **Parts (probes):** disjoint read intervals.
- **Exact occurrence of a part in A:** all its bases are matches, and no
  deletion lies between two of them. A then maps the part, unchanged, onto a
  reference substring at a single diagonal.
- **Occurrence of a part in G:** a position where the reference substring
  equals the part, with no `N`.

## 2. The index

**Lemma I (enumeration).** Let the index hold every N-free q-mer of G that
starts at a multiple of s. Then every occurrence of a part of length
≥ q + s − 1 is found, exactly once, by looking up the part's q-mer at
offset i, for i in [0, s).

*Proof.* An occurrence at position p covers p … p + ℓ − 1, with ℓ ≥ q + s − 1.
Exactly one i in [0, s) makes p + i a multiple of s. The q-mer at p + i lies
inside the occurrence (i + q ≤ ℓ), is N-free, and is indexed. Its hits h give
p = h − i. The other offsets i′ give positions p + i′ that are not sampled. ∎

*Code:* `count_parts` and `enumerate_parts` (core.h), `Checker::occurrences`
(check.cpp). *Mutants:* E3, C4.

**Proposition V0 (index validation).** Suppose:
- every index entry's key equals the reference q-mer at its position, and
  that position is a multiple of s;
- no position appears twice;
- the number of entries equals the number of N-free q-mers at multiples of s;
- within each directory bucket, keys are sorted.

Then the index is exactly the set of sampled N-free q-mers, and a range
lookup (bucket, then binary search) returns exactly the entries with the
queried key.

*Proof.* The entries are distinct sampled N-free positions, by the first
two conditions. They number as many as all such positions, so they are all
of them. Sorting within the bucket fixed by the key's high bits makes the
equal keys contiguous, and the binary search finds that run. ∎

*Code:* `validate_index` (check.cpp), run as `certa check-index`.
*Tests:* three valid indexes accepted, five corruptions rejected.
*Mutants:* C8, C9.

## 3. End-to-end certificates (S0, S1, SR)

Here d(r, ℓ) is the semi-global unit-cost edit distance of the whole read
at locus ℓ, with `N` costing 1. These proofs were first written in
FINDINGS §2; they are restated briefly.

**Theorem 1 (pigeonhole), Theorem 1′ (q-gram threshold t).**
- *Statement:* let S be parts (each ≥ q + s − 1) with all occurrences
  enumerated. Every locus with d ≤ |S| − t has ≥ t parts of S occurring
  exactly there, at diagonals within ±d of the locus's start diagonal.
- *Proof:* an edit touches at most one part, so d edits leave
  ≥ |S| − d ≥ t parts untouched, each occurring exactly. Its diagonal
  differs from the start diagonal by the indels before it, at most d.
  Lemma I enumerates it. ∎
- *Radius:* R = |S| − t.
- *Mutants:* E1, E3.

**Theorem 2 (band).** Banded DP over the cluster's diagonals ± 2R finds
d(r, ℓ) exactly for every such locus. *Proof:* its alignment stays within
±R of its start diagonal, which is within ±R of a member diagonal. ∎
*Mutant:* E4.

**Theorem 3 (best two, in any order).** Clusters are verified keeping the
two smallest distances b1 ≤ b2, each with limit b2 − 1.
- *Statement:* b1 is the minimum. When b1 < b2 ≤ R, b2 is the exact second
  minimum over other clusters.
- *Proof:* a pruned cluster has distance ≥ the b2 of that moment, and b2
  only decreases. ∎
- *Mutants:* E5, E6. *Test:* a metamorphic test reverses the order.

**Theorem 4 (SR).** Two distinct exact full-length matches prove d1 = 0 and
n_opt ≥ 2, under any scoring. *Mutant:* E2.

## 4. The affine local objective

**Lemma L.** Let the parts be disjoint, each of length ≥ 3, and let N(A) be
the parts inside A's span that do not occur exactly in A. Then
λ(A) ≥ 5·|N(A)|.

*Proof.* Charge each event's cost to the parts it touches, with every part
in N(A) receiving ≥ 5. A part touched by an event is one that contains one
of its read bases, or both read bases around it.
- **A mismatch** lies in at most one part: charge 5.
- **A deletion run** of k lies between two consecutive read bases, so it
  touches at most the one part containing both: charge 6 + k ≥ 7.
- **An insertion run** of k consecutive read bases touches parts
  p₁, …, pⱼ, in read order.
  - The parts strictly between p₁ and pⱼ lie wholly inside the run, so
    k ≥ 2 + (j − 2)·ℓmin for j ≥ 2.
  - Its cost 6 + 2k ≥ 10 + 2(j − 2)·ℓmin ≥ 5j holds whenever ℓmin ≥ 3
    (10 + 6(j − 2) = 6j − 2 ≥ 5j for j ≥ 2). For j = 1 the cost is ≥ 8.
  - Split the cost evenly: ≥ 5 per part.

A part in N(A) has a mismatch, an inserted base or a deletion between two
of its bases, so it receives ≥ 5. Each event's charges sum to at most its
cost. ∎

- *Sharpness:* parts of length 2 fail. A 4-base insertion over three parts
  costs 14 < 15; `verify_lemmas` finds it.
- *Verified* exhaustively for parts ≥ 3: all layouts, all transcripts of
  9-base reads, 52.8 M cases.
- *Mechanized* in Lean 4 for transcripts of any length
  ([`lean/Certa/LemmaL.lean`](lean/Certa/LemmaL.lean), `lemmaL`), with a
  simpler charging: every base is credited with its own loss, and an
  insertion run's opening cost of 6 is split as 3 to its first base and 3 to
  its last. The credits sum to at most the loss, and every inexact part of
  length ≥ 3 collects at least 5. The length-2 counterexample is
  `lemmaL_length_two`.
- *Corollary (U_t):* an alignment with fewer than t exactly occurring parts
  (of a set S of m parts) scores at most

  U_t = max over spans [x, y) of (y − x) − 5·max(0, #parts of S inside − (t − 1)) − 5[x > 0] − 5[y < L].

  *Code:* `local_bound` (core.h), `span_bound` (pair.h), `lemma_l`
  (check.cpp). *Mutants:* L1, L2, L4, L5, C1.

**Lemma S (span candidates).** In U_t it suffices to try x in
{0} ∪ {off + 1} and y in {L} ∪ {end − 1}, because the value changes only
where a part enters or leaves the span. For parts ≥ 5, the full span [0, L)
alone attains the maximum.
- *Proof of the second claim:* starting the span just inside part k
  excludes k + 1 parts and saves at most 5(k + 1). It loses
  off_k + 1 ≥ k·ℓmin + 1 bases plus the 5-point clip. With ℓmin ≥ 5 the loss
  is the larger. Symmetrically for the end. ∎
- *Verified:* 181,860 cases (layouts on L = 24, every counted set,
  extra ≤ 2). `span_bound` equals the brute-force maximum on all of them.
  Parts of length 4 give 3,885 counterexamples.
- *So:* mutant L6, which shifts the span start candidates, is
  **equivalent**. CERTA's parts are ≥ 17 bases.
- *Used:* `span_bound` (pair.h) returns the closed form L − 5 per counted
  part whenever the counted parts are disjoint, in order, in [0, L) and
  ≥ 5 long, and tries every span candidate otherwise (mutant P15, the
  closed form for parts of 4, is killed by `verify_lemmas`).
- *Mechanized* (the second claim) in Lean 4 for reads of any length
  ([`lean/Certa/LemmaS.lean`](lean/Certa/LemmaS.lean)): `lemmaS`, every
  span's bound at most the full span's for disjoint parts of length ≥ 5,
  any counted set and any extra; `lemmaS_length_four`, three parts of 4
  where a shorter span wins. The proof counts the parts dropped at each
  end: at most (x + 4) / 5 start before x.

**Lemma W (gaps).** An alignment with gap total G ≥ 1 has λ ≥ 6 + G. So an
alignment scoring above b has G ≤ pad(b) = max(0, L − b − 7).
- *Proof:* a gap run of k costs ≥ 6 + k. ∎
- *Note:* pad(b) = 0 when b ≥ L − 7: such an alignment has no gap at all.
  The window search in pair.h missed this clamp (§5).
- *Mechanized* in Lean 4 ([`lean/Certa/LemmaW.lean`](lean/Certa/LemmaW.lean),
  `lemmaW`): on Lemma E's alignments, a score above b gives gap total
  ≤ max(0, L − b − 7).

**Lemma C (chains).** Let all occurrences of the parts of S be enumerated,
sorted by diagonal, and split where consecutive diagonals differ by more
than pad. Let A have ≥ 1 exactly occurring part and gap total ≤ pad. Then:
- all of A's exact parts' occurrences lie in one chain C;
- every aligned base of A lies on a diagonal in [min C − pad, max C + pad];
- σ(A) ≤ ub(C), the Lemma L bound with C's parts free.

*Proof.* A's exact parts lie on diagonals within G ≤ pad of each other, so
no split falls between them. A's diagonals stay within G of each exact
part's. A's exactly occurring parts are among C's. ∎

*Mechanized* in Lean 4 ([`lean/Certa/LemmaC.lean`](lean/Certa/LemmaC.lean)):
its geometric core, `diag_drift` (two aligned bases of A lie on diagonals at
most G apart; a diagonal moves only at gaps), `diag_within_pad` (with Lemma
W, at most pad(b) apart when A scores above b) and `no_split_within` (no
split falls between two sorted diagonals at most pad apart); and its bound,
`lemmaC_bound` (if every exactly occurring part of A is among C's parts, A
loses at least 5 for each part of its span that is not C's: Lemma L with C's
parts free, so σ(A) ≤ ub(C)).

*Mutants:* L7, C3.

**Theorem L (tier SL, certified local optimum).**
- *Setup:*
  - floor = U_t;
  - pad = pad(floor);
  - chains built as in Lemma C;
  - every chain with ub > floor is evaluated by banded DP over its diagonals
    ± pad, keeping the best two scores s1 ≥ s2 (initially floor);
  - a chain with ub ≤ s2 is skipped.
- *Statement:* if s1 > floor, then OPT = s1.
- *Proof:* let A score above s1.
  - A has ≥ t exact parts (else σ ≤ floor) and gap total ≤ pad (Lemma W),
    so it lies in a chain C (Lemma C).
  - If C was evaluated, the DP saw A, so σ(A) ≤ s1.
  - If C was skipped, σ(A) ≤ ub(C) ≤ s2 ≤ s1. ∎
- *Mutants:* L1–L8.

**Corollary L-long (long reads, Aim 6, `tools/lr_ceiling --certify`).**
- *Setup:*
  - the parts are Theorem G's partition at τ, with t = 1. Every occurrence
    is enumerated (Lemma R), and stretches where no part starts are left
    uncovered;
  - lb is the score of a known alignment, and lb > U₁;
  - the chains are built with pad = L − lb − 7;
  - every chain whose Lemma L bound (its parts free) exceeds the best score
    so far is evaluated by banded DP over its diagonals ± pad.
- *Statement:* OPT = max(lb, the best DP score).
- *Proof:* Theorem L with floor = lb.
  - An alignment above lb has an exact part, because lb > U₁.
  - It has gap total ≤ pad (Lemma W), so it lies in a chain (Lemma C).
  - Either the DP saw it, or it scores at most the chain's bound, which is ≤
    the best found.
  - Lemma L needs only disjoint parts of length ≥ 3, not a tiling: uncovered
    stretches are simply not charged. ∎
- *Checks* (`--verify`):
  - every part's occurrence count is recounted from a second anchor, and U₁
    recomputed by brute force;
  - a full DP around the known locus never exceeds the proven OPT.

**The band DP row (`band_row`; both long-read certificates).** One row of the
banded Gotoh DP, W cells at a time (W = 8 with AVX2, else 4; NEON on
aarch64). With h0 = max(m, ins), the deletion state unrolls to
D[b] = max over b' < b of (h0[b'] + b') − b − 6, because a deletion never
follows its own opening better than it extends it.
- Taking b' ≤ b adds h0[b] − 6, which never beats h0[b]. So
  H[b] = max(h0[b], max over b' ≤ b of (h0[b'] + b') − b − 6): an inclusive
  prefix maximum, computed by log₂ W shifts within a vector and a carry
  across vectors.
- H, I and the row's best m equal the scalar recurrence's values exactly.
  Every h0 is a real score (a piece's start value is ≥ −5), so the two can
  differ only in unreachable sentinel values of D.
- *Tests:* exact equality with `band_row_scalar` on 20,000 random rows
  (sentinels, N bases and first rows included), on x86 with AVX2, x86 with
  SSE and aarch64. The certificates' outputs on 500 HiFi reads are
  identical with and without it.
- *Mutants:* LR28–LR31 (a missing shift step, a deletion 1 too cheap in the
  vector body and in the tail, the carry kept per lane).
- *Mechanized* in Lean 4 ([`lean/Certa/BandRow.lean`](lean/Certa/BandRow.lean),
  `vec_eq_scalar`): the scalar recurrence and the vector form (the inclusive
  prefix maximum with a carry) give the same row, for rows of any width,
  whenever every h0 is at least the sentinel.

**The traceback DP, vectorized (`affine_fill_vec`; tiers PR and SL, the
reported CIGAR).** `affine_align` fills a band of B diagonals with
traceback bits, then traces back. Its fill runs 4 cells at a time (32-bit
lanes) when the padded band lies inside the reference, and must give the
scalar fill's (`affine_fill_scalar`) values and traceback bytes in every
cell, sentinel cells included, so the traceback, and therefore every CIGAR,
is unchanged.
- *E and the match term* read only the previous row, so they are
  computed lane-wise as in the scalar loop.
- *F, the deletion state,* is F[b] = max(hc[b−1] − 7, F[b−1] − 1) with
  hc = max(h0, F) and F[0] = NEG. Opening from an F cell
  (F[b−1] − 7) is dominated by extending it (F[b−1] − 1), so
  F[b] = max(h0[b−1] − 7, F[b−1] − 1). Unrolled, this is exactly
  F[b] = max(NEG + 6, max over b' < b of (h0[b'] + b')) − b − 6 for every
  integer input, so unlike the band DP row above no "real score" premise is
  needed. It is an exclusive prefix maximum: two shifts within a vector and
  a carry that starts at NEG + 6.
- *The bits.* "F extends F" at b is F[b−1] − 1 > hc[b−1] − 7, read from the
  previous lane. The H source compares the same values in the same order as
  the scalar loop. The best end cell is taken lane by lane, in order, so ties
  go to the same, first cell.
- *Memory:* each vector stores its four traceback bytes whole. Past a row's
  end they land on the first cells of the next row, which is written later;
  after the last row they land past every cell the traceback reads, still
  inside `tb` ((LMAX + 1) · BMAX bytes; a band whose width is not a multiple
  of 4 is narrower than BMAX = 96).
- *Tests:* every traceback byte and the best cell equal the scalar fill on
  4,000 random bands per seed (reads with N, start ranges, bands up to
  BMAX). On 2 M HG002 pairs the SAM and fallback files are byte-identical.
- *Mutants:* E10–E13 (a deletion one point too cheap, the F-extension bit
  off by one point, ties to the last cell, the start range one cell short).
  E9, the same start-range fault in the scalar fill, is caught by the same
  test. An inclusive prefix (b' ≤ b) would change F but no cell value and no
  bit: it adds h0[b] − 6 to F[b], which never exceeds h0[b] (so neither hc[b]
  nor the H source changes), and when it raises F[b] the next lane's
  "F extends F" test, F[b] + 6 > hc[b] = h0[b], stays false. A first mutant
  of that form survived the tests for this reason and was replaced.
- *Mechanized* in Lean 4 ([`lean/Certa/AffineFill.lean`](lean/Certa/AffineFill.lean),
  `vec_eq_scalar`): within a row of any width, the scalar loop and the
  vector form give the same F, the same cell value and the same "F extends
  F" bit in every cell, for any integers.
- *16-bit lanes (`affine_fill_vec16`, 8 cells at a time; used whenever the
  band padded to 8 lies inside the reference).* The same operations with
  NEG = −16000 instead of −2²⁰. Every value the fill compares is either a
  real score, above −9 LMAX − 5 (a cell loses at most 9 points a row), or
  NEG plus an offset in (−9 LMAX − BMAX, LMAX); so no lane wraps, real scores
  stay above every NEG-derived value, and two NEG-derived values compare as
  their offsets do, whichever NEG they start from. Every comparison, bit and
  best cell is therefore the same; a best derived from NEG (no allowed
  start) is returned shifted back to the 32-bit NEG. *Tests:* the same 4,000
  bands per seed, every byte. *Mutants:* E22 (the prefix maximum without its
  4-lane step), E23 (NEG returned in 16-bit form).

**The edit-distance band (`band_distance_vec`; tiers S0/S1 and S2's
search).** `band_distance_scalar` is the unit-cost semi-global DP of a read
on a band of B ≤ 48 diagonals, stopping with limit + 1 once a whole row
exceeds the limit. Two vector forms give its result for every input.
- *16-bit lanes, 8 cells:* match and insertion read the previous row; the
  deletion step cur[b] = min(v[b], cur[b−1] + 1) unrolls to
  cur[b] = b + min over b' ≤ b of (v[b'] − b'), an inclusive prefix minimum.
  The row minimum is min v (every cell is at least it, and its cell keeps it).
- *8-bit lanes, 16 cells, capped at cap = limit + 1* (limit ≤ 60, the
  band's window inside the reference; every caller's limit is ≤ S2_MAX).
  Each recurrence term is non-decreasing in its argument and adds ≥ 0, so
  min(D, cap) satisfies the recurrence with every term capped (cells outside
  the band count as cap). Row minima never decrease (each cell is at least
  some cell of the row above), so the scalar result is
  min(final row minimum, limit + 1), and so is the capped one. It stops once
  every cell is at the cap.
- *`hamming_vec`:* 16 positions at a time, the count clamped to limit + 1
  as the scalar loop leaves it.
- *Tests:* 20,000 random bands and windows per seed against the scalar
  forms: reads with N, bands past either end of the reference, reads copying
  a window across a contig's N padding (N facing N costs 1), reads with a
  5–13-base deletion starting on diagonal 1 or 2 (the deletion chain crosses
  most lanes of one vector), limits up to 14.
- *Mutants:* E14 and E15 (a prefix-minimum step dropped), E16 (the cap one
  too high), E17 (lanes past the band not reset), E18 and E19 (N matching N).

**The score-only band (`band_score_vec`; pair windows, chains, tier R).**
Each row in one 8-lane pass: the same int16 values as `band_score_scalar`
(substitution, start range, insertion, band_del16's deletion maximum, the
row's best), lanes past the band masked out of the best. With a threshold
θ, only scores above θ must be exact: an alignment continuing from row i
gains at most a point per remaining base (≤ H + L − i), and one starting
later starts at −kClip (≤ L − i − kClip), so the DP stops once both bounds
are ≤ θ, returning the exact score when it is above θ and the best so far
(≤ θ) otherwise. `window_search` passes θ = max(best band so far, its
threshold t): a band at or below θ can change neither its result nor the
band it reports, which callers read only when the result is above t.
*Tests:* equality with the scalar form, and the threshold contract, on
4,000 random bands per seed. *Mutants:* E20 (a prefix-maximum step
dropped), E21 (stopping one point early).

**Other kernels, same results.** `repeat_detail::match_run` (tier R's
exact-copy checks) counts matching bases 8 at a time with a word test that
cannot carry between bytes (bases are 0–4); tested against the scalar count
on 20,000 windows, mutant E24. `lookup_bucket` scans a bucket of at most 64
entries linearly for the first key ≥ k and the first > k (the keys are sorted
within a bucket), mutant E25. Tier R's lookups go through `lookup_batch`
(the same ranges, prefetched), and `map_cpu` prefetches the next reads'
directory entries and buckets.

**Pruning the band (`band_best` with a floor; Corollary L-long).** The DP
needs only alignments above the best score so far, the floor. A cell (or
insertion state) of value v at read base f leads at most to
v + (L − 1 − f), since each later base adds at most 1; if that is ≤ the
floor, the cell is dead.
- Every cell on the best path to a live cell is live: its value is lower by
  at most the bases between them, and the bound grows by as many. So
  dropping dead cells leaves every live cell's value exact, and the end of
  any alignment above the floor is live. The DP's result is exact when it
  exceeds the floor, else at most the floor; max(floor, result) is
  unchanged.
- *The rows computed:* while a new piece can still lead above the floor
  (start + 1 + (L − 1 − f) > floor), the whole band. After that, a live
  cell's best path holds at least two aligned pairs (one pair and a gap is
  already below a new piece). Removing its last aligned pair moves the rest
  of the path up one row on the same diagonals and loses at most 1. So the
  cell above is live, and row f's live cells lie in the previous row's live
  range [c0, c1). Each row is computed on that range, and the dead cells at
  its ends are dropped.
- On HiFi, 35.5 % of the band's cells are live, and 35 % are computed.
- *Tests:* against the full DP on bands of random widths and offsets around
  each planted piece, with floors around the best (equal above the floor,
  at most the floor otherwise); every certificate on 242,534 HiFi reads is
  unchanged.
- *Mutants:* LR34–LR37 (the live test 1 too strict; the row range one cell
  short on either side; new pieces cut off 3 rows early).
- *Mechanized* in Lean 4 for alignment paths
  ([`lean/Certa/Prune.lean`](lean/Certa/Prune.lean)): `prefix_live` (every
  cell on a path to a live cell is live) and `above_live` (in a row where no
  new piece can lead above the floor, a live cell's path has two aligned
  pairs, and dropping one, `remove_pair`, gives a live path to the cell
  above). The scores are band_best's: start 0 or −5, pairs +1 / −4, a gap
  opening −7 and extending −1. That the DP's cell value is the best path
  ending there is Gotoh's recurrence, and is not mechanized.

**Jump (split) alignments (Aim 6 follow-up; `certa/longread.h`).** A long
read that crosses a structural variant has no good single local alignment.
- *The objective:* a jump alignment is a sequence of pieces A₁ … A_k, each a
  gapped local alignment of a read interval [a_j, b_j) on its own strand and
  locus, in read order (b_j ≤ a_{j+1}). It scores
  J = Σ σ(A_j) − P (k − 1) − 5 [a₁ > 0] − 5 [b_k < L], where σ has no clip
  penalties and read bases outside the pieces score 0.
  - k = 1 is the local objective.
  - P = 10 is how BWA-MEM scores a split read: two local alignments, each
    clipped at the junction.
  - The default is P = 30, so indels up to 24 bases stay gaps: a gap of k
    costs 6 + k, a jump P.
- *Lemma L-J:* if no piece holds an exact occurrence of a part of its strand,
  J ≤ U_J, the maximum over decompositions of
  Σ [(b_j − a_j) − 5 (parts of strand s_j inside [a_j, b_j))] − P (k − 1) −
  the end clips.
  - *Proof:* each piece is a local alignment of its read interval. Lemma L,
    applied to it without clip penalties, charges ≥ 5 for each inexact part
    inside it, disjointly. Summing over the pieces gives a decomposition's
    value, which is at most U_J. ∎
- *Computing U_J:* `lemma_l_jump` computes U_J exactly in one pass over the
  read. Its states are: inside a piece of strand s (with or without its
  current part covered from its first base), or between pieces.
  - *Verified:* against a brute force over all decompositions, on 210,000
    random two-strand part layouts and penalties
    (`tests/verify_lemmas.cpp`).
  - *How it compares with U₁:* U_J ≥ U₁. They differ only where the strands'
    partitions skip different stretches; on HiFi reads U_J − U₁ is 0.4 on
    average, up to 628.
- *A limit:* a novel insertion of length g, absent from the reference, can
  be a piece of its own that no part witnesses. Its bound is about
  g (1 − 5/ℓ) − 2P, where ℓ is the part length. So a long novel insertion
  cannot be certified exactly by exact parts; such a read keeps an
  interval.
- *Measured (bound side):* see proposal §6.6.

**Theorem L-J (the split optimum; `certify_jump`).**
- *Setup:*
  - lb ≤ OPT_J and lb > U_J;
  - *raising lb:* each round first combines the evaluated chains' pieces
    alone, which gives a real alignment. lb is raised to its score, so it
    stays ≤ OPT_J. Everything below uses the current lb, with
    pad = L − lb − 7;
  - chains of part occurrences are built per strand with the first pad. A
    later, smaller pad only coarsens that grouping, so Lemma C still holds
    for it. A *piece
    of chain C* holds an exact occurrence of a part of C, so its read
    interval contains that part. A piece holding no exact part of its strand
    is *unwitnessed*;
  - chains are *evaluated* in rounds. An evaluated chain C gets a read
    window W_C. Its pieces are then computed exactly: a banded DP over its
    diagonals ± pad, on the rows of W_C, in which a piece starts from the
    best prefix value minus P (or −5 / 0 for the first piece);
  - strand-1 chains run on the forward read against the reverse-complemented
    reference (diagonal N − L − d), so every chain advances with the read. A
    pass over the reversed read runs against the reversed reference and its
    complement, on the mirrored parts.
- *The values of a round:*
  - pre[a]: the best decomposition whose pieces all end at or before a. The
    evaluated chains' pieces are exact on their windows. The unevaluated
    chains' pieces are bounded by Lemma L with the union of their parts free
    (one type per strand). Unwitnessed pieces are bounded by Lemma L;
  - suf₂[b]: the same over pieces starting at or after b, from the same pass
    over the reversed read;
  - suf₁[b]: that pass without the unevaluated chains;
  - with st(a) = max(−5 [a > 0], pre[a] − P), en_i(b) = max(−5 [b < L],
    suf_i[b] − P) and V_C(a, b) = Lemma L for [a, b) with C's parts free.
- *Pruning:* pot(C) = max over [a, b) containing a part of C of
  st(a) + V_C(a, b) + en₁(b). Each round evaluates chains with pot(C) > lb:
  in the first round the one with the best single-piece bound, then the best
  by pot. The rounds end when no unevaluated chain has pot(C) > lb.
- *Windows:* when C is evaluated, W_C is the smallest interval holding every
  [a, b) that contains a part of C and has st(a) + V_C(a, b) + en₂(b) > lb.
- *The static path:* the first round's pre and suf₂ bound every alignment,
  because no chain is evaluated yet. So pot from them, with both-sided
  values, prunes in one step, and windows from them hold.
  - All chains with that pot above lb are then evaluated together, and the
    final pass follows.
  - It is taken when the rounds would cost more than half of it.
  - For reads where most chains need a DP, it costs the old single pass,
    instead of rounds that redo the evaluated chains each time.
- *The final pass:* the evaluated chains' exact pieces, with (upper) or
  without (lower) the unwitnessed pieces' Lemma L bounds.
- *Statement:* OPT_J ∈ [max(lb, lower), max(lb, upper)]; equal ends prove it.
  If a cell cap stops the rounds, OPT_J ∈ [lb, max(lb, U)], where U is the
  full value of the last pre pass.
- *Proof:* let A be a jump alignment above lb.
  - Each piece's gap total is ≤ pad: its loss is part of A's, which is
    < L − lb (Lemma W). A piece holding exact parts of its strand has them all
    in one chain C (Lemma C). So it lies on C's diagonals ± pad, and Lemma L
    with C's parts free bounds it.
  - *Windows, by induction over the rounds.* Claim: every piece of A of an
    evaluated chain lies in that chain's window. Let C be evaluated in round
    t, and p = [a, b) a piece of C in A.
    - The pieces of A before p are of three kinds: of chains evaluated
      earlier (in their windows by induction, so bounded by the exact DP); of
      chains unevaluated at t, C among them (Lemma L with the union free); or
      unwitnessed. So pre[a] bounds them.
    - Likewise suf₂[b] bounds the pieces after p.
    - So st(a) + V_C(a, b) + en₂(b) ≥ A > lb, and [a, b) ⊆ W_C.
  - *The pruned chains.* Suppose A uses a chain left unevaluated. Let p =
    [a, b) be A's last piece of such a chain C.
    - pre[a] bounds the pieces before p, as above.
    - The pieces after p are of evaluated chains (in their windows) or
      unwitnessed, so suf₁[b] bounds them.
    - So A ≤ pot(C) ≤ lb, a contradiction.
  - So every piece of A is either a piece of an evaluated chain in its window
    or unwitnessed, and A ≤ upper. The lower layer combines real pieces
    only, so it is attained.
  - Under the cap, the last pre pass bounds every piece of A in the same way,
    so A ≤ U.
  - *Raising lb* only shrinks the set of alignments to bound. A window set
    at an earlier, lower lb still holds every piece of an alignment above
    the current one.
  - *The static path:* every piece around a piece of C is bounded by the
    first round's values, which need no exact chain. So the window and the
    pruning arguments apply directly, with no induction. ∎
- *Why two suffix values:*
  - pruning bounds an alignment at its *last* piece of an unevaluated chain,
    so the unevaluated chains are needed on one side only. suf₁ ≤ suf₂, so it
    prunes more;
  - a window needs both sides, because an alignment may hold unevaluated
    pieces after a piece of C.
- *Verified* against brute force on 20,000 random layouts
  (`tests/verify_lemmas.cpp`): pre, suf₁ (through the reversed pass on the
  mirrored parts), pot (containment included) and the window each equal
  their definition.
- *Tests* (`tests/test_longread.cpp`):
  - a synthetic genome with planted deletions, small deletions, novel and
    diverged insertions, inversions, translocations, duplications, ragged
    ends, three loci (three chains in the optimum, the middle one
    reverse-complemented), and two loci joined by a diverged copy;
  - against an exact brute-force DP of the jump objective over the whole
    genome, both strands: every proven OPT_J equals it, and every interval
    contains it;
  - this holds from the known lb, from the weakest admissible lb (U_J + 1)
    and from lb = OPT_J − 1 (the narrowest band). Each is run three ways:
    with the default batch; with one chain per round, so that every pruning
    decision rests on the exact values of the chains evaluated before it;
    and on the static path;
  - under a small cell cap, the capped interval contains OPT_J;
  - every proven local optimum (Corollary L-long) equals the brute force
    with P = ∞.
  - *Mutants:* LR1–LR27. They cover pads; pruning; unwitnessed pieces; the
    strand-1 targets of both passes; the clip and jump penalties; Lemma L's
    charge and free parts; chain splitting; the prefix union type; the
    reversed pass's band, rows, parts and free flags; the suffix index; the
    windows (one-sided values, one base short at either end); the containment
    state; the capped interval; lb raised to a bound instead of a real
    alignment; and the static path's pot, set and windows.
- *Limits:*
  - an unwitnessed piece that pays off (a novel or diverged insertion) leaves
    an interval;
  - a read with a long lossy stretch, such as an insertion of repeat
    sequence:
    - while any chain covering the stretch is unevaluated, the union bound
      lets the prefix cross it almost free, so many chains look viable;
    - the unwitnessed bound (about 0.7 per base for parts of about 16 bases)
      leaves an interval about 0.7 × the stretch wide;
  - measured: proposal §6.6.

**SL's SUB bound (Aim 3).** SUB ≤ max(s2, the best alignment in the chosen
chain that starts before the reported one or at/after its end).
- *Proof:* an alignment at another locus outside the chosen chain is bounded
  as in Theorem L, by s2. One inside the chosen chain does not overlap the
  reported interval [p, p + span), so it starts < p or ≥ p + span. The start
  ranges cover both cases; the first is a superset. ∎
- *Mutants:* A4, A10, L1 (the planted floor reads put SUB exactly at the
  floor).

**Lemma E (edit distance to local score).** Let A be a local alignment whose
end-to-end completion has d ≥ 1 edits. Its completion aligns the clipped
bases diagonally, each counted as an edit in the worst case. Then
σ(A) ≤ edit_bound(L, d) = L − d − (d ≥ 2 ? 5 : 4).

*Proof.* L − σ(A) is the sum of event costs. Compare each event's cost to the
edits it adds:
- a mismatch: 5 for 1 edit;
- an insertion of k: 6 + 2k for k edits;
- a deletion of k: 6 + k for k edits;
- a clipped end of c bases: c + 5 for ≤ c edits.

So L − σ ≥ d + Σ(extras), with every extra ≥ 4. If d = 1 there is an event,
so L − σ ≥ d + 4. If d ≥ 2, either ≥ 2 events (extras ≥ 8), or one event
carrying ≥ 2 edits, which is a gap or a clip (extra ≥ 5). ∎

- *Tight:* one mismatch (d = 1), or a clipped end of d bases (d ≥ 2).
- *Verified* exhaustively: 8.1 M clipped transcripts of 8-base reads,
  attained for d = 1 … 4.
- *Mechanized* in Lean 4 for reads of any length
  ([`lean/Certa/LemmaE.lean`](lean/Certa/LemmaE.lean)): `lemmaE` (with the
  score defined as BWA-MEM's and shown equal to span − loss, `span_eq`),
  and `lemmaE_attained`, the bound attained for every d ≥ 1 by an alignment
  with an aligned base. The proof carries the "extra" through the
  transcript: ≥ 4 beyond the edits for one edit, ≥ 5 for two or more, and
  1 per further base of an insertion run already open.
- *Mutants:* A1, A5.

**Lemma O (overlap inside one band).** Let A and B be alignments whose
aligned bases all lie on diagonals of one band of w diagonals, and suppose
their reference intervals do not overlap. Then σ(B) ≤ L − σ(A) + w − 1.
- *Proof:* say B lies right of A. B's first aligned base x_B maps at or past
  A's last, so x_B + d_B ≥ y_A + d_A, with both diagonals in the band. Hence
  x_B ≥ y_A − (w − 1). B spans ≤ L − x_B ≤ L − y_A + w − 1 read bases, and
  y_A ≥ σ(A). Symmetrically on the left. ∎
- *In the code:* `process_read` uses the looser L − score + 2w + 6.
- *Why no real alignment reaches past Lemma E for L ≥ 98.* A certified S0/S1
  alignment with R ≤ 5 scores ≥ L − 8R ≥ 58 (an edit costs at most 8, a
  1-base insertion), and w ≤ 48. So the true in-band maximum
  L − score + w − 1 is below edit_bound(L, R + 1) ≥ L − R − 6.
- *So:* lowering the code's term (mutant A2) can tighten SUB_hi, but never
  below a real alignment. A2 is **equivalent** for these reads. It survives
  every test, as it should.

**Proposition E′ (S0/S1 interval and SUB).** Let a read be certified at
radius R, with best distance d1, tie count n_best and second distance d2.
Let d* = d1 if tied, else d2 if d2 ≤ R, else R + 1. Let w be the reported
cluster's band width. Then:
- OPT ≤ max(score, edit_bound(L, d*));
- SUB ≤ max(edit_bound(L, d*), min(score, L − score + 2w + 6)).

*Proof.*
- **Alignments outside the reported cluster's band.** Every locus outside
  the reported cluster has end-to-end distance ≥ d* (Theorems 1–3). Any
  alignment there scoring above edit_bound(L, R + 1) completes within R edits
  (Lemma E), so it lies in some other cluster, of distance ≥ d*.
- **Alignments inside the band** score ≤ score, the band optimum. A
  non-overlapping one there is bounded by Lemma O.
- **The rest:** Lemma E and monotonicity of edit_bound. ∎

*Mutants:* A1, A3, A5.

## 5. Pairs (Theorem P) and S2

**Objective.** σ(A₁) + σ(A₂) over proper pairs: opposite strands, with the
reverse mate's start 0 … D after the forward mate's.

**Lemma P0 (pair floor).** A proper pair with no exact enumerated part in
either mate scores ≤ U₁ + U₂ (Lemma L per mate). The pair floor is
max(U₁ + U₂, L₁ + L₂ − 7 − 47, L₁ + L₂ − max_loss). Pairs above it have
gap total ≤ pad = max(0, L₁ + L₂ − floor − 7) in each mate (Lemma W).

*Mutants:* P1 (family A plants a pair exactly at U₁ + U₂), P2 (family B
plants a 46-base deletion at pad 47).

**Lemma WS (window search).** Let y be a mate of length L with start range
[smin, smax]. Window parts: ~12-base pieces of y, scanned directly in the
reference. Let t = max(span_bound(all window parts), L − 7 − 47, need) and
pad_w = max(0, L − t − 7). Then window_search returns the best score of y
over alignments starting in [smin, smax], and every alignment it does not
return scores ≤ t.

*Proof.* Let A score above t.
- It has ≥ 1 exactly occurring window part (Lemma L, t = 1).
- Its gap total is ≤ pad_w (Lemma W).
- It ends before smax + L + pad_w, so the scan sees all its exact window
  parts.
- They lie in one chain (Lemma C), whose bound is ≥ σ(A), and the chain's
  bands cover A.
- Chains are processed best bound first and stop once a bound is ≤ the best
  found. So A is either found, or bounded by a found alignment. ∎

**Bug found (fixed).** The code computed pad_w = L − t − 7 without the clamp.
For t > L − 7 (a high `need`) the band was empty, so even an exact partner
was missed.
- *How it surfaced:* the new S2 oracle (seed 2, pair 110). The union window
  around a mate placement reported "nothing above 145" while a perfect mate
  (150) lay inside it.
- *Mutant:* P14 (the clamp removed). Inside Theorem P a partner's own
  anchor usually covers what the clamp-free search misses, so the pair
  oracle alone does not reliably catch it. The lemma is therefore also
  tested directly: windows around each threshold-family mate (centred on
  its best start, ending there, starting there), thresholds up to
  need = L − 2, against a brute-force oracle of the best alignment starting
  at each position.

**Theorem P (certified proper pair).** Statement and proof as in FINDINGS §2.
The proof uses only Lemmas P0, C and WS:
- every pair above the floor has an exact part in some mate x, so it passes
  through an anchor (a chain of x) and a partner in that anchor's union
  window;
- the anchors' bounds, and the `up` values of the evaluated ones, cover all
  pairs;
- certification requires every bound other than R1's own to be ≤ R1.

*Mutants:* P1–P3, P5–P7, P9, P10, P14. Families A, B, D and E put pairs at
the floor, the band pad, the window threshold and the window chain skip.
The direct Lemma WS test catches P14.

**Proposition P2 (S2).** S2 bounds every proper pair P′ that differs from R1's
pair P* = (A₁*, A₂*) in some mate, i.e. that mate does not overlap.

*Proof.* Let P′ score above the floor. It passes through some anchor C, a
chain of mate x.
- **C is not one of R1's anchors.** Then P′ ≤ the value C contributed to
  S2. That is its bound, if skipped; its dynamic bound; its `up`, if
  evaluated; or, if C held an earlier R1, the `upb` added when R1 moved on.
  The dynamic bound leaves out partners in finished chains, whose own
  evaluation bounded those pairs.
- **C is one of R1's anchors** (it realized P*). Then C's `up` ≥ P′, but
  `up` is about R1 itself and is not in S2. Split by mate x:
  - *x overlaps A_x\*:* then y does not overlap A_y\*. y starts in C's union
    window, before A_y\*'s start or at/after its end.
    P′ ≤ s_x + (the window search over those start ranges), where s_x is C's
    best for x.
  - *x does not overlap A_x\*:* then x starts in C before A_x\* or after it.
    P′ ≤ (the chain search over those start ranges) + max(fu, thr), which
    bounds any partner in the union window.
  - *Pairs at or below the anchor's `relevant`:* its bands were narrowed to
    that level, and S2 ≥ relevant covers them.
- **Finally,** S2 is capped at R1 = OPT. ∎

**Bug found (fixed).** Before Aim 5, S2 left out R1's anchors entirely. On
tandem-repeat pairs, a second partner in the winning window, often an exact
tie, went unbounded. For example, S2 was 246 while another pair also scored
275 = R1. The pair MAPQ was overstated.
- *How it surfaced:* the S2 oracle (best other proper pair, whole
  reference) failed 9 times on the old code.
- *Mutants:* P11, P12, P13.

## 6. Repeats

**Theorem G (greedy partition).** Let F_τ be the intervals of the read
occurring ≤ τ times. Scanning left to right and taking at each position s
the shortest interval of F_τ that starts there gives the maximum number of
disjoint intervals of F_τ.

*Proof.* F_τ is closed under extension: occurrences of a superstring are a
subset. So e(s) is non-decreasing. Interval scheduling by earliest end is
optimal, and the earliest-ending interval starting at or after s starts at
s. ∎

*Verified* against a DP optimum on 20,000 random instances.

*Mechanized* in Lean 4 ([`lean/Certa/TheoremG.lean`](lean/Certa/TheoremG.lean)):
- `greedy_max`: for any family of intervals of [0, L) closed under
  extension, the greedy picks are disjoint intervals of the family, and no
  such family is larger. The proof: the first interval [i, j) of any family
  starting at or after s extends to [s, j), so the greedy interval at s
  ends at or before j, and the rest follows by induction.
- `greedy_max_tau`: the same for F_τ itself, with occurrences counted in a
  reference text. A superstring occurs at most as often as the string
  (`occ_super`): extending to the right removes matches at each position,
  and one character to the left maps each match to the next position.

**Lemma R (anchor extension, tier R's implementation of Theorem G).** The
occurrences of [s, b) are exactly the hits h of an anchor q-mer at c, with
s ≤ c and c + q ≤ b, whose exact match extends left to s and right to b.
- So occ(b) for every b comes from one pass over the anchor's hits, each
  extended once.
- A part containing a count-0 q-mer (or an `N`) has no occurrence: it is
  *absent*.

*Proof.* An occurrence of [s, b) contains the anchor at offset c − s, so it
is a hit with both extensions. Conversely, such a hit is an occurrence. ∎

*Mutants:* R4, R5. R4 (the extension stops one base short of the read
end) is caught by planted reads (j):
- the last stretch of the read occurs > τ times, up to the read end, so the
  partition ends before it;
- one base short, that stretch looks like a rare part with no occurrence,
  and is charged as an absent part;
- U then falls to 130, below OPT = 135.

**Theorem R-int (tier R interval).** Let U be the Lemma L bound of tier R's
parts (t = 1), U_eff = max(U, L − 7 − 47), and pad = pad(U_eff). Chains are
evaluated best bound first until a bound is ≤ max(lo, U_eff). Then
OPT ∈ [lo, max(lo, U_eff)], where lo is the best score computed.

*Proof.* As Theorem L, with floor = U_eff, plus:
- the rare parts' occurrences are complete (Lemma R; a part whose
  occurrences do not fit ends the partition before it, so fewer parts, still
  sound);
- absent parts have no occurrence;
- lo is a real alignment's score. ∎

*Mutants:* R1, R2, R3.

**Theorem R (ties).** If the reference string under A* occurs at another
non-overlapping position (or on the other strand), the same transcript gives
an alignment of equal score at a different locus. *Mutant:* R6.

**Tier R's SUB interval (Aim 3).**
- **Lower end:** sub_lo is the score of an evaluated band lying wholly at
  another locus, or lo with a copy.
- **Upper end, for an exact read:** sub_hi bounds SUB as in Theorem R-int.
  Away chains are bounded by their best score. The own chain is split by
  start: before p, or at/after p + span.
- **Upper end, open intervals:** sub_hi = hi, because SUB ≤ OPT ≤ hi.

*Mutants:* A6–A9, A11, A12.

**Theorem R′ (window collapse).** Identical reference windows have identical
band maxima. This is not built: it was not needed (proposal §6.2).

## 7. Interval output and XQ

**Proposition Q (XQ).** BWA-MEM's single-end MAPQ formula
mapq_se(score, sub, L) is non-increasing in sub and non-decreasing in score.
- **XQ** = min(MAPQ, mapq_se(score, SUB_hi, L)) is therefore ≤ the formula
  at the true SUB, for every record whose SUB_hi is proven.
- *Verified* on 2.3 M (L, score, sub) triples.
- *The formula's calibration is not a theorem.* §6.3 of the proposal
  measures it.

**What is not a certificate.** The pair placement of interval mode (MR
records) chooses among proven alignments. Its MAPQ and placement are
heuristic, and its records say so (XT:Z:MR, XQ 0).

## 8. The checker (Theorem V)

**Theorem V (as built in `src/check.cpp`).** Assume an index validated by
Proposition V0. A record is accepted only if:
- AS rescores from the CIGAR within one contig;
- an alignment scoring at least XL is shown;
- for each bound b ∈ {XU, XH} below L, the witness (probes P, support t)
  satisfies:
  - (i) the probes are disjoint, in the read, and each ≥ q + s − 1 long;
  - (ii) U_t(P) ≤ b;
  - (iii) every chain of P's occurrences with ≥ t distinct probes and
    ub > b, built with pad(b) and searched by DP over its diagonals ± pad(b),
    holds no alignment above b. For XH, only alignments wholly left or right
    of the record's interval count;
- XQ ≤ MAPQ, and XQ ≤ mapq_se(AS, XH, L).

Every accepted claim is then true.

*Proof.* For XU:
- an alignment above b with fewer than t exact probes scores ≤ U_t ≤ b
  (Lemma L);
- otherwise it has gap total ≤ pad(b) (Lemma W) and lies in a chain with
  ≥ t probes, inside its diagonal band (Lemma C);
- the checker's own enumeration finds every occurrence (Lemma I, valid
  index), so that chain was built;
- if its bound exceeds b, the DP saw the alignment; otherwise σ ≤ ub ≤ b.

For XH, the same argument applies to alignments at other loci. The rest is
by recomputation. ∎

- *Witness sufficiency (why honest records pass):* each tier's bound was
  derived from the same probes by the same lemmas. Every U_t the checker
  computes is ≤ the mapper's bound (Lemma E ≥ Lemma L for S0/S1;
  sub_hi ≥ U_eff ≥ U for tier R; s2 ≥ floor for SL).
- *Tests and mutants:* §6.4 of the proposal, and C1–C9.

**Theorem V, PR mates (closes the Aim 4 gap).** A PR mate's XQ follows
BWA-MEM's pairing rule. It is accepted only if XQ ≤ MAPQ and
XQ ≤ max(q_se, min(q_pe, q_se + 40)), where:
- **q_se** = mapq_se(σ(A), XH, L). A is the mate's single-end alignment (XV if
  present, else the record's own), and XH is proven as above about loci not
  overlapping A. Without XH, q_se = 0.
- **q_pe** is the raw MAPQ of AS + AS(mate) − XS. It counts only if:
  - both records rescore and form a proper pair: one contig, opposite
    strands, the reverse mate's first aligned base 0 … D after the forward
    mate's, with D from the pair witness XY;
  - the pair witness (each mate's parts per strand, support t) proves
    S2 ≤ XS, where S2 is the best proper pair that differs from the shown
    one in some mate (Proposition P2's quantity). With
    pad = L₁ + L₂ − XS − 7, the witness must satisfy:
    - **(P1)** U_t(mate 1) + U_t(mate 2) ≤ XS;
    - **(P2–P3)** for each mate x, take each chain C of its parts'
      occurrences (split where diagonals are > pad apart) with ≥ t distinct
      parts. x's first aligned base then lies in [d_lo, d_hi + L_x − 1],
      where [d_lo, d_hi] are C's diagonals ± pad.
      - *Skip:* C may be skipped when ub(C) + B_y ≤ XS, or when
        x_C + B_y ≤ XS, with B_y a bound on y in the window: the larger of
        U_t(y) and the bounds of y's chains whose start range meets it, and
        x_C the best x found in C's segment.
      - *Otherwise:* the exact best proper pair through C must be ≤ XS. Pair
        x's best score per first aligned base f (an exhaustive DP over C's
        segment) with y's best per first aligned base r (one over the
        window), over r within D of f on the other strand (r ∈ [f, f + D]
        for a forward x, [f − D, f] for a reverse one).
      - *The exclusion:* when C lies on the shown x's strand, the value is
        max(x outside A_x paired with any y, any x paired with y outside A_y),
        "outside" meaning wholly left or right of the record's alignment.

Every accepted XQ is then at most the pairing rule at the true inputs: q_se
at (σ(A), SUB(A)), and q_pe at the shown pair's true margin over S2.

*Proof.* The rule max(a, min(b, a + 40)) is non-decreasing in a and b.
- q_se is at most its true value: SUB(A) ≤ XH (the XH case above), and
  mapq_se falls as SUB rises.
- q_pe is at most its true value once S2 ≤ XS. Suppose a proper pair P′ =
  (B_x, B_y) differs from the shown pair in some mate and σ(P′) > XS.
  - *Neither mate has t exactly occurring parts:* then
    σ(P′) ≤ U_t + U_t ≤ XS (Lemma L, P1).
  - *So some mate x has.* σ(B_y) ≤ L_y, so B_x loses less than
    L₁ + L₂ − XS, and its gap total is ≤ pad (Lemma W). The same holds for
    B_y.
  - *B_x is in a chain.* Its exact parts lie in one chain C (Lemma C), with
    ≥ t distinct parts, and σ(B_x) ≤ ub(C). Its first aligned base f lies in
    [d_lo, d_hi + L_x − 1], and its cells inside C's segment. So x's
    per-start best at f is ≥ σ(B_x).
  - *B_y is in the window.* B_y, on the other strand, starts at r within D
    of f (properness). Its cells, gap total ≤ pad, lie inside the window. So
    y's per-start best at r is ≥ σ(B_y).
  - *The skip is sound.* σ(B_y) ≤ U_t(y) if B_y has fewer than t exact
    parts. Otherwise B_y lies in a chain of y whose start range contains r.
    Either way σ(B_y) ≤ B_y's bound used in the skip.
  - *The exclusion.* If C is not on A_x's strand, B_x does not overlap A_x.
    Otherwise the shown pair is proper, so A_y lies on B_y's strand. Then
    "differs in some mate" means B_x lies wholly left or right of A_x, or
    B_y wholly left or right of A_y: one of the two terms counts it.

  In every case σ(P′) is at most a checked value, which is ≤ XS. ∎

- *Witness sufficiency (honest claims pass):*
  - The pair witness is exactly the parts Lemma P0's floor used
    (support 1), so (P1) holds: U₁ + U₂ ≤ floor ≤ S2.
  - Every pair that (P3) counts is a real proper pair differing from the
    shown one in some mate. The one exception is a pair across a contig
    boundary within D. So (P3)'s value is ≤ S2 ≤ XS.
  - So an honest XS is always proven, up to the work caps.
  - *An earlier version* added the best x in a chain to the best y anywhere
    in its window. It combined placements that are not a proper pair, and
    failed to prove 12 honest HG002 pairs; exact local oracles confirmed
    that their S2 held. Pairing per start fixed it.
- *Tests:* every certified pair of the PR oracle tests and of the threshold
  families goes through the check. Honest claims are accepted. False ones
  are rejected: XS one below the oracle's best other pair, the pair witness
  dropped, D below the pair's distance, XQ one above the rule, a mate moved.
- *Mutants:* C10–C19.
  - Families G–J make each of (P1), the pad and the two exclusion terms the
    only thing that refutes a false XS.
  - The skip bounds are tight only at ties. Ties have XQ 0 honestly, so the
    tests also claim XS = R1 − 1 there, with XQ raised to what that margin
    gives.

## 9. Bounds proposed but not built

**Lemma H (one-mismatch probes).** Let the parts have length ≥ 4. A part is
witnessed when an ungapped occurrence with ≤ 1 mismatch lies on A's path.
Then every unwitnessed part inside the span costs ≥ 7, disjointly.

*Proof.* Charge as in Lemma L, with these changes:
- **An insertion run of k** touching parts p₁ … pⱼ, holding k₁ … kⱼ of its
  bases:
  - charge 3 + 2kᵢ to each end part (i = 1, j);
  - charge 2kᵢ to each part it covers wholly;
  - the charges sum to 6 + 2k.
- **A wholly covered part** has kᵢ ≥ 4, so ≥ 8.
- **An unwitnessed part with no other event** is one of:
  - substitutions only: it needs ≥ 2 of them, ≥ 10;
  - a deletion inside it: ≥ 7, and a deletion touches no other part;
  - an insertion wholly inside it: ≥ 8;
  - an end share of a single base: its other bases lie on one diagonal, so
    the ungapped placement has ≤ 1 mismatch and the part *is* witnessed.

  So an unwitnessed part with an end share has kᵢ ≥ 2 (charge ≥ 7), or
  another event (≥ 5 more on top of ≥ 5). ∎
- *Verified* exhaustively: every layout at L = 10, and three parts of 4 at
  L = 12, 502 M cases.
- *Parts of length 3 fail:* a 7-base insertion over three parts costs
  20 < 21.
- *Mechanized* in Lean 4 for transcripts of any length
  ([`lean/Certa/LemmaH.lean`](lean/Certa/LemmaH.lean)): `lemmaH`, with
  `Witnessed` read off the transcript as `verify_lemmas`' `witnessed1`, and
  `lemmaH_length_three`, the counterexample above. Lemma L's credits
  suffice: every unwitnessed part of length ≥ 4 collects at least 7.
- *Witnessing:* the check reads it off the transcript. Coincidental extra
  witnesses only help.

**Proposition A (the ceiling for error-tolerant probes).** Proof in the
proposal §4.3. *Verified:* for e = 0 … 3 the constructed part costs e + 7,
and the nearest reference substring is e + 1 edits away. For e = 0, every
ungapped placement has ≥ 2 mismatches.

**Theorem C (completeness by refinement).** Proof in the proposal §4.6. Its
terminal step (BWT-SW on a suffix trie) is not built. What the implementation
guarantees is the interval [XL, XU] at whatever stage it stops.

## 9a. Mutants that survive by design

A mutant that changes no claim on any input cannot be killed by a test, so
it must be argued instead. Six of the 111 survive.
- **L6** (span start candidates shifted): Lemma S, §4.
- **A2** (Lemma O's term lowered): §4.
- **P2, P8, P13** (band pad − 2, narrowed band − 3, `relevant` term dropped)
  and **P4** (window scan without pad): Theorem P is redundant at those
  points, as follows.

**Proposition P3 (the band edge is covered twice).** Let P = (A_x, A_y)
score above the floor. Suppose A_x has gap total G_x > pad − 3, the last two
diagonals of mate x's band, and every mate has m ≥ 1 enumerated parts. Then:
- A_y is exact;
- the anchor of y's chain bounds P, whatever happened to x's chain.

*Proof.*
- **A_y has an exact part.** If it did not, σ_y ≤ U_y. With
  floor ≥ U_x + U_y and σ_x ≤ L_x − 6 − G_x, P > floor would need
  G_x < pad + 1 − (L_y − U_y) ≤ pad − 4 (L_y − U_y = 5m_y ≥ 5).
- **A_y is exact.** From P > floor and σ_x ≤ L_x − 6 − (pad − 2): σ_y ≥ L_y − 2,
  and any event costs ≥ 5.
- **So P passes through y's chain C_y,** with A_x in C_y's union window.
  - *If C_y is evaluated:* the window search for x (Lemma WS) uses its own
    pad_w, which the mutation does not touch, so it finds or bounds A_x.
  - *If C_y is skipped by its dynamic bound:* that bound replaces x-chains
    already done (possibly by the mutant, missing A_x) with U_x. So it is
    ≥ L_y + U_x. Meanwhile P ≤ L_y + L_x − 4 − pad ≤ L_y + U_x, because
    pad ≥ 5m_x − 4: with floor = U_x + U_y, pad = 5(m_x + m_y) − 7; with
    the band-term floor, pad = 47 and m_x ≤ 10. ∎

The same argument covers **P8**: there the edge is the narrowed band of an
anchor evaluated after R1 is known.

**P4 (scan without pad).** A union window's start range extends past the
anchor's diagonals by pad + L_x − 1. A pair above the floor has σ_x > L_x − 54,
so x is clipped by < 49 bases. Its partner therefore starts ≥ L_x − 50 > pad_w
before the end of the start range, and ends before smax + L_y, so the
union searches never need the pad.
- *The exact window* only yields realized pairs. If it misses a partner, the
  union bound makes upb > R1 and nothing is certified.
- *The left exclusion window of S2* only needs alignments that end before
  A_y\*'s start, which end before its scan end.

**P13 (`relevant` term).** S2's window term adds s_x + max(v, t) with
t ≥ need = relevant − s_x, so S2 ≥ relevant whenever one of the two start
ranges is non-empty. One always is: the union window spans ≥ D = 500 start
positions, more than a mate covers.

## 10. Traceability

| Statement | Code | Oracle tests | Mutants | Exhaustive |
|---|---|---|---|---|
| Lemma I, V0 | core.h `enumerate_parts`; check.cpp | all; index validation | E3, C4, C8, C9 | — |
| Thm 1, 1′, 2, 3, 4 | core.h `process_read`, `evaluate_clusters` | edit-distance oracle | E1–E7 | — |
| Lemma L, S, W, C | core.h `local_bound`; pair.h `span_bound`; check.cpp | local oracle | L1–L5, L7, C1, C3; L6 equivalent (Lemma S) | L, S |
| Thm L, SL SUB | core.h `certify_local` | local + SUB oracles, planted reads e–i | L1–L5, L7, L8, A4, A10 | — |
| Lemma E, O, Prop. E′ | core.h `edit_bound`, `process_read`; interval.h | [XL, XU] and SUB oracles, planted reads | A1, A2 (equivalent), A3, A5 | E |
| Thm P, Lemma WS, Prop. P3 | pair.h | pair oracle, threshold families A–F, Lemma WS against a start oracle | P1, P3, P5–P7, P9, P10, P14; P2, P4, P8 equivalent (§9a) | — |
| Prop. P2 (S2) | pair.h | other-pair oracle, family F | P11, P12; P13 equivalent (§9a) | — |
| Thm G, Lemma R, R-int, Thm R | repeat.h | local + SUB oracles, planted reads d and j | R1–R6, A6–A9, A11, A12 | G |
| Prop. Q (XQ) | interval.h `mapq_se` | — | C7 | XQ |
| Thm V (incl. PR mates) | check.cpp | falsified claims rejected (single records and pairs), real-data fault injection | C1–C19, mapper mutants K1–K6 | — |
| Lemma H, Prop. A | not built | — | — | H, A |

## 11. For a reviewer

Points worth checking first:
1. **Lemma L's insertion case** (ℓmin ≥ 3) and Lemma H's charging scheme.
   These are the two places where a counting slip would make a bound
   unsound. Both are exhaustively verified on small instances, but the
   general argument is the proof.
2. **Lemma C's diagonal argument with insertions.** The claim is that an
   alignment's cells stay within G of each exact part's diagonal.
3. **Proposition P2's case split.** It is new in Aim 5 and replaces a bound
   that was wrong.
4. **The equivalence arguments** for the mutants that survive by design
   (L6, A2), and their stated preconditions (parts ≥ 5; reads ≥ 98 bp).
5. **Mechanization.** Lemmas L, E, H, S (its second claim, which makes
   mutant L6 equivalent) and W, Lemma C's geometric core, Theorem G, the
   band pruning's two path facts and the vector DP row's identity are
   machine-checked in Lean 4 ([`lean/`](lean/README.md)). Theorem L-J's
   pruning (rounds, windows) is not.
