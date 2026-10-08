# Exact, Provably Optimal Read Alignment: a Proposal (CERTA-X)

*Prepared October 2026. This continues [RESEARCH_ROADMAP.md](RESEARCH_ROADMAP.md) Part 3. It builds on the
measurements and theorems in [certa/FINDINGS.md](certa/FINDINGS.md), and on the literature in
[SWOT_ANALYSIS.md](SWOT_ANALYSIS.md). CERTA-X is a working name.*

> **Status of the claims.** Each result is labelled:
> - **[P+T]** proved, implemented in `certa/` and checked against a brute-force oracle;
> - **[P]** proved in this document;
> - **[S]** proof sketch, to be written out formally;
> - **[H]** hypothesis or projection, to be measured (§7).
>
> Numbers marked *measured* come from the HG002, HG001 and HG005 runs in `certa/bench/results/`.

## Contents

1. [Summary](#1-summary)
2. [Why exactness, and where the current prototype stops](#2-why-exactness-and-where-the-current-prototype-stops)
3. [Problem statement and output contract](#3-problem-statement-and-output-contract)
4. [Theory](#4-theory)
5. [The algorithm](#5-the-algorithm)
6. [Cost model and projections from our measurements](#6-cost-model-and-projections-from-our-measurements)
7. [Evaluation plan and aims](#7-evaluation-plan-and-aims)
8. [Long reads](#8-long-reads)
9. [Relation to prior work](#9-relation-to-prior-work)
10. [Risks](#10-risks)
11. [References](#11-references)

---

## 1. Summary

**Goal.** A read aligner in which every output is a *mathematical claim with a proof*. For each read
it returns an alignment and an interval [ℓ, u] that provably contains the optimal score OPT over the
entire reference. When ℓ = u, the alignment is the proven global optimum. When u < T, the read is
proven unalignable. No output is ever a heuristic guess presented as an answer.

**What is already established** (CERTA prototype, `certa/`):
- exact-part certificates for end-to-end, local and paired alignment (Theorems 1–4, L and P);
- an oracle-tested implementation with byte-identical CPU, ARM, x86 and GPU output;
- **95.4 % of HG002 reads** with a proven optimum;
- **variant-calling F1 above BWA-MEM2** in 11 of 12 comparisons across HG001, HG002 and HG005;
- **4.0×** BWA-MEM2's speed in paired mode.

**What this proposal adds:**
1. **A completeness theorem.** Refinement drives the interval to a single point for every read
   (Theorem C). It is phrased as an interval certificate, so the per-read work can be capped without
   ever overclaiming.
2. **A sharp design result for affine scoring (Proposition A).** Error-tolerant seeds raise the
   per-part bound only from 5 to at most e + 7, because one long deletion destroys a part cheaply.
   *More, shorter exact parts* is the effective lever: +5 points per part at near-linear cost.
3. **An optimal, linear-time part choice** on an FM-index or a second hash index (Theorem G). It
   moves the certifiable ceiling past the q = 22 barrier measured in FINDINGS §4.
4. **Count-based proofs for repeats.** Tie proofs need no enumeration (Theorem R), and verification
   is collapsed across identical copies.
5. **Verifiable output.** Each certificate is a witness that a small, independent checker can
   validate (Theorem V), so a laboratory can audit a run without trusting the mapper.

**The honest boundary.** Conditional lower bounds rule out a uniformly cheap exact aligner (§4.1). The
cost of exactness must grow with how far a read is from the reference. CERTA-X is designed so that
this cost is paid only by the few reads that need it, and by those reads in proportion to their
distance.

## 2. Why exactness, and where the current prototype stops

### 2.1 Measured weaknesses of heuristic aligners

On 1.89 M certified HG002 reads, every tested aligner places some reads at a locus that is
provably at least one edit worse than the optimum. These are *measured* counts
(`certa/FINDINGS.md` §1):

| Tool | Provably misplaced | …with MAPQ ≥ 20 | Exact multi-copy reads given MAPQ ≥ 20 |
|---|---|---|---|
| BWA-MEM2 | 177 | 13 | 0 |
| minibwa | 420 | 38 | 0 |
| strobealign | 4,984 | 800 | 704 |
| minimap2 `-x sr` | 2,237 | 54 | 0 |
| bowtie2 | 921 | 61 | 0 |

The faster the tool, the more it misses. That is the speed/rigour trade-off an exact method
removes for the reads it certifies.

### 2.2 How far the certificate reaches today (measured)

| Setting | Reads with a proven optimum | Remaining |
|---|---|---|
| Single-end, HG002 (NovaSeq X) | 94.4 % (95.65 % with pass 3 and SL) | heuristic fallback |
| Paired-end, HG002 | 95.4 % | heuristic fallback |
| Paired-end, HG001 / HG005 (NovaSeq 6000, 151 bp) | 92.9 % / 92.9 % | heuristic fallback |

**Why the rest fail.** After the local-alignment tier (single-end HG002, 2 M reads), 5.3 % of
reads stay uncertified:
- **3.6 % are repeat-limited:** fewer than 6 parts fit the hit budget.
- **1.8 % are high-loss:** all 6 parts were searched, but the best alignment loses ≥ 25 points.

BWA-MEM2's own alignments of those uncertified reads lose this much (cumulative share of the
106,868 reads):

| loss < | 10 | 25 | 30 | 40 | 45 | 50 | 60 | 100 | ≥ 100 | unmapped |
|---|---|---|---|---|---|---|---|---|---|---|
| share | 17.2 % | 34.0 % | 39.6 % | 52.9 % | 57.4 % | 61.4 % | 67.8 % | 88.1 % | 7.6 % | 4.3 % |

**The q = 22 barrier.** A 150 bp read has at most 6 disjoint 22-mers. With 5 points per destroyed
part, no exact-part certificate can prove an alignment that loses ≥ 30 points (≥ 25 with t = 2).
Shorter parts need a second index.

**Paired-end.** Of the pairs sent to the fallback:
- 30 % are improper or have an unmapped mate;
- among the proper ones with small loss, 83 % have a mate with MAPQ < 20, so they sit in repeats.

**Where the time goes (measured, paired-end HG002).**
- **Fallback cost per read:** the reads sent to the fallback are also the expensive ones. They are
  12.9 % of reads (both mates of every pair with an uncertified mate), yet cost 46 % of minibwa's
  total CPU time: about 0.32 ms per read, against minibwa's 0.088 ms average.
- **CPU per read today:** minibwa 88 µs; BWA-MEM2 423 µs; CERTA's own tiers (certificate and pair
  pass, without the fallback) 30 µs.
- **Break-even budget:** an exact replacement handles only the 4.6 % uncertified reads. CERTA's
  total CPU stays at or below minibwa's while 30 µs + 0.046·x ≤ 88 µs, that is x ≲ **1.2 ms per
  residual read**. With x ≲ 2.4 ms it stays 3× below BWA-MEM2. (These compare CPU time; CERTA also
  uses the GPU.)

## 3. Problem statement and output contract

**Notation.**
- **Reference G:** the contigs concatenated, separated by runs of `N`. `N` never matches.
- **Read r** of length L, aligned on either strand.
- **Alignment A:** maps the read span [x, y) to a reference interval with an edit transcript.

**Score** (BWA-MEM defaults, clipping included; this is the objective CERTA already certifies,
Theorem L):

σ(A) = (#matches) − 4·(#mismatches) − Σ_gaps (6 + length) − 5·[x > 0] − 5·[y < L]

**Loss:** λ(A) = (y − x) − σ(A) − 5·([x > 0] + [y < L]) ≥ 0. A perfect, unclipped alignment has
λ = 0.

**Optimum:** OPT(r) = max σ(A) over all alignments, both strands, the whole reference.

**Locus:** two alignments are at the same locus if they lie on the same strand and their reference
intervals overlap.
- **SUB** is the best score at any locus that does not overlap the reported one.
- **n_opt** is the number of pairwise non-overlapping loci that attain OPT.

**Threshold T:** the minimum reportable score (BWA-MEM's default `-T 30`).

**Output contract.** For every read, CERTA-X returns:
1. an alignment A* together with an **interval [ℓ, u] ∋ OPT**, where ℓ = σ(A*);
2. a **tie statement**: n_opt = 1 proven, n_opt ≥ 2 proven, or unknown;
3. an **interval for SUB**, from which MAPQ gets a certified *lower* bound (MAPQ falls as SUB
   rises);
4. a **witness ω** that an independent checker can verify (Theorem V).

**Exact** means ℓ = u. **Proven unalignable** means u < T. **Pairs** use the pair objective of
Theorem P: σ(A₁) + σ(A₂) over proper pairs. BWA-MEM's pairing rule compares the best proper pair
with the sum of the mates' individual optima minus an unpaired penalty, and becomes exactly
decidable once both quantities are certified.

## 4. Theory

### 4.1 What cannot be done, and the design principle it implies

- **No subquadratic edit distance.** Edit distance cannot be computed in strongly subquadratic time
  unless SETH fails (Backurs & Indyk 2015).
- **No cheap indexing with errors.** For text *indexing* with mismatches or edits, conditional lower
  bounds exclude indexes that answer queries with k errors in time and space polynomial in m and k
  for k = Θ(log n) (Cohen-Addad, Feuilloley & Starikovskaya 2019).

So **an exact aligner cannot be uniformly fast.** Its cost must grow with the number of differences
between the read and the reference.

**Design principle.** Make the cost *distance-sensitive and adaptive*:
- pay almost nothing for the ≥ 93 % of reads within a few edits;
- for the rest, pay in proportion to the loss λ that the read actually has;
- learn λ while searching, by narrowing [ℓ, u].

### 4.2 The framework: branch and bound over witness classes **[P]**

**Probes and witnesses.** A *probe* is a read interval p with a matching rule, such as exact
occurrence or ≤ 1 mismatch. Its occurrence set Occ(p) is the set of reference positions where it
matches. An alignment A is *witnessed* by p if p lies inside A's span and occurs, under its rule, at
the reference position that A assigns to it.

**Theorem C1 (interval certificate).** Let Π be a set of probes. Suppose:
- (i) every occurrence of every probe in Π is enumerated;
- (ii) every *unwitnessed* alignment satisfies σ(A) ≤ U(Π);
- (iii) every witnessed alignment with σ(A) > U(Π) lies inside the verification band of its
  witness's occurrence.

Let W be the best score found by exact band DP over the witnessed candidates. Then:
- if W > U(Π), OPT = W exactly;
- otherwise OPT ∈ [ℓ, U(Π)] for any realized alignment score ℓ.

*Proof.* Every alignment is witnessed or not. A witnessed alignment scoring above U lies in a
verified band, so it scores ≤ W. An unwitnessed alignment scores ≤ U by (ii). ∎

Theorems 1, L and P of the prototype are instances of C1. The rest of this section supplies better
bounds U, cheaper witnesses, and a refinement that always terminates.

### 4.3 Bounds for unwitnessed alignments

**Lemma L, general parts [P+T].** Let the parts be disjoint read intervals, each of length ≥ 3.
For any alignment A, every part lying inside A's span that does not occur exactly in A adds ≥ 5 to
λ(A), disjointly. Hence:

U_t(Π) = max over spans [x, y) of (y − x) − 5·max(0, #{parts inside [x, y)} − t + 1) − 5·([x > 0] + [y < L])

bounds every alignment with fewer than t exactly matching parts. The prototype proves this for
equal-length 22-base parts (`include/certa/core.h`, `local_bound`); the argument only uses part
length ≥ 3, so it covers parts of any length.

**Lemma H, one-mismatch probes [S].** Let the parts have length ≥ 4, and let a part be witnessed
when an *ungapped* occurrence with ≤ 1 mismatch lies on A's path. Then every unwitnessed part inside
the span adds ≥ 7 to λ(A), disjointly.

*Sketch.* An insertion of k bases touching j parts costs 6 + 2k (its bases also lose their +1).
Charge it as 3 + 2kᵢ to each of the two end parts, which hold kᵢ of its bases, and 2kᵢ to any part
it covers entirely (≥ 8 for parts of length ≥ 4). The charges add up to its cost. Now take an
unwitnessed part with no other event:
- *substitutions only:* it needs ≥ 2 of them (≥ 10);
- *a deletion strictly inside it:* ≥ 7, and a deletion touches no other part;
- *an insertion wholly inside it:* ≥ 8;
- *an insertion whose share is a single base at the part's end:* the other bases lie on one
  diagonal, so that ungapped placement has ≤ 1 mismatch and the part is witnessed. Unwitnessed
  therefore means kᵢ ≥ 2 (charge ≥ 7) or another event in the part (≥ 10).

FINDINGS §4 states this bound; the full case analysis still needs to be written out.

**Proposition A, the affine ceiling for error-tolerant probes [P].** For probes that tolerate up to
e edits, *no* valid per-part bound exceeds e + 7. For substitution-only (Hamming) probes of any
tolerance, none exceeds 7.

*Proof (by construction).* A bound is valid only if *every* unwitnessed part in *every* alignment
costs at least that much, so one counterexample settles it.
- **The construction:** take a part p = (AC)^m with m ≥ 2(e + 1). Let the reference hold its left
  half, then G^{e+1}, then its right half, with T-runs on both flanks. Let A align p across that
  region with one deletion of the e + 1 G's in the middle and no other events. The part costs
  6 + (e + 1) = e + 7.
- **Why the part is unwitnessed:** any reference substring overlapping the region is at edit
  distance > e from p. If it contains the whole G-run, each G costs an edit, because p has no G. If
  it does not, it misses half of p's AC-pattern, or swaps it for T's.
- **Hamming probes:** take e = 0, a single 1-base deletion, cost 7. Every ungapped placement shifts
  one half of the alternating pattern by one, giving about m mismatches. ∎

| Probe type | Per-part bound | Enumeration cost per part |
|---|---|---|
| exact | 5 (tight) | one lookup |
| ≤ 1 mismatch (Hamming) | 7 (Lemma H; tight by Proposition A) | about 3ℓ lookups, or 2 search-scheme passes |
| ≤ e edits | at most e + 7 | grows exponentially in e (search schemes) |

**Consequence.** Under BWA-MEM's affine scoring, tolerating more errors buys at most one point per
extra tolerated error, at exponential cost. Adding parts buys five points per part at near-linear
cost. **The effective lever is more, shorter, exact parts.** This is the main design result of the
proposal, and it explains why the 1-mismatch idea in FINDINGS §4 helps only a little.

**Absent-part bound [P].** A part with Occ(p) = ∅ cannot occur exactly anywhere. It costs nothing
to enumerate and still counts in U.
- A sequencing error inside a unique region usually makes the part containing it absent.
- Adapter or contaminant reads have almost all parts absent, so U collapses and OPT < T is proven
  at once.

**Pair bound (Theorem P) [P+T].** A proper pair with no exact enumerated part in either mate scores
at most U₁ + U₂. The partner's window search uses ~12-base parts scanned directly in the reference.
Together these double the effective number of parts.

### 4.4 Choosing parts optimally **[P]**

Let occ(s) be the number of occurrences of string s in G (both strands). Fix a per-part occurrence
limit τ. The family F_τ = {[i, j) : occ(r[i..j)) ≤ τ} is closed under extension on either side,
because occurrences of a superstring are a subset.

**Theorem G (optimal greedy partition).**
- **The algorithm:** scan left to right; at position s, take the shortest prefix r[s..e(s)) with
  occ ≤ τ, then continue at e(s).
- **Optimality:** this yields the maximum number of disjoint parts from F_τ.
- **Cost:** O(L) backward or forward extension steps on a (bidirectional) FM-index.

*Proof.* e(s) is non-decreasing in s, because extending an interval to the left keeps it in F_τ.
Picking the interval that ends earliest is optimal for interval scheduling, and here that is always
the interval starting at the current position. Each extension step is charged to the read position
it consumes. ∎

Absent parts are picked automatically, since occ = 0 ≤ τ. τ is the refinement knob: a larger τ
gives shorter parts, more parts and a lower U, at the price of more occurrences to enumerate.

**Proposition B, the ceiling-cost trade-off [P under a random-sequence model; H on real data].**
With equal parts of length ℓ on a 150 bp read, the certifiable loss ceiling is
5·(⌊150/ℓ⌋ − t + 1). The expected occurrences are ⌊150/ℓ⌋ · 2|G| / 4^ℓ, with |G| = 2.95 × 10⁹:

| ℓ | parts | ceiling, t = 2 | ceiling, t = 1 | expected occurrences per read |
|---|---|---|---|---|
| 22 (today) | 6 | 25 | 30 | 0.002 |
| 18 | 8 | 35 | 40 | 0.7 |
| 16 | 9 | 40 | 45 | 12 |
| **15** | **10** | **45** | **50** | **55** |
| 14 | 10 | 45 | 50 | 220 |
| 13 | 11 | 50 | 55 | 970 |
| 12 | 12 | 55 | 60 | 4,200 |

**The sweet spot is ℓ ≈ 15:** the ceiling nearly doubles, from 25 to 45, at about 55 occurrences
per read.
- With t = 2, chains hit by a single part need no DP at all (Theorem 1′). This is what keeps the
  random hits of short parts cheap.
- On the real genome, the k-mer spectrum is skewed by repeats. Theorem G then lengthens parts there
  automatically, and the measured occurrence counts are Aim 1's first output.

### 4.5 Repeats: count, do not enumerate

**Theorem R, tie proof from a count [P].** Let A* be a reported alignment and ref(A*) the reference
string it covers. If ref(A*) occurs at two or more non-overlapping positions (counting either
strand), then n_opt ≥ 2 whenever A* is optimal. This holds for every scoring scheme.

*Proof.* Each occurrence supports the same transcript, so it gives an alignment with the same
score at a different locus. ∎

- **Cost:** one FM-index count of ref(A*), plus locating two occurrences.
- **Relation to SR:** this generalizes tier SR, which required an exact full-read match, to reads
  with edits.

**Theorem R′, verification collapse [P].** Two candidate windows with identical reference sequence
have identical band maxima. So the DP runs once per *distinct* window string; hashing a window costs
about 200 byte reads, against about 15,000 DP cells for verifying it. This lets the per-read hit
budget grow by up to two orders of magnitude in near-identical repeats (segmental duplications,
recent L1s) at the same DP cost [H]. Diverged copies, such as old Alus, do not collapse. But they
produce chains with few parts, which the bound U_t skips.

### 4.6 Completeness and monotone refinement **[P]**

**Refinement steps.** Each step either certifies the read or strictly enlarges the probe set:
- **R1** raise τ (shorter parts, Theorem G);
- **R2** add the mate (Theorem P);
- **R3** add one-mismatch probes (Lemma H);
- **R4** raise the hit budget with window collapse (Theorem R′);
- **R5 (terminal)** exact local alignment against the whole reference by branch and bound on the
  suffix trie (BWT-SW: Lam et al. 2008), seeded with the current best ℓ.

**Theorem C (completeness).**
- **Correctness at every stage:** the interval [ℓ, u] is valid after every step.
- **Monotonicity:** ℓ never decreases, u never increases, and a certified read stays certified.
- **Completeness:** R5 terminates with ℓ = u. Every alignment is a path in the suffix trie
  starting at the root, and a branch is pruned only when an admissible bound shows it cannot beat ℓ.
  An admissible bound is: the current DP column's best value plus one point per remaining read base,
  lowered further by absent parts.
- **Worst case:** O(L·|G|), as SETH requires.
- **Practical cost:** grows with ℓ's distance from L. For reads with large loss, R5 behaves like
  approximate search with about λ/5 errors, which is exponential, consistent with §4.1.

**Policy consequence.** Run refinements in order of cost. Stop when the interval is a point, when
u < T, or when a per-read budget is reached. In the last case, report the interval [ℓ, u]. The
reported alignment is the best found, and the output says *exactly* how much better an unseen
alignment could be. This replaces today's unbounded heuristic fallback with a bounded statement.

### 4.7 Verifiable output **[P]**

CERTA-X is designed as a *certifying algorithm* in the sense of McConnell, Mehlhorn, Näher &
Schweitzer (2011).

**The witness ω:**
- the probe set (offsets, lengths, rules) and each probe's occurrence count or index range;
- the threshold t and the floor U;
- the list of chains whose bound exceeds U, with their band coordinates and DP maxima;
- the transcript of A*;
- for ties, two occurrence positions of ref(A*).

**Theorem V.** A checker that:
- (a) rescores A* from its CIGAR;
- (b) re-enumerates every probe's occurrence positions from an index of the reference, and rebuilds
  the chains;
- (c) recomputes U and every chain's bound from the probe layout;
- (d) recomputes the band maximum of every chain whose bound exceeds U;

and accepts only if every claim matches, never accepts a false interval, given a correct index for
(b).

*Proof.* Steps (a)–(d) re-establish exactly the hypotheses of Theorem C1. Rebuilding the chains in
(b) is what prevents a faulty mapper from omitting a chain: counts alone would not. ∎

**Trust base.** The checker is a few hundred lines with no heuristics or part selection. Its cost is
about that of CERTA's own enumeration and verification, so it can re-check every read, or a random
sample. A sample of n reads detects a certificate-error rate p with probability 1 − (1 − p)ⁿ. That
is a concrete QA procedure for clinical audits. The index can itself be validated against the
reference independently, once per reference build.

**Stretch goal: machine-checked proofs.** Lemma L, Proposition A and Theorem G are finite
combinatorial statements, suitable for mechanization in Lean 4.

### 4.8 Certified MAPQ

MAPQ in BWA-MEM's formula falls as SUB rises. A certified upper bound on SUB, obtained from the
second-best chain and U, therefore gives a certified *lower* bound on MAPQ.
- **Today:** the prototype reports BWA-style MAPQ from the SUB it found. That choice matched
  BWA-MEM2's variant-calling behaviour.
- **Under CERTA-X:** SUB becomes an interval, and Aim 3 measures how calibrated the conservative
  end is.
- **Known gap:** the pair tier's second-best score (FINDINGS §7) must also bound alternative
  partners inside the winning anchor's window.

## 5. The algorithm

**Tiers.** Each one is a refinement in the sense of §4.6, ordered by cost:

| Tier | Index | Probes / method | Status |
|---|---|---|---|
| **A** | hash, q = 22 (GPU) | 6 exact parts, Theorems 1/1′/2/3/4, L | [P+T], production prototype |
| **P** | A + window scan | mate anchors, ~12-base window parts (Theorem P) | [P+T] |
| **F** | second index, q ≈ 15 hash *or* FM-index | greedy optimal parts (Theorem G), absent parts, t = 2 | new [H] |
| **R** | FM-index counts | tie proofs (Theorem R), window collapse (Theorem R′), larger budgets | new [H] |
| **H** | F | one-mismatch probes (Lemma H), only where they lower U enough | new [S] |
| **X** | FM-index (suffix trie) | budgeted branch and bound (R5), search-scheme ordering | new; complete [P] |

**Per-read control flow:**
1. Tier A on the GPU, as today (about 93–95 % of reads end here, exact).
2. For pairs, run tier P.
3. For the rest:
   - **3a.** Compute the greedy τ-partition (tier F) and enumerate. If W > U, the read is exact.
   - **3b.** In repeats, prove ties by counting (tier R) and collapse identical windows.
   - **3c.** If the interval is still open and the budget allows, run tier H, then tier X.
4. Emit A*, [ℓ, u], the tie flag, the SUB interval and the witness. No read ever goes to a heuristic
   mapper.

**Implementation shortcut for tier F.** The prototype's index builder already accepts any q. A
second hash index at q = 15, s = 1 needs about 18 GB of 16-bit keys and positions, plus a 2 GB
directory. It gives
O(1) range lookups and positions directly, avoiding FM-index locate costs. Variable-length parts
(Theorem G proper) need an FM-index later; fixed ℓ = 15 tests the ceiling now. Two code limits must
be lifted first:
- `PMAX = 8`, because `pmask` is 8 bits;
- `KMAX = 5`.

## 6. Cost model and projections from our measurements

**Measured building blocks** (galaxy: Neoverse-N1 and A100):

| Component | Cost |
|---|---|
| Tier A, single-end, 20 M reads, A100 + 16 threads | 11.0 s wall (9× faster than BWA-MEM2's whole run) |
| `band_score` kernel (16-bit, vectorized), 95 × 150 band | 31 µs (~2.2 ns/cell) |
| `affine_align` (scalar, with traceback), same band | 99 µs |
| Pair pass (tier P) | ~0.46 ms CPU per attempted pair; 11.9 % of pairs attempted |
| CERTA's tiers A + P, paired, without fallback | 30 µs CPU per read |
| Reads in the minibwa fallback (paired) | ~0.32 ms CPU each |
| shmap-rs, 24 kb HiFi read, single thread | ~5.9 ms per read (`a2:~/shmap-rs/profiling/`) |

**Projected effect of tier F [H].** This uses the single-end HG002 loss distribution of §2.2, in which
5.34 % of reads remain uncertified after tier SL. The projection is approximate, and two effects
pull in opposite directions:
- BWA-MEM2's loss can only overstate the true optimum's loss, so more reads than counted here
  have a true loss under each ceiling;
- repeats can stop parts from fitting the budget, so fewer of them get certified.

| Loss ceiling | Uncertified reads below it | Change in exact share | Exact share (from 94.66 % after SL) |
|---|---|---|---|
| 25 (today, t = 2) | 34.0 % | — (these are mostly repeat-limited today) | 94.66 % |
| 45 (ℓ = 15, t = 2) | 57.4 % | about +1.25 points | about 95.9 % |
| 50 (ℓ = 15, t = 1) | 61.4 % | about +1.47 points | about 96.1 % |
| 60 (ℓ = 12, t = 1) | 67.8 % | about +1.80 points | about 96.5 % |

The extra points come from the reads with loss between today's ceiling (25) and the new one.

**Repeats (tier R).** The 3.6 % repeat-limited reads are the larger pool.
- **Evidence for budgets:** pass 3 (budget 4,096, t = 2) already certified +0.88 points, at 2.8×
  tier A's single-end map time. Window collapse (Theorem R′) targets that cost.
- **Evidence for tie proofs:** 83 % of the low-loss proper pairs in the fallback have an ambiguous
  mate. These are exactly the reads that Theorem R turns into proven ties.

**Projected endpoint [H].** About 97–98 % of reads exact. The remaining 2–3 % (chimeras,
contaminants, very high loss) get intervals, most of them narrow, or a proof that OPT < T. Aim 1
replaces these projections with measurements.

**Speed budget** (CPU per read, paired-end, from §2.2). The exact tiers F, R, H and X replace the
heuristic fallback and see only the 4.6 % uncertified reads:
- ≲ 1.2 ms CPU per residual read keeps CERTA's total CPU at or below minibwa's;
- ≲ 2.4 ms keeps it 3× below BWA-MEM2's.

Tier F at ℓ = 15 needs about 10 lookups and about 55 occurrences per read. Those are mostly
single-part chains, which t = 2 skips without DP, so it should fit comfortably [H]. Tier X has no
such guarantee, so it runs under a per-read budget and returns an interval when the budget runs out.

## 7. Evaluation plan and aims

| Aim | Deliverable | Measurement | Go / no-go |
|---|---|---|---|
| **1. Ceiling** (no new data needed) | q = 15 second index; PMAX/KMAX lifted; tier F pass | exact share, occurrence counts, CPU per residual read on HG002, HG001, HG005 | Exact share +1 point at ≤ 1 ms CPU per residual read |
| **2. Repeats** | FM-index counts; Theorem R ties; Theorem R′ collapse | exact share and tie proofs among the 3.6 % repeat-limited reads; DP calls saved | Half of the repeat-limited reads exact or proven ties |
| **3. Interval output and MAPQ** | [ℓ, u], SUB interval, conservative MAPQ, SAM tags (`XL`, `XU`) | interval-width distribution; MAPQ calibration on GIAB stratified regions | Variant F1 non-inferior to the current pipeline with **no** heuristic fallback |
| **4. Checker** | independent checker (Theorem V); sampling audit | check time per read; zero disagreements | Checker accepts 100 % of certificates; injected faults all detected |
| **5. Theory** | written proofs of Lemma H, Theorems C, G, R, V; optional Lean 4 mechanization | — | Proofs reviewed; mutation tests for each bound |
| **6. Long reads** | HiFi feasibility on a2's existing data (§8) | share of reads whose loss is below the certifiable ceiling | Go if ≥ 90 % of HiFi reads fit |

**Ground rules:**
- **Datasets:** the three GIAB samples already processed (HG001, HG002, HG005) and the existing
  protocol (chr20, GATK and bcftools, rtg vcfeval). New data only after approval, per
  RESEARCH_ROADMAP Part 4.
- **Tests:** every new bound gets the same treatment as Theorems L and P: a whole-reference brute-force
  oracle, adversarial reads placed at each threshold, and mutation testing. This also closes the
  open threshold-mutation gap (FINDINGS §7).

## 8. Long reads

**Bounds.** With P parts, the ceiling is about 5·P points of loss. A 15 kb HiFi read at
0.1–1 % error has 15–150 edits and about 680 disjoint 22-mers. If even a third of those parts are
unique, U leaves a margin of about 1,000 points against a loss of a few hundred. **Certifying the
locus is easy for HiFi** [H, Aim 6].
- ONT at Q20+ (1–2 % error) is borderline.
- Older noisy reads (5–15 %) are out of reach of exact-part certificates, by the same arithmetic.

**Verification cost.** This is the real problem. A band of ±1,000 diagonals over 15 kb is too
expensive. Within a locus, the optimal pairwise alignment should instead be computed by an exact,
near-linear aligner: A*PA/A*PA2 (Groot Koerkamp & Ivanov 2024), or WFA for small scores. Exact part
matches serve as its seeds. Certification then costs one exact pairwise alignment per candidate
locus, plus the bound.

**Objective.** Long reads cross structural variants, so the objective needs split or two-piece
affine alignment. A certified split alignment generalizes Theorem P (k pieces under a chaining
penalty). That is new theory, and it is listed as future work.

**Relation to shmap-rs.** shmap-rs maps in about 5.9 ms per HiFi read without base-level alignment.
A certifying long-read aligner will not beat it on raw speed. The useful combination: shmap-rs
places all reads, and CERTA-X certifies only the reads that shmap-rs flags as ambiguous (low
containment margin, segmental duplications). The proof is spent where mappers disagree.

## 9. Relation to prior work

| Approach | Exact? | Objective | Covers all reads? | Bounds when it stops? | Verifiable witness? |
|---|---|---|---|---|---|
| Heuristic mappers (BWA-MEM2, minibwa, strobealign, minimap2) | no | affine, heuristic | yes | no | no |
| Lossless k-error filters (RazerS 3, Yara, GEM, Columba; search schemes) | yes, within k edits | edit distance | only reads within k | no | no |
| BWT-SW (Lam et al. 2008) | yes | local alignment | yes, slowly | no | no |
| A*PA / A*PA2, WFA | yes | pairwise (global) | one pair at a time | no | no |
| CERTA prototype | yes, for certified reads | edit distance (S0/S1), affine local (SL), pairs (PR) | 93–95 % | no (heuristic fallback) | implicit |
| **CERTA-X (proposed)** | yes, or a proven interval | affine local and pairs, one objective | **yes** | **yes, [ℓ, u]** | **yes (Theorem V)** |

**Claimed novelty.** We are not aware of prior work that combines:
- (i) a per-read optimality certificate under an affine, clipping-aware objective, with
  repeat-adaptive part selection;
- (ii) a bounded interval answer when the per-read budget runs out;
- (iii) a checkable witness;
- (iv) Proposition A's ceiling, which explains why error-tolerant seeds give way to more exact
  parts under affine scoring.

The certificate idea itself parallels A*PA's seed heuristic, which is also an admissible "one cost
per unmatched seed" bound, applied here genome-wide instead of between two sequences.

## 10. Risks

- **Tier F cost on the real k-mer spectrum.** Real 15-mer counts are heavier-tailed than the model.
  Theorem G adapts part lengths, but the occurrence budget may bind earlier than projected. Aim 1
  measures this first.
- **Repeats.** Window collapse helps only for near-identical copies. Diverged multi-copy families
  may leave wide intervals. That is acceptable under the contract, but it limits the exact share.
- **The tie between edit distance and affine score.** Tiers S0/S1 certify edit distance and report
  an affine alignment at that locus (FINDINGS §6). CERTA-X unifies everything on the affine
  objective. Lemma L already supports this for the floor, but Aim 1 must confirm the cost.
- **Proof debt.** Lemma H and Theorem C need full write-ups, and the pair tier's second-best score
  needs its fix. The project's rule stands: no tier ships before its oracle and mutation tests.
- **Speed vs minibwa** stays the weakest claim. Exactness is the contribution, and speed is
  reported, not promised.

## 11. References

(Methods already reviewed are in [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md); only new or key citations
are listed here.)

- Backurs, A., & Indyk, P. (2015). Edit distance cannot be computed in strongly subquadratic time
  (unless SETH is false). *STOC 2015*, 51-58. https://doi.org/10.1145/2746539.2746612
- Cohen-Addad, V., Feuilloley, L., & Starikovskaya, T. (2019). Lower bounds for text indexing with
  mismatches and differences. *SODA 2019*, 1146-1164. https://doi.org/10.1137/1.9781611975482.70
- Groot Koerkamp, R., & Ivanov, P. (2024). Exact global alignment using A* with chaining seed
  heuristic and match pruning. *Bioinformatics*, 40(3), btae032.
  https://doi.org/10.1093/bioinformatics/btae032
- Groot Koerkamp, R. (2024). A*PA2: up to 19x faster exact global alignment. *WABI 2024*, LIPIcs
  312, 17. https://doi.org/10.4230/LIPIcs.WABI.2024.17
- Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018). Optimum search schemes
  for approximate string matching using bidirectional FM-index. arXiv:1711.02035.
- Kucherov, G., Salikhov, K., & Tsur, D. (2016). Approximate string matching using a bidirectional
  index. *Theoretical Computer Science*, 638, 145-158. https://doi.org/10.1016/j.tcs.2015.10.043
- Lam, T. W., Sung, W. K., Tam, S. L., Wong, C. K., & Yiu, S. M. (2008). Compressed indexing and
  local alignment of DNA. *Bioinformatics*, 24(6), 791-797.
  https://doi.org/10.1093/bioinformatics/btn032
- McConnell, R. M., Mehlhorn, K., Näher, S., & Schweitzer, P. (2011). Certifying algorithms.
  *Computer Science Review*, 5(2), 119-161. https://doi.org/10.1016/j.cosrev.2010.09.009
- Navarro, G. (2001). A guided tour to approximate string matching. *ACM Computing Surveys*, 33(1),
  31-88. https://doi.org/10.1145/375360.375365
- Marco-Sola, S., Moure, J. C., Moreto, M., & Espinosa, A. (2021). Fast gap-affine pairwise
  alignment using the wavefront algorithm. *Bioinformatics*, 37(4), 456-463.
  https://doi.org/10.1093/bioinformatics/btaa777
