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

**Lemma H, one-mismatch probes [P].** Let the parts have length ≥ 4, and let a part be witnessed
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

The full case analysis is written out in [certa/THEORY.md](certa/THEORY.md) §9 (Aim 5). It is
also checked exhaustively: every transcript against every part layout at L = 10, plus the densest
layout at L = 12, with no violation in 502 M pairs. With 3-base parts the bound fails, so the
length condition is needed.

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
sample. (Built and measured in Aim 4, §6.4: 540 lines, 89–134 µs per record.) A sample of n reads detects a certificate-error rate p with probability 1 − (1 − p)ⁿ. That
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
- **Measured (§6.3):** the conservative end is biased toward the reference and fails as a caller
  input. Lemma E costs one point per edit beyond the radius, a SNV five. It stays an audit tag.
- **Gap closed (Aim 5, §6.5):** the pair tier's second-best score now also bounds alternative
  partners inside the winning anchor's window, and placements of the mate elsewhere in its own
  chain (Proposition P2 in [certa/THEORY.md](certa/THEORY.md)).

## 5. The algorithm

**Tiers.** Each one is a refinement in the sense of §4.6, ordered by cost:

| Tier | Index | Probes / method | Status |
|---|---|---|---|
| **A** | hash, q = 22 (GPU) | 6 exact parts, Theorems 1/1′/2/3/4, L | [P+T], production prototype |
| **P** | A + window scan | mate anchors, ~12-base window parts (Theorem P) | [P+T] |
| **F** | second index, q ≈ 15 hash *or* FM-index | greedy optimal parts (Theorem G), absent parts, t = 2 | new [H] |
| **R** | hash index (s = 1); FM-index optional | long rare/absent parts by anchor extension (Theorem G), copy sampling, tie proofs (Theorem R) | measured (§6.2); window collapse (Theorem R′) not yet needed |
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

### 6.1 Aim 1 result (measured, October 2026)

**Setup.** A second index at q = 15 (s = 1): 2.95 × 10⁹ entries, 21.9 GB, built in under 2 minutes
at 84 GB peak memory. Two code limits were lifted:
- up to 16 parts per read, with a 4-bit part number in the candidate encoding;
- 16-bit part masks.

Local passes (Theorem L) may now use every part that fits the budget. The reads that the q = 22
pipeline leaves uncertified are re-run through the q = 15 index on the CPU (16 threads). The data
are random 2 M single-end samples of HG002 (NovaSeq X) and HG005 (NovaSeq 6000). The raw logs are
in [`certa/bench/results/aim1_2026-10/`](certa/bench/results/aim1_2026-10/summary.tsv).

| Sample | q = 22 default | + q15, t = 2 | q = 22 large budget | + q15, t = 2 | + q15, t = 2, large budget |
|---|---|---|---|---|---|
| HG002 | 94.66 % | 95.57 % (0.24 ms) | 95.65 % | **96.29 %** | 96.83 % (2.2 ms) |
| HG005 | 92.58 % | 93.98 % (0.25 ms) | 93.53 % | **94.66 %** | 95.31 % (1.9 ms) |

Times are CPU per residual read for the q = 15 pass.

**Go / no-go.** The criterion was +1 point at ≤ 1 ms CPU per residual read: **go**.
- *Cheap q = 15 pass alone* (t = 2, budget 256): +1.40 points on HG005 and +0.91 on HG002, at
  0.25 ms.
- *Combined with q = 22's large-budget pass:* +1.63 points on HG002 and +2.08 on HG005, at about
  0.6 ms per residual read. That pass's CPU cost is estimated from the increase in map time on 16
  threads.

**Correctness on real data.** The new certificates were checked against BWA-MEM2's alignments of the
same reads, as in FINDINGS §3: 39,649 checked on HG002 and 50,913 on HG005.
- 0 unexplained violations;
- every flagged case (30 and 33) sits in a window with reference `N`s, which BWA-MEM2 replaces with
  random bases.

**Against the projection.** §6 projected about +1.25 points from the default for ceiling 45. The
measurement gives +0.91 at budget 256 and +2.02 at budget 4,096, so the projection was in range.
One assumption was wrong: the random-sequence model of Proposition B predicted about 55 candidates
per read. The real genome gives about 200 at budget 256, four times more, because of its repeats.
The budget therefore binds earlier than the model says, as §10 anticipated, and this is where
tier R (counts and window collapse) and Theorem G's adaptive part lengths come in.

**Side finding: a regression fixed.** Since the local tier was added, the host's wider band (96
diagonals, for local and pair bands) had also widened the end-to-end certificate's
cluster-width limit. About 78 reads in 2 M, in tandem repeats, were then certified S1 on the CPU
but SL on the GPU: both valid proofs, but no longer byte-identical. A shared `BAND_E2E = 48`
restores identity; SAM bodies and uncertified reads are byte-identical again on 2 M HG002 reads.

### 6.2 Aim 2 result (measured, October 2026)

**What was built.** Tier R ([`certa/include/certa/repeat.h`](certa/include/certa/repeat.h); command
`certa repeat`). It implements Theorem G with the existing hash index (s = 1):
- **Parts:** a part [s, b) grows until it occurs at most τ times, or not at all.
- **Exact counts:** the occurrences of a long part are the hits of one anchor 22-mer inside it whose
  exact match extends to both ends of the part. One pass over the anchor's hits gives the count for
  every b.
- **Absent parts:** a part containing a 22-mer that occurs nowhere (one read error or SNP is
  enough) is absent and needs no enumeration.
- **The interval:** the upper bound comes from Lemma L plus branch and bound over chains. The lower
  bound comes from aligning the read to the copies that match it longest.
- **Ties:** Theorem R finds a second exact copy of the winning reference window.

**Tests.**
- *Oracle:* on the test genome (50-copy element, tandem repeat, diverged copies, decoys) the oracle
  checks 981 intervals per seed, on three seeds (ARM and x86). Every lo ≤ OPT ≤ hi; every exact
  claim equals OPT; every claimed copy is verified string by string.
- *Mutation testing:* all 6 injected faults are caught: a bound lowered by 3, a narrower band, a
  looser skip test, an off-by-one in hit extension, rare 22-mers treated as absent, and an
  overlapping copy counted as a tie.

**The pre-registered cohort** (§7, Aim 2). These are the reads the default pipeline leaves
uncertified while searching fewer parts than they have (|S| < P). They are followed through every
tier: q = 22 large budget, q = 15 (Aim 1) and tier R. Raw output is in
[`certa/bench/results/aim2_2026-10/`](certa/bench/results/aim2_2026-10/prereg_cohort.txt).

| Cohort | HG002 | HG005 |
|---|---|---|
| repeat-limited reads (share of all reads) | 70,908 (3.55 %) | 79,767 (3.99 %) |
| **exact or proven ambiguous**, tier R τ = 256 (0.5 ms/read) | **50.5 %** | **43.8 %** |
| the same, tier R τ = 1,024, larger limits (1.0–1.2 ms/read) | 53.2 % | 46.0 % |
| … of which near-perfect copies (BWA-MEM2 score ≥ 120) | 86.3 % | 82.9 % |
| … of which lossy reads (score < 120, or unmapped) | 3.0 % | 1.0 % |

Tier R's own share of the cohort, at τ = 256, is 14.8 points on HG002 and 13.0 on HG005. On both
samples, tier R proves ties for about 5–7 % of the cohort: two loci with the same optimal score,
hence MAPQ 0, proven.

**Validity on real data.** BWA-MEM2's alignments of the same reads never score above tier R's upper
bound, nor above the claimed optimum for exact reads: 0 violations among 180,804 residual reads.

**Go / no-go: partial go.** The criterion ("half of the repeat-limited reads exact or proven ties")
is met on HG002 and missed on HG005.
- **Where it fails:** not the repeat mechanism. On near-perfect copies, the reads the tier was
  designed for, it reaches 79–86 %.
- **What remains:** 40–45 % of the cohort are *lossy* reads (BWA-MEM2 score < 120: chimeras,
  low-complexity sequence, reads with many differences) that happen to fall in repeats. Exact
  certification of these needs excluding full-length alignments with 30–100 points of loss. That
  is beyond any exact-part certificate (Proposition A) and expensive for any exact method (§4.1).
- **The honest output for them:** an interval. The next lever is not a repeat tier but the
  interval contract itself (Aim 3): report [lo, hi] and a conservative MAPQ.

### 6.3 Aim 3: interval output (criterion pre-registered 2026-10-11, before the HG005 run)

**Deliverable.** `certa map --intervals` ([`certa/include/certa/interval.h`](certa/include/certa/interval.h)).
No read goes to a heuristic mapper. Every read that the certified tiers leave open goes through
tier R. Each SAM record carries:
- `XL`, `XU`: OPT ∈ [XL, XU];
- `XQ`: a certified lower bound of BWA-MEM's MAPQ formula, from a proven upper bound on SUB;
- `XT`: the tier (R exact, IV interval, MR placed by the mate's window, UN unmapped).

Pairs with an uncertified mate are placed from the alignments the tiers found plus an
exhaustive window search for each mate near the other. MAPQ follows BWA-MEM's pairing rule.

**Go criterion (fixed before the run).** HG005, chr20, the protocol of §7: GATK and bcftools,
SNVs and indels, rtg vcfeval against GIAB v4.2.1.
- **The criterion:** the interval pipeline's F1 is at least the current pipeline's
  (CERTA + minibwa fallback) minus 0.0010, in each of the 4 comparisons.
- **Policies:** the run is evaluated twice, from the same alignments. The criterion applies to
  MAPQ = the point estimate (BWA-MEM's formula with the best other locus found). MAPQ = XQ (fully
  certified) is reported alongside.
- **Data:** HG001's and HG002's FASTQs are no longer on disk, so HG005 is the only sample.

**Result (measured 2026-10-11): go.** F1 ([raw](certa/bench/results/aim3_2026-10/hg005_f1.tsv)):

| HG005 chr20 | CERTA + minibwa | **intervals, no fallback** | Δ | BWA-MEM2 | intervals, MAPQ = XQ |
|---|---|---|---|---|---|
| GATK SNVs | 0.9926 | **0.9933** | +0.0007 | 0.9917 | 0.3077 |
| GATK indels | 0.9938 | **0.9938** | 0.0000 | 0.9925 | 0.5051 |
| bcftools SNVs | 0.9926 | **0.9937** | +0.0011 | 0.9916 | 0.7849 |
| bcftools indels | 0.9472 | **0.9464** | −0.0008 | 0.9454 | 0.5842 |

- **All 4 comparisons are within the margin, and all 4 are above BWA-MEM2.** They still are after
  the Aim 5 fixes to the pair tier (re-run, §6.5).
  - GATK SNVs: 35 fewer false positives and 70 fewer false negatives than the current pipeline.
  - bcftools indels: the only loss, 20 more false negatives.
- **What the run produced** (725.8 M reads):
  - tier R on 7.14 % of reads;
  - records: 0.99 % exact by tier R, 3.87 % open intervals, 0.96 % placed in the mate's window,
    1.58 % unmapped;
  - 84.85 % of all records proven exact (XL = XU). XU − XL is 1–5 for 6.22 %, 6–10 for 1.76 %,
    11–25 for 2.19 %, 26–50 for 2.00 % and > 50 for 1.39 %.
- **The unmapped reads** are mostly low-quality or satellite reads without an error-free 22-mer.
  On HG002, BWA-MEM2 aligns about 60 % of them, mostly with score < 60 or MAPQ 0. The q = 15
  index of Aim 1 would seed them.
- **Cost:**
  - HG005: 2,647 s wall (64 threads + A100), against 2,409 s for CERTA + minibwa.
  - 2 M HG002 pairs: 343 CPU-s against 303 (+13 %), about 1.2 ms per uncertified read, which is the
    break-even budget of §2.2. BWA-MEM2 needs about 5× more.
- **Exactness by tier** (HG002, 2 M pairs):
  - S0, SR and tier R exact: 100 %;
  - PR mates: 70 %;
  - S1: 73 %. Edit-distance certificates are not affine certificates: another locus where a clipped
    alignment needs a few more edits can still score higher (Lemma E).

**MAPQ = XQ fails, and the reason is a finding.** The certified bound on SUB costs one point per
edit beyond the search radius, because a locus that matches all but the clipped last d bases
scores L − d − 5 (Lemma E). A SNV costs 5 points.
- **The consequence:** a one-SNV read certified at radius 2 (where most reads stop) gets XQ ≈ 10. At
  radius 5 it reaches XQ ≈ 21. Each further SNV needs about 5 more edits of radius.
- **Reads that match the reference** keep XQ 26–40.
- **The effect on callers:** a caller's MAPQ ≥ 20 filter keeps the reference allele and drops the
  variant (GATK SNV sensitivity 0.18).
- **The conclusion:** XQ is an audit quantity, not a caller input. It would become usable with a SUB
  bound from Lemma L instead of Lemma E. The clipped alignments that drive Lemma E's bound contain
  exact parts, so they are enumerated, and a SUB bound built the way tier SL builds its floor costs 5
  points per part (open).

**MAPQ calibration** (1 M simulated pairs with known origin, [raw](certa/bench/results/aim3_2026-10/sim_calibration.txt)).
Records misplaced with MAPQ ≥ 20:
- **Totals:** intervals 539, CERTA + minibwa 504, BWA-MEM2 31, minibwa 26.
- **Where they come from:** 460 of the 539 are certified reads that both CERTA pipelines share:
  - S0/S1 reads at radius 0–1, where `mapq_of` assumes no second locus;
  - PR pairs.
- **The new placements:** tier R and the window search add 79, against 3 for minibwa on the same
  reads.
- **XQ ≥ 20** covers 57 % of reads with 173 misplaced. All of them are PR pairs whose score-optimal
  placement is not the origin: 150 tandem shifts, and 23 segmental duplications where BWA-MEM2 picks
  the same pair. XQ makes the formula's inputs exact; it cannot make BWA-MEM's formula calibrated.

**Tests.**
- *Oracle:* every certified [XL, XU] and every SUB bound is checked against the whole-reference
  oracle, about 3,200 and 660 checks per seed, on 4 seeds, on ARM and x86. Tier R's SUB intervals
  and the pair placement are checked the same way.
- *Planted reads* make Lemma E and the own-chain searches tight.
- *Mutation testing:* 11 of 12 injected faults are caught. The one that survives weakens a term no
  alignment can reach for reads ≥ 100 bp.
- *Regression:* default output is byte-identical to the previous version.

### 6.4 Aim 4 result: the independent checker (measured, October 2026)

**What was built.**
- **The checker:** `certa check` ([`certa/src/check.cpp`](certa/src/check.cpp), about 540 lines). It shares
  no alignment, bound or lookup code with the mapper: it has its own index lookup, Lemma L, banded
  DP and SAM parser.
- **The witness:** `certa map --intervals --witness` adds what the checker needs to each record:
  - `XW`: the probes (read intervals) and the support t that the bounds were proven from;
  - `XH`: the SUB bound behind XQ;
  - `XV`: the alignment attaining XL, when it is not the record's own.
- **The verification**, per record (Theorem V, in its strongest form: the checker rebuilds the chains
  itself instead of trusting a list):
  1. It rescores AS from the CIGAR.
  2. It checks that an alignment scoring XL is shown.
  3. For XU, it recomputes U_t (Lemma L), enumerates every exact occurrence of every probe, chains
     them with pad = L − XU − 7, and runs a DP over every chain whose Lemma L bound exceeds XU.
  4. It checks XH the same way, counting only loci that don't overlap the record's.
  5. It checks XQ against BWA-MEM's formula at (AS, XH).
- **The index:** `certa check-index` proves, once per reference build, that the index holds exactly
  the N-free q-mers at every s-th position (all keys match, no position repeats, the count is right)
  and that the stored reference equals the FASTA. For GRCh38 (2.95 G entries) this takes 34 s on 64
  threads.

**Results** ([raw](certa/bench/results/aim4_2026-10/)):

| Test | Result |
|---|---|
| Honest records, HG002 2 M pairs | **4,000,000 of 4,000,000 accepted**: 920 k non-trivial XU and 3.28 M XH bounds proven |
| Honest records, HG005 2 M pairs | **4,000,000 of 4,000,000 accepted**: 1.16 M XU and 3.11 M XH bounds proven |
| Faults false by construction, real SAM | **23,733 of 23,733 rejected**: AS + 1, last base clipped, XU < XL, XQ + 10, witness removed |
| Faults of unknown truth, real SAM | moved alignment 99.87 %, XH − 10 99.81 %, moved probe 57 % rejected (the rest stay true and provable) |
| Untouched records rejected | **0** |
| Mapper mutants on real data (6) | **all 6 flagged**, 20 to 143,329 records each. One is a bug in the mapper's own Lemma L, caught by the checker's independent one. |
| Oracle unit tests (4 seeds, ARM and x86) | every honest claim accepted; every claim falsified against the oracles rejected (7,000–8,500 per seed, a third of them caught by the bound search itself) |
| Corrupted indexes and FASTA | all rejected |

**Cost.** 89 µs CPU per record on HG002 and 134 µs on HG005, about the mapper's own 86 µs. That is
11.5 s on 32 threads for 4 M records.
- **Sampling audits:** a uniform sample of 100,000 records takes 0.6 s. When it is all accepted, it
  excludes an error rate above 0.003 % at 95 % confidence.
- **A real bug found:** the checker's first real-data run rejected 2 of 4 M records. Their
  alignments ran off the end of the circular chrM, and the mapper had dropped the alignment but kept
  XL and XQ. The fix makes such records claim nothing they cannot show.

**Go / no-go: go, with one declared gap.** Every claim the checker covers is verified, and every
injected fault that falsifies one is detected.
- **The gap:** the XQ of PR records (about 4 % of records) rests on the pair certificate (Theorem P),
  which the checker does not re-derive yet. These are reported as unchecked; their XL and XU are
  checked.

**Update 2026-10-14: the gap is closed.** `certa check` now verifies the XQ of PR mates as well
(Theorem V for PR mates, [certa/THEORY.md](certa/THEORY.md) §8;
[raw](certa/bench/results/aim4_pr_2026-10/)).
- **What a PR record now carries** (`--witness`):
  - `XY`: the proper-pair distance D and the pair certificate's own parts for each mate
    (support 1);
  - `XH` for the single-end term of its XQ, with `XV` showing the single-end alignment it is
    about.
- **How the checker proves XS,** the bound on every other proper pair, from the pair's two
  records:
  - pairs with no exact part in either mate are bounded by Lemma L per mate;
  - for the rest, each chain of a mate's part occurrences is paired, by an exhaustive DP, with every
    placement of the partner within D, excluding the shown pair.
- **A first version was too coarse.** It added the best x in a chain to the best partner anywhere in
  its window, and left 24 honest HG002 records unproven. Exact local oracles confirmed that the
  mapper's S2 held in every case.
  - *The fix:* pairing per start. Every pair counted is then a real proper pair, so an honest XS is
    always provable, up to work caps.

| Test | Result |
|---|---|
| Honest records, HG002 2 M pairs | **4,000,000 accepted**; XS proven for all 148,503 PR mates whose XQ rests on it |
| Honest records, HG005 2 M pairs | **4,000,000 accepted**; XS proven for all 155,179 |
| Faults false by construction, real HG002 SAM | **24,306 of 24,306 rejected**: the Aim 4 faults, plus PR XQ + 10 (261), pair witness removed (442), D below the pair's distance (421) |
| Untouched records rejected | **0** |
| Oracle tests, 4 seeds | every certified pair accepted. All false pair claims rejected: 520–637 per seed in the PR tests, 113–123 more in the threshold families. They are XS one below the oracle's best other pair (ties included), witness removed, D too small, XQ one above the rule, a mate moved |
| Mutants of the pair check (C10–C19) | **all 10 killed**. Families G–J and a rebuilt H each isolate one bound. With them the suite has 66 mutants: 60 killed, 6 proven equivalent (§6.5) |

**Cost:** 181–232 µs CPU per record, up from 89–134 µs; the pair proofs account for the
difference. Checking 4 M records takes 23–30 s on 32 threads.

**Speed budget** (CPU per read, paired-end, from §2.2). The exact tiers F, R, H and X replace the
heuristic fallback and see only the 4.6 % uncertified reads:
- ≲ 1.2 ms CPU per residual read keeps CERTA's total CPU at or below minibwa's;
- ≲ 2.4 ms keeps it 3× below BWA-MEM2's.

Tier F at ℓ = 15 needs about 10 lookups and about 55 occurrences per read. Those are mostly
single-part chains, which t = 2 skips without DP, so it should fit comfortably [H]. Tier X has no
such guarantee, so it runs under a per-read budget and returns an interval when the budget runs out.

### 6.5 Aim 5 result: written proofs and mutation tests (October 2026)

**What was built.**
- **Written proofs:** [`certa/THEORY.md`](certa/THEORY.md) states and proves every bound the
  implementation uses. Each statement names the code, the oracle tests, the mutants that break it
  and, for finite statements, the exhaustive check. It covers:
  - the index (Lemma I, Proposition V0) and the end-to-end certificates (Theorems 1–4);
  - the affine local objective (Lemmas L, S, W, C, E, O; Theorem L);
  - pairs (Lemmas P0 and WS, Theorem P, Proposition P2 for S2);
  - repeats (Theorem G, Lemma R, Theorems R-int and R);
  - XQ (Proposition Q) and the checker (Theorem V, as built);
  - the bounds proposed but not built: Lemma H in full (it was a sketch, §4.3), Proposition A,
    Theorem C.
- **Exhaustive checks** of the finite statements
  ([`certa/tests/verify_lemmas.cpp`](certa/tests/verify_lemmas.cpp), in CI, 33 s).
- **Mutation tests** ([`certa/tests/mutants.tsv`](certa/tests/mutants.tsv),
  `tests/run_mutants.sh`): 56 one-line faults. Each weakens one bound by 1–3 points or drops one
  case. Every mutant runs through the oracle tests on two seeds, and it is killed when a check
  fails.

**Results** ([raw](certa/bench/results/aim5_2026-10/)):

| Test | Result |
|---|---|
| Lemma L, parts ≥ 3: every transcript against every part layout, L = 9 | 52.8 M pairs, **0 violations**. With 2-base parts it fails (17 counterexamples). |
| Lemma H, parts ≥ 4: all layouts at L = 10, the densest at L = 12 | 502 M pairs, **0 violations**. With 3-base parts it fails. |
| Lemma E (L = 8), Lemma S, Proposition A, Theorem G, XQ | all hold. Lemma E is attained for d = 1–4; greedy = DP optimum on 20,000 instances. |
| Mutants (end-to-end 7, SL 8, PR and S2 14, tier R 6, intervals 12, checker 9) | **50 killed; the other 6 proven equivalent** |
| Oracle tests, 4 seeds | all pass |

**The 6 equivalent mutants** change no claim on any input, so no test can kill them. THEORY.md §9a
proves each one:
- **L6 and A2** weaken a term that is never the maximum (Lemma S), or never binds for reads
  ≥ 98 bp (Lemma O).
- **P2, P4, P8 and P13** sit where Theorem P is redundant:
  - a pair at the band edge is also bounded through its partner's anchor (Proposition P3);
  - the union window already covers the scan pad;
  - S2's window term implies the `relevant` term.

**What writing the proofs found: two bugs in tier PR's S2, both fixed.**
1. **The second partner** (the gap declared in §4.8). S2 bounded neither a second partner inside
   the winning anchor's window nor a placement of the mate elsewhere in its own chain (tandem
   repeats). A new test oracle for the best other pair found 9 such pairs. S2 now searches both,
   excluding R1's own alignments (Proposition P2).
2. **The window pad.** For a high threshold, the window search's band pad went negative, and an
   exact partner was missed. The same oracle found it. The pad is now clamped at 0.

Killing every mutant also forced new adversarial inputs:
- pairs placed at each PR threshold (families A–F), which closes the threshold-mutation gap of
  FINDINGS §7;
- planted reads at tier SL's thresholds and at Theorem G's read end;
- a direct test of Lemma WS;
- two checker tampers. The checker now also rejects a record whose own alignment scores above its
  XU.

**Effect on real data.**
- **2 M HG002 pairs** ([raw](certa/bench/results/aim5_2026-10/pr_fix_effect_hg002.txt)): 10,901
  records (0.29 %) change MAPQ, almost all of them PR. 4,419 fall below 20, 23 rise to 20 or above,
  and 751 placements change.
- **The Aim 3 criterion still holds** on the fixed code. HG005 was re-run with the same protocol
  and pre-registered MAPQ policy ([raw](certa/bench/results/aim5_2026-10/hg005_rerun_f1.tsv)):

| HG005 chr20 | CERTA + minibwa | intervals, Aim 3 | **intervals, Aim 5** | Δ vs CERTA + minibwa | BWA-MEM2 |
|---|---|---|---|---|---|
| GATK SNVs | 0.9926 | 0.9933 | **0.9933** | +0.0007 | 0.9917 |
| GATK indels | 0.9938 | 0.9938 | **0.9935** | −0.0002 | 0.9925 |
| bcftools SNVs | 0.9926 | 0.9937 | **0.9938** | +0.0012 | 0.9916 |
| bcftools indels | 0.9472 | 0.9464 | **0.9467** | −0.0005 | 0.9454 |

All 4 are within the 0.0010 margin and above BWA-MEM2. Records with MAPQ ≥ 20 fall by 1.10 M
(0.15 % of 725.8 M reads), and wall time rises by 0.5 %.

**Go / no-go: half met.**
- **Mutation tests for each bound: met.** Every bound the implementation uses has a mutant, and
  every mutant is killed or proven equivalent.
- **Proofs reviewed: not yet.** THEORY.md has not been reviewed by a second person; its §11 lists
  where to start.
- **Lean 4 (optional): started (2026-10-15).** Lemmas L, E, H and S and Theorem G are
  machine-checked for reads of any length ([certa/lean](certa/lean/README.md); Lean 4.35, no
  Mathlib), with the counterexamples showing the part lengths are needed (2 for Lemma L, 3 for
  Lemma H, 4 for Lemma S), and Lemma E's bound attained for every d. Lemma S makes mutant L6
  equivalent by proof. Theorem G is proved for F_τ with occurrences counted in a reference text.
  The band pruning that brought HiFi certification under minimap2's cost is proved for alignment
  paths, as are Lemma W (the pad), Lemma C's diagonal geometry and the vectorized DP row's
  identity. Theorem L-J's pruning remains.

### 6.6 Aim 6: long reads (criterion pre-registered 2026-10-13, before the run)

**Question.** Does the exact-part certificate reach HiFi reads? §8 argues that it does: a 13 kb
read has hundreds of disjoint parts, so the ceiling (about 5 points of loss per part) lies far above
the loss of a HiFi alignment. This aim measures that on real reads before anything is built.

**Data.** The 242,534 real HG002 HiFi reads in `hifi_1x.fa` (GIAB PacBio CCS 15 kb, Sequel
Q20 reads; all reads of cell m54238_180901_011437 and the first 103,846 of m54238_180902_013549,
about 1× coverage; md5 `18fec02a…`), already on a2 and galaxy. Reference: GRCh38, the same FASTA
as all CERTA runs. No new data.

**Definitions (fixed before the run).**
- **Best known alignment:** minimap2 2.31 `-x map-hifi -c --eqx`, primary alignment. Its score
  under CERTA's scores (match +1, mismatch −4, gap 6 + k, 5 per clipped end) is σ ≤ OPT.
- **Parts:** for each strand of the read, disjoint parts that each occur at most τ = 16 times in
  GRCh38. They are chosen greedily left to right (Theorem G) and counted from the q = 22, s = 1
  index by anchor extension (Lemma R). Absent parts (no occurrence) are included. Stretches where
  no part of ≤ τ occurrences starts are skipped.
- **Ceiling:** U₁, the Lemma L bound with t = 1 over those parts: an alignment with no exact part
  occurrence scores at most U₁. The larger value over the two strands is used.
- **A read fits** when σ > U₁. Its optimum then lies among the alignments through enumerated
  occurrences, so certifying it needs no heuristic step. What that costs is measured separately.

**Go criterion.** At least 90 % of all 242,534 reads fit. Unmapped reads count as not fitting.

**Also reported, not part of the criterion:**
- the margin σ − U₁;
- why reads miss: unmapped, clipped, or repeat-limited;
- the enumerated occurrences per read, as a proxy for the verification cost;
- the same at τ = 256.

The measuring tool is debugged on 500 reads of the 10× set that are not in the sample.

**Result (measured 2026-10-13): go.** Tool: [`certa/tools/lr_ceiling.cpp`](certa/tools/lr_ceiling.cpp);
[raw](certa/bench/results/aim6_2026-10/).

| HG002 HiFi, 242,534 reads | **τ = 16 (pre-registered)** | τ = 256 |
|---|---|---|
| **Reads that fit (σ > U₁)** | **226,253 = 93.29 %** (95 % CI 93.19–93.39 %) | 227,364 = 93.75 % |
| Median loss L − σ / median ceiling L − U₁ | 247 / 2,605 | 247 / 2,605 |
| Median margin σ − U₁ (mapped reads) | 2,289 | 2,291 |
| Parts on the alignment strand, median | 521 (41 per kb) | 521 |
| Enumerated occurrences per read: median, 99th percentile | 653, 2,425 | 1,295, 13,645 |

- **The checks:** every part's occurrence count was recounted from a second anchor, and every U₁
  recomputed by brute force: all equal. Every primary CIGAR's match and mismatch claims were
  checked against the reference: none failed.
- **§8's estimate holds: certifying the locus is easy.** A typical read's loss (median 247 points
  for a 12.8 kb read) is a tenth of its ceiling (median 2,605).
- **The 6.7 % that miss** fall into the classes §8 anticipated (descriptive classes, defined in
  `bench/aim6_summary.py`):
  - *clipped*, 2.32 % of reads: 80 % have split alignments (SV breakpoints, chimeras);
  - *high loss*, 2.06 %: 59 % contain a gap of ≥ 500 bases. That is an SV inside the alignment,
    which a single affine alignment charges about 2 points per inserted base;
  - *repeat-limited*, 1.81 %: median MAPQ 1. Long parts in multi-copy repeats lower the ceiling;
  - *unmapped*, 0.53 %.
  - By region: half of the chrY reads miss, and nearly all reads on the satellite contig
    KI270438.1.
- **What it does not show.** "Fit" is the bound half of a certificate. A fitting read still needs
  its candidate chains verified: the true locus by one exact long-read alignment (A*PA or WFA, §8),
  and every other chain whose Lemma L bound exceeds OPT.
  - Enumeration is cheap: a median of 653 occurrences per read, mostly one per part at the true
    locus.
  - The alignment cost is measured in the follow-up below.
  - σ comes from minimap2, so the share of misses is an upper estimate: OPT can only be higher.
- **Next:** the split-alignment objective of §8 addresses the clipped and high-loss classes
  (4.4 % of reads).

**Follow-up (2026-10-14): the verification cost, and the complete certificate.**
[raw](certa/bench/results/aim6_2026-10/certify.txt). The bound half needed no exact long-read
aligner, and neither does the rest.
- **The idea:** a read that fits has a known alignment scoring lb > U₁. Any alignment scoring more
  has gap total ≤ L − lb − 7 (Lemma W), about 250 for a median read. So the band is ±250, not
  ±2,600.
- **The method:** a banded DP over every chain of part occurrences whose Lemma L bound exceeds the
  best so far proves OPT exactly (Corollary L-long in [certa/THEORY.md](certa/THEORY.md); tool
  `lr_ceiling --certify`).

| HG002 HiFi, 242,534 reads | Result |
|---|---|
| **OPT proven** | **226,330 reads (93.32 %)**: every read that fits (with the trimmed path) |
| Chains needing a DP per read | median 1 (the true locus), 99th percentile 2, max 83 |
| DP cells per read | median 2.1 M computed, of a 6.1 M band |
| DP time per read (one thread) | median 3.1 ms, mean 5.0 ms, 99th percentile 30 ms (the first, scalar DP: 20, 35, 217 ms) |
| Whole run, CPU | 3,575 s, about 15 ms per read (the first version: 11,629 s). minimap2's mapping and alignment take 4,171 s: **certifying costs less than minimap2 (0.86×)**; the first version cost 2.8× |
| OPT above minimap2's alignment (CERTA's scores) | 17.8 % of certified reads, median by 1 point. 333 reads by ≥ 50 against every minimap2 line, e.g. an optimum through bases minimap2 clipped |
| Check on 2,005 reads (every 121st) | parts recounted, U₁ by brute force, and a full DP over ±1 kb around minimap2's locus: never above the proven OPT, equal to it for 1,876 of 1,879 |

**From 2.8× to 0.86× minimap2's CPU (2026-10-15)**
([raw](certa/bench/results/aim6_2026-10/certify_cost.txt),
[kernel](certa/bench/results/aim6_2026-10/vector_dp.txt)). Three changes; none changes a
certificate.
- *The DP row is vectorized* (`band_row`). The deletion chain of a Gotoh row is a prefix
  maximum, computed with a few in-register shifts per vector. The kernel alone is 2.3× faster on
  galaxy's NEON (4 lanes) and 8× on x86 with AVX2 (8 lanes). This took the cost to 1.6×.
- *Batched index lookups.* A profile then showed half the CPU in partitioning the read into
  parts, mostly in the index lookups: two dependent random reads into a 23 GB index per q-mer.
  They are now pipelined in groups, with prefetches.
- *The DP pruned to cells that can beat lb* (THEORY.md, "Pruning the band"). A cell's value plus
  the bases after it bounds every alignment through it. Only 35.5 % of the band's cells pass, and
  35 % are computed. The 3 reads that used to hit the cell cap now fit under it (OPT = lb).

**Follow-up 2 (2026-10-14): reads that cross structural variants, bound side**
([raw](certa/bench/results/aim6_2026-10/split.txt); exploratory, not pre-registered).
- **The objective:** the split (jump) objective lets a read align as pieces in read order, each on
  its own locus and strand, at a cost P per jump ([certa/THEORY.md](certa/THEORY.md), "Jump (split)
  alignments"). P = 10 matches how BWA-MEM scores a split read.
- **Lemma L-J:** it bounds every split alignment with no exact part by U_J, which one pass over the
  read computes. It equals a brute force on 210,000 random cases.
- **The lower bound:** lb_J is assembled from minimap2's own lines. Each line is also split at gaps
  longer than P − 6.

| Share of the 242,534 reads that fit | P = 10 | P = 30 | P = 60 |
|---|---|---|---|
| Single piece (the criterion, as above) | 93.29 % | 93.29 % | 93.29 % |
| **Split objective** | **94.58 %** | **94.39 %** | **94.30 %** |
| Clipped misses that now fit | 1,419 of 5,627 | 1,307 | 1,211 |
| High-loss misses that now fit | 1,005 of 4,987 | 837 | 753 |

- **The misses that stay:** unmapped reads, and reads in centromeric satellites, chrY and
  diverged duplications. There minimap2's identity is 0.87–0.95, so the reference does not hold the
  read closely enough for any exact-part certificate.
- **A long novel insertion can never be certified exactly by exact parts.** It is a piece that no
  part can witness, so such a read keeps an interval.
- **The exact split certificate (follow-up 3, 2026-10-14/15)**
  ([raw](certa/bench/results/aim6_2026-10/split_certify.txt); Theorem L-J in
  [certa/THEORY.md](certa/THEORY.md)).
  - *How it works:* each evaluated chain gets an exact banded DP in which a piece may start from the
    best prefix value minus P. Reverse-strand chains run against the reverse-complemented reference.
    An upper layer adds Lemma L bounds for pieces that no part witnesses.
  - *Which chains:* they are evaluated in rounds, pruned against exact prefix and suffix values of
    the chains evaluated before (the suffix from a pass over the reversed read).
    - Each chain's DP runs only on a read window where its pieces can matter.
    - lb is raised to the best alignment found, which shrinks the band.
    - When the rounds would cost more, one static step prunes against the first round's values.
    - A read over the cell cap keeps an interval.
  - *Tested:* against an exact brute force of the split objective on a synthetic genome with
    planted SVs. Every proven optimum equals it and every interval contains it, also from the
    weakest and the tightest admissible lb, with one chain per round, on the static path and under
    a small cap. The pruning values are checked against brute force on 20,000 random layouts.
    Mutants LR1–LR27 are all killed.
  - *On the 2,674 HiFi reads that fit only the split objective* (the first version, with static
    pruning only, in brackets):

    | | Result |
    |---|---|
    | **OPT_J proven** | **1,551** (1,544): 0.64 % of all reads, so exact certification covers 93.96 % of reads. Same OPT_J on every read both prove |
    | An interval | 1,043 (991), median width 112; diverged or novel insertions keep intervals by design |
    | Over the cell cap | 80 (139), each still with an interval |
    | Time per read (one thread) | median **1.4 s** (45 s), mean 19 s (62 s), with the vectorized DP row |
    | CPU in all | 66,783 s (170,706 s) |
    | Chains evaluated | median 48 of 154 (145 of 148 kept) |

  - *On 500 reads (all 476 that fit the split objective):* OPT_J proven for 442, a median of
    0.14 s per read, one chain evaluated.
  - *What stays costly:* reads with a long lossy stretch, such as an insertion of repeat sequence.
    While a chain covering the stretch is unevaluated, the bound lets the prefix cross it almost
    free, so many chains look viable. The unwitnessed bound (about 0.7 per base) then leaves an
    interval about 0.7 × the stretch wide.

## 7. Evaluation plan and aims

| Aim | Deliverable | Measurement | Go / no-go |
|---|---|---|---|
| **1. Ceiling** (no new data needed) | q = 15 second index; PMAX lifted to 16; tier F pass | exact share, occurrence counts, CPU per residual read on HG002, HG001, HG005 | Exact share +1 point at ≤ 1 ms CPU per residual read. **Done for HG002 and HG005: go (§6.1).** HG001 needs its FASTQs re-downloaded. |
| **2. Repeats** | FM-index counts; Theorem R ties; Theorem R′ collapse | exact share and tie proofs among the 3.6 % repeat-limited reads; DP calls saved | Half of the repeat-limited reads exact or proven ties. **Partial go (§6.2): 50.5–53.2 % on HG002, 43.8–46.0 % on HG005; 79–86 % of near-perfect copies.** Built with the hash index (no FM-index needed); window collapse (R′) not needed so far. |
| **3. Interval output and MAPQ** | [ℓ, u], SUB interval, conservative MAPQ, SAM tags (`XL`, `XU`) | interval-width distribution; MAPQ calibration on GIAB stratified regions | Variant F1 non-inferior to the current pipeline with **no** heuristic fallback. **Go (§6.3): HG005 F1 −0.0008 to +0.0011 against CERTA + minibwa, above BWA-MEM2 in all 4.** Calibration was measured on simulated pairs instead (the GIAB stratification files are not downloaded). XQ as MAPQ fails: it is biased toward the reference. |
| **4. Checker** | independent checker (Theorem V); sampling audit | check time per read; zero disagreements | Checker accepts 100 % of certificates; injected faults all detected. **Go (§6.4): 8 M real records accepted; all 23,733 false-by-construction faults and all 6 mapper mutants detected; 89–134 µs per record.** The gap (pair-level XQ, about 4 % of records) **is closed (§6.4, 2026-10-14): all 303,682 XS claims on 8 M records proven; 181–232 µs per record.** |
| **5. Theory** | written proofs of Lemma H, Theorems C, G, R, V; optional Lean 4 mechanization | — | Proofs reviewed; mutation tests for each bound. **Half met (§6.5): every bound's mutants killed (50) or proven equivalent (6); exhaustive lemma checks pass; two S2 bugs found and fixed.** Review by a second person pending; Lean started: Lemmas L, E, H, S and Theorem G machine-checked. |
| **6. Long reads** | HiFi feasibility on a2's existing data (§8) | share of reads whose loss is below the certifiable ceiling | Go if ≥ 90 % of HiFi reads fit. **Go (§6.6, definitions pre-registered): 93.29 % of 242,534 HG002 HiFi reads fit (93.75 % at τ = 256).** The misses span SVs (4.4 %), sit in repeats (1.8 %) or are unmapped (0.5 %). **Follow-up: the certificate completed for 93.32 % (OPT proven), at 0.86× minimap2's CPU; with the split objective, 93.96 % exact (§6.6).** |

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
locus is easy for HiFi.** Measured in Aim 6 (§6.6): 93.3 % of real HG002 HiFi reads fit, with a
median margin of 2,289 points. Most of the rest span structural variants.
- ONT at Q20+ (1–2 % error) is borderline.
- Older noisy reads (5–15 %) are out of reach of exact-part certificates, by the same arithmetic.

**Verification cost.** This looked like the real problem: a band of ±1,000 diagonals over 15 kb is
too expensive.
- **The proposed route:** compute the optimal pairwise alignment within a locus with an exact,
  near-linear aligner (A*PA/A*PA2, Groot Koerkamp & Ivanov 2024, or WFA for small scores), seeded
  by the exact part matches.
- **The measured route (§6.6):** it turned out simpler. A known alignment of score lb shrinks the
  band to L − lb − 7, about ±250 for a median HiFi read. Usually one chain needs it. Certification
  then costs a median of 3.1 ms of (vectorized, pruned) DP per read, and 0.86× minimap2's CPU in
  all.
- **Where an exact near-linear aligner would still help:** in making that DP cheaper.

**Objective.** Long reads cross structural variants, so the objective needs split or two-piece
affine alignment. A certified split alignment generalizes Theorem P (k pieces under a chaining
penalty).
- **Defined and bounded (§6.6):** the jump objective, with pieces in read order and P per jump.
  Lemma L-J extends the exact-part ceiling to it, computed in one pass and verified exhaustively.
- **Measured:** on HiFi it lets 1.1 more points of reads fit at P = 30, about a fifth of the
  clipped and high-loss reads.
- **The exact certificate of the split optimum:** built and oracle-tested. It is proven for 1,551
  more HiFi reads, at a median of 1.4 s per SV read. The mean, 19 s, comes from reads with long
  lossy stretches (§6.6).

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
