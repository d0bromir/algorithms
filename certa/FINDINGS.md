# CERTA: weaknesses of current short-read aligners, and a certified alternative

This is a paper-oriented summary of what this repository has established:
what is measured, what is proven, and what is not yet claimed. Raw data
are in [`bench/results/hg002_2026-09/`](bench/results/hg002_2026-09/); the
methods are in [`bench/RESULTS.md`](bench/RESULTS.md) and the code in
[`include/certa/core.h`](include/certa/core.h).

> **Status:** single-end reads, one GIAB sample (HG002, NovaSeq X), chr20
> for variant accuracy. This is a pilot, not the multi-sample protocol of
> `RESEARCH_ROADMAP.md` Part 4. Accuracy and timing tables are marked with
> the version they were measured on.

## 1. The problem with current aligners

All production short-read aligners (BWA-MEM2, minibwa, strobealign,
minimap2, bowtie2) are **heuristic**. They search a limited set of seeds
and chains, then extend with banded DP. They are usually right, but:

- **W1. No per-read guarantee.** A read can be placed at a locus that is
  worse than another one in the genome, and nothing in the output says so.
- **W2. Confidence is a heuristic too.** MAPQ comes from the hits the
  seeding happened to find. A missed paralog gives a confident MAPQ to an
  ambiguous read, and a real tie can be reported as unique.
- **W3. Speed is bought with search effort.** The faster tools search less,
  and on HG002 they make more provable errors (below). There is no knob that
  is safe by construction.
- **W4. Reference `N` handling.** BWA-family indexes replace `N`/IUPAC
  reference bases with random nucleotides. A read can then "match" a base
  the reference does not contain (observed on HG002, section 4).

### Measured: provable errors of each aligner on real data

CERTA certifies the optimal locus of a read, and whether it is unique or
tied, by proof (section 2). That makes it possible to count, for real reads,
where each other aligner is *provably* wrong. Scores are recomputed from
CIGAR and NM under BWA-MEM's default scheme for every tool, so each tool is
judged by the same scale. The data are 2 M HG002 NovaSeq X reads; 1.89 M are
certified.

| Tool | Placed at a locus ≥ 1 edit worse than the proven best | …of which with MAPQ ≥ 20 | Unmapped although a ≤ 5-edit alignment exists | Exact multi-copy reads given MAPQ ≥ 20 |
|---|---|---|---|---|
| BWA-MEM2 2.3 | 177 (0.009 %) | 13 | 0 | 0 |
| minibwa 0.7 | 420 (0.022 %) | 38 | 0 | 0 |
| strobealign 0.18 | 4,984 (0.264 %) | 800 | 11 | 704 (1.01 % of ties) |
| minimap2 2.31 `-x sr` | 2,237 (0.118 %) | 54 | 524 | 0 |
| bowtie2 2.5.5 | 921 (0.049 %) | 61 | 134 | 0 |
| **CERTA, certified reads** | **0 by construction** | **0** | **0** | **0** |

(`bench/weaknesses.py`, `results/hg002_2026-09/weaknesses_2M.txt`.)

- **Per genome:** BWA-MEM2's rate means about 40,000 provably misplaced
  reads per 22× genome. strobealign, one of the fastest tools, has 28× that
  rate, and 1 % of reads that exist in ≥ 2 identical copies get MAPQ ≥ 20.
- **Speed vs errors:** there is a clear speed/accuracy trade-off among the
  heuristic tools (timings in `bench/RESULTS.md` §2): minibwa is 3.8× faster
  than BWA-MEM2 with 2.4× its provable errors.

## 2. The algorithm and what it guarantees

CERTA maps a read in two tiers: **certified** and **fallback**. The
certified tier makes the claims below. Reads it cannot certify go to a
conventional aligner and carry no claim.

Notation:
- **r:** a read of length L.
- **d(r, ℓ):** the unit-cost semi-global edit distance of r at reference
  locus ℓ. An `N` in the read or the reference always costs 1.
- **Parts:** the read is split into P disjoint parts, each of length at
  least q + s − 1 (q = 22).
- **Index:** it stores the q-mer at every s-th reference position (s = 1 or
  8).

**Theorem 1 (lossless enumeration).**
- *Statement:* let S be a set of parts. Suppose every index hit of every
  part in S is enumerated, at all s shifts. Then every locus ℓ with
  d(r, ℓ) ≤ |S| − 1 yields a candidate diagonal within ±d(r, ℓ) of ℓ's start
  diagonal.
- *Proof:* d edits touch at most d parts. Because |S| > d, some part of S is
  untouched and occurs exactly at ℓ. Its diagonal differs from ℓ's start
  diagonal by the number of indels before it, which is at most d. Any exact
  occurrence of length ≥ q + s − 1 contains a q-mer starting at a sampled
  position, at one of the s shifts, so it is in the index. ∎
- *Consequence:* S may be **any** set of parts. CERTA picks the rarest parts
  that fit a hit budget, using index range sizes that are known before
  anything is enumerated. So the certified radius **R = |S| − 1 adapts per
  read**. This is the difference from fixed pigeonhole filters, which fail
  when any fixed part is repetitive. On HG002 that change took
  certification from 83.4 % to 94.5 % of reads (§3).

**Theorem 1′ (q-gram-lemma threshold).**
- *Statement:* with S enumerated as in Theorem 1 and any t ≥ 1, every locus
  ℓ with d(r, ℓ) ≤ |S| − t yields one cluster that contains hits of at
  least t distinct parts of S. So clusters supported by fewer than t parts
  can be skipped without loss, and the certified radius becomes
  R = |S| − t.
- *Proof:* d ≤ |S| − t edits leave ≥ t parts of S untouched, each occurring
  exactly at ℓ with a diagonal in [δ − d, δ + d], where δ is ℓ's start
  diagonal. Clusters join consecutive sorted diagonals whose gap is ≤ 4R
  (pad 2R), and these diagonals are ≤ 2d ≤ 2R apart, so all of them, and
  every candidate between them, fall in one cluster. ∎
- *Why it matters:* when every part is repetitive, the budget is spent on
  thousands of single-part hits. With t = 2 these are skipped. In the oracle
  test, t = 2 and t = 3 meet every certificate check.

**Theorem 2 (exact verification).**
- *Statement:* for a cluster of candidate diagonals [dmin, dmax], banded
  semi-global DP over diagonals [dmin − 2R, dmax + 2R] computes d(r, ℓ)
  exactly for every locus ℓ with d(r, ℓ) ≤ R whose witness part lies in the
  cluster.
- *Proof:* an alignment with ≤ R edits stays within ±R of its start
  diagonal, which is within ±R of the witness diagonal (Theorem 1). ∎

**Theorem 3 (order-independent best and second-best).**
- *Statement:* clusters are evaluated keeping the two smallest distances
  b1 ≤ b2, each cluster with limit b2 − 1. The final b1 is the minimum
  distance over all clusters. If b1 < b2, the final b2 is the second
  minimum (≤ R) or exceeds R. Neither depends on evaluation order.
- *Proof:* a cluster pruned with limit b2 − 1 has distance ≥ b2 at that
  time, and b2 never increases. ∎
- *Tested* by a metamorphic check that reverses the order.

**Corollary (the certificate).** For a certified read CERTA reports:
- **`XE` = d1:** the minimum edit distance over the whole reference, both
  strands;
- **`XB`:** whether the optimum is unique or tied (a lower bound when tied);
- **`XD`:** the exact second-best distance when it is ≤ R, otherwise "none
  within R";
- **`XR`:** the radius R itself.

**Theorem 4 (certified repeats, tier SR).** If ≥ 2 distinct exact full-read
matches are found among sampled hits, then d1 = 0 is optimal and the read
has ≥ 2 optimal loci. This holds under any scoring, so MAPQ 0 is proven.

**Theorem L (certified local alignment, tier SL; opt-in `--local t`).** The
objective here is BWA-MEM's own: the local score with soft clipping (+1 /
−4 / gap −6−1k, −5 per clipped end), which `AS` reports. The *loss* of an
alignment A of read span [x, y) is (y − x) − score(A) − its clip penalties.
- *Lemma L:* every enumerated part fully inside [x, y) that does not occur
  exactly in A adds ≥ 5 to the loss, disjointly.
  - A mismatch costs 5.
  - A deletion lies between two read bases, so it touches ≤ 1 part, for
    ≥ 7.
  - An insertion of k bases touching j ≥ 2 parts covers the j − 2 parts
    between them, so k ≥ 2 + (j − 2)·(part length), and 6 + 2k ≥ 5j.
- *Consequences:*
  - An alignment with fewer than t exact parts of S scores at most
    floor = max over spans of (y − x) − 5·(#parts inside − t + 1) − clips.
  - Its total gap length is G ≤ pad = L − floor − 7, because the loss is
    ≥ 6 + G and < L − floor.
  - So its exact parts lie on diagonals within G of each other, in one
    chain C of candidate diagonals joined at gaps ≤ pad. The alignment
    stays within ±pad of each of them.
  - It scores ≤ ub(C), the same bound with C's part mask treated as free.
- *Statement:* chains with ub > floor are evaluated by banded affine DP in
  any order, keeping the best two scores s1 ≥ s2 (initially floor) and
  skipping chains with ub ≤ s2. If s1 > floor, then s1 is the maximum local
  score over the whole reference, both strands. The reported alignment
  attains it.
- *Proof:* a better alignment either lies in an evaluated chain's band, so
  it would have been found, or in a skipped chain, so it scores ≤ s2 ≤ s1. ∎
- *Tags:* `XT:Z:SL`, `AS` (the certified maximum, clip penalties included),
  `XF` (the floor it exceeds).
- *Scope:* this certifies the score's optimality, not uniqueness. Ties
  between chains set MAPQ 0.

**Theorem P (certified proper pair, tier PR; paired input).** The objective
is the pair score σ(A1) + σ(A2), with σ the Theorem L local score. It is
taken over *proper* pairs: the mates lie on opposite strands, and the
reverse mate's first aligned base is 0 to D bases (`--max-dist`, 850) after
the forward mate's.
- *Floor:* if no enumerated part of either mate occurs exactly, Lemma L
  bounds the pair by U1 + U2. The floor is
  max(U1 + U2, L1 + L2 − 7 − 47), where the second term keeps each mate's
  gap total within the 96-diagonal host band.
- *Anchors:* every proper pair scoring above the floor has an exact part
  in some mate x, so mate x lies in one of x's chains C (an *anchor*).
  - Its partner y starts inside C's insert window.
  - The partner is found exhaustively there. The window is about 1 kb, so
    y's ~12-base parts are scanned directly in the reference, with no
    index. Lemma L then applies with up to 12 parts, so any partner not
    found scores at most the window's floor.
- *Bounds:*
  - Anchors are evaluated best bound first, with
    bound(C) = ub_x(C) + (the best bound of y near C, or U_y).
  - A partner in a chain whose pairs are already accounted for is covered,
    so the bound tightens as chains are finished.
  - Each evaluation gives `up`, a bound on every pair through C, and a pair
    actually *realized*: x at its best place in C, and y at its best in that
    place's exact window.
- *Statement:* let R1 be the best realized pair. If R1 > floor, every anchor
  of R1's pair has up = R1, and every other bound (other anchors' up,
  skipped anchors' bounds, the floor) is ≤ R1, then R1 is the maximum pair
  score over all proper pairs in the reference. ∎
- *Also reported:* S2 is the largest of those other bounds and drives a
  BWA-MEM-style pair MAPQ. Since Aim 5, S2 bounds every proper pair that
  differs from R1's pair in some mate, including a second partner inside
  the winning anchor's own window, and the mate placed elsewhere in its own
  chain (THEORY.md, Proposition P2). Before that fix, such a pair was
  bounded only by R1.
- *Tags:* `XT:Z:PR`, `AS` (each mate's score), `XF` (floor), `XS` (S2).
- *Pipeline:* when the pair cannot be certified, a mate keeps its
  single-end certificate. Both mates still go to the fallback mapper so
  that it can pair them, and its record of the certified mate is dropped.

**What is deliberately not claimed**
- **Optimal ≠ true origin:** the certificate is about the optimization
  problem, not the read's true origin. A read from a diverged repeat copy
  can align better to another copy. That is unavoidable for any aligner,
  but here it is at least *visible*: XD/XB say how close the alternative is.
- **Reported alignment:** the reported CIGAR is the affine-optimal alignment
  *at the certified locus*, with BWA-MEM's default scoring (+1/−4, gap
  −6−1k, clip −5), so downstream callers see the usual representation. The
  locus and d1 are certified; the CIGAR representation is a convention.

## 3. Evidence that the implementation meets the theorems

- **Brute-force oracle test (`tests/test_certa.cpp`):** Sellers' DP over the
  entire reference, on a genome with diverged repeat copies, a 50-copy exact
  repeat, tandem repeats and `N` runs. It checks soundness, optimality of d1,
  losslessness, uniqueness claims, SR claims (≥ 2 exact loci in the oracle),
  order independence, thread determinism, and in-memory vs memory-mapped
  identity. It runs 16 configurations on multiple seeds, including host
  budgets up to 4096 and q-gram thresholds t = 2 and 3.
- **Mutation testing:** every injected bug in the certificate logic is
  caught. The bugs were an overclaimed radius, a one-copy SR proof, a
  dropped sampling shift, a narrowed band, a wrong verification limit, a
  lost previous best, and a wrong key-suffix mask.
- **Tier SL (Theorem L):**
  - *Second oracle:* full Gotoh DP of the affine local objective over the
    whole reference, both strands.
  - *Read sets:* reads with adapter-like heads and tails, chimeras,
    adversarial reads, and every read in local-only configurations. The
    adversarial reads have a mismatch in most parts plus a planted, slightly
    worse decoy, or a long indel near a read end.
  - *Result:* about 460 SL certificates per seed are checked, and all equal
    the oracle optimum.
  - *Mutation testing* (3 seeds each): overclaims of ≥ 2–3 points are
    caught. These are a floor lowered by 3, a band narrowed by 2, a chain
    bound lowered by 2, one extra free part, and a skip test off by 2.
    Two 1–2-point overclaims survive, because they need an alignment
    exactly at the boundary: a floor lowered by 1, and a span start shifted
    by 2.
- **Tier PR (Theorem P):**
  - *Third oracle:* the best proper pair over the whole reference. It takes
    each mate's best local score per start position (local DP on the
    reversed read, both orientations), then a sliding-window maximum over
    the allowed distance.
  - *Pairs:*
    - 90 simulated pairs from unique sequence, from inside a 50-copy
      element, from diverged copies and from a tandem repeat, with 0–8
      edits and adapter tails.
    - 40 adversarial pairs. At the true locus mate 1 has 1–4 mismatches in
      distinct parts; at a planted decoy mate 1 is exact and mate 2 has 1–6
      mismatches. So the two pair scores are close, either way round.
  - *Result:* about 85 PR certificates per seed are checked, and all equal
    the oracle optimum.
  - *Fast DP kernel:* `band_score`, the 16-bit vectorized twin of
    `affine_align`, equals it on 4,000 random bands per seed.
  - *Mutation testing* (2 seeds each): the decoy pairs catch both
    partner-bound bugs (a missing nearby-partner bound, and treating
    unfinished chains as finished), and the kernel test catches a deletion
    recurrence off by one. Seven threshold mutations, off by 2–3 points in
    the floor, band pad or window threshold, survived these tests. Aim 5
    closed this gap (§7): pairs placed exactly at each threshold, a direct
    test of the window search, and an oracle for S2. Every PR mutant is now
    killed, or proven to change no claim (THEORY.md §9a).
- **Real data:**
  - 1,672,120 certified HG002 reads were compared with BWA-MEM2's unclipped
    alignments. BWA-MEM2 never has fewer edits than the certified minimum,
    except 3 reads whose windows contain a reference `N` (W4).
  - Output is byte-identical on ARM64, x86-64 and an A100. One regression
    broke this. From the local-tier commit (`3664ae6`) on, the host's wider
    band also widened the end-to-end cluster limit, so about 78 reads in 2 M
    (tandem repeats) were certified S1 on the CPU but SL on the GPU. Both
    were valid proofs, but the output was not identical. Fixed in Aim 1
    with a shared `BAND_E2E = 48`; SAM bodies and uncertified reads are
    byte-identical again on 2 M HG002 reads (CPU vs A100).
  - CI covers native x86-64 and ARM64.

## 4. End-to-end results (pilot, HG002 NovaSeq X, single-end)

**Certification (v0.6, full index, A100).** 94.46 % of reads are
certified:
- 73.6 % exact;
- 19.4 % with 1–5 edits;
- 1.5 % as certified repeats.

The remaining 5.5 % go to minibwa, which runs concurrently on the CPU.

**Optional pass 3 (CPU, `--budget3`).** It re-runs, with budgets up to 4096
hits per strand, the reads whose rarest parts did not all fit pass 2's
budget (3.6 % of reads). Measured on 2 M HG002 reads (A100 + 16 threads):

| Pass 3 | Certified | Map time |
|---|---|---|
| off (default) | 94.46 % | 0.76 s |
| budget 1024, t = 2 | 94.62 % | 1.01 s |
| budget 4096, t = 2 (Theorem 1′) | 95.34 % | 2.16 s |
| budget 4096, t = 1 | 95.95 % | 9.61 s |

The t = 2 threshold makes the large budget affordable (4.4× less map time
than t = 1) at the cost of one unit of radius. Pass 3 certifies reads but
spends more CPU on them than minibwa would, so it is off by default and
offered as a maximum-certification mode. On 20 M reads with minibwa
running concurrently (median of 3,
`results/hg002_2026-09/galaxy_timing_v07_pass3.tsv`):

| Pipeline | Certified | Wall | CPU time | vs BWA-MEM2 |
|---|---|---|---|---|
| default (pass 3 off) | 94.38 % | 23.40 s | 716 s | 4.27× |
| maximum certification (budget 4096, t = 2) | 95.28 % | 33.69 s | 859 s | 2.96× |

The default re-run matches v0.6 (22.94 s) within 2 %.

**Optional certified local alignment (tier SL, `--local t`, CPU pass).**
Measured on 2 M HG002 reads (A100 + 16 threads):

| Setting | Certified | SL | Map time |
|---|---|---|---|
| defaults | 94.46 % | — | 0.76 s |
| `--local 2` | 94.66 % | 0.20 % | 1.23 s |
| `--local 1` | 94.82 % | 0.36 % | 13.39 s |
| `--budget3 4096 --local 2` | 95.65 % | 0.32 % | 5.58 s |

Real-data check of the SL certificate (`bench/compare_sam.py`, 6,306 SL
reads with `--budget3 4096 --local 2`): BWA-MEM2 never scores above the
certified maximum, except 4 reads whose window contains reference `N`s
(BWA replaces them with random bases). BWA-MEM2 reports a lower-scoring
alignment at another locus for 575 SL reads. Of these, 570 are
low-complexity (NovaSeq poly-G artifacts), and 1 complex read has
MAPQ ≥ 20. So SL widens certification, but it does not reveal a new
BWA-MEM2 weakness.

**Why SL rescues so few reads.** The diagnostic `uncertified by parts
searched |S|` shows that after `--local 2`, 5.3 % remain uncertified:
- **3.6 % are limited by repeats.** Fewer than 6 parts fit the budget
  (0.47 % have none), so the floor is close to L.
- **1.8 % searched all 6 parts.** Their best alignment loses ≥ 25 points:
  long clips, chimeras, or many edits.

Both limits come from q = 22. A 150 bp read has at most 6 disjoint 22-mers,
and Lemma L charges 5 points per destroyed part, so no exact-part
certificate can prove an alignment that loses ≥ 30 points. Moving that line
needs:
- shorter parts (a second, smaller-q index);
- error-tolerant parts (1-mismatch neighbourhoods give ≥ 7 points per
  part, but cost ~66 lookups per part);
- or paired-end information.

**Paired-end (Theorem P; HG002 R1 + R2, the mates of the same reads).**
Certified reads, 20M pairs:

| | Single-end (R1) | Paired-end |
|---|---|---|
| reads certified | 94.38 % | **95.38 %** |
| … by the pair certificate (PR) | — | 6.99 % |
| pairs with both mates certified | — | 93.53 % |

The pair pass runs only for pairs with one certified mate and one
uncertified or tied mate (11.9 % of pairs).
- *Success rate:* on a 100k-pair sample it certifies 74 % of the tied pairs
  and 46 % of the one-uncertified pairs.
- *Pairs with both mates uncertified* are skipped: under 1 % of them can be
  certified.
- *Why the rest fail:* most fallback pairs are improper (30 %) or
  ambiguous. Among fallback pairs that minibwa pairs properly with low
  loss, 83 % have a mate with MAPQ < 20.

**Variant-calling accuracy, paired-end** (all 444,524,794 pairs mapped, chr20
vs GIAB v4.2.1, 44×; same protocol as single-end;
`results/hg002_2026-10/giab_chr20_vcfeval_pe.tsv`). Precision / recall / F1:

| Pipeline | GATK SNV | GATK indel | bcftools SNV | bcftools indel |
|---|---|---|---|---|
| BWA-MEM2 PE | 0.9900 / 0.9948 / 0.9924 | 0.9935 / 0.9932 / 0.9934 | 0.9872 / 0.9943 / 0.9907 | 0.9311 / 0.9203 / 0.9257 |
| **CERTA PE + minibwa** | **0.9937** / 0.9925 / **0.9931** | **0.9954** / 0.9926 / **0.9940** | **0.9910** / 0.9926 / **0.9918** | **0.9346 / 0.9238 / 0.9291** |

- *F1:* higher than BWA-MEM2's for SNVs and indels with both callers
  (+0.0007, +0.0006, +0.0011, +0.0034).
- *Where the gain comes from:* the opposite of single-end. Precision is
  higher (GATK SNV false positives 452 vs 719), and recall is slightly
  lower (false negatives 533 vs 369).
- *Certification on the whole sample:* 95.38 % of the 889 M reads.

**Three samples (paired-end, chr20 vs GIAB v4.2.1).** These are HG002
(NovaSeq X) and two more GIAB samples on another instrument: HG001 and HG005,
NovaSeq 6000 2×151 PCR-free, study PRJNA734598 (SRR14724533, SRR14724528).
Every pair was mapped (722 M and 726 M reads), with the same protocol as
above. For HG001 and HG005 GIAB provides only the standard benchmark BED,
not HG002's "noinconsistent" variant.
(`results/hg002_2026-10/giab_chr20_vcfeval_hg001_hg005_pe.tsv`)

F1, CERTA / BWA-MEM2 (Δ):

| Sample | Certified | GATK SNV | GATK indel | bcftools SNV | bcftools indel |
|---|---|---|---|---|---|
| HG002 | 95.38 % | 0.9931 / 0.9924 (+0.0007) | 0.9940 / 0.9934 (+0.0006) | 0.9918 / 0.9907 (+0.0011) | 0.9291 / 0.9257 (+0.0034) |
| HG001 | 92.93 % | 0.9925 / 0.9921 (+0.0004) | **0.9909 / 0.9922 (−0.0013)** | 0.9921 / 0.9915 (+0.0006) | 0.9469 / 0.9456 (+0.0013) |
| HG005 | 92.86 % | 0.9926 / 0.9917 (+0.0009) | 0.9938 / 0.9925 (+0.0013) | 0.9926 / 0.9916 (+0.0010) | 0.9472 / 0.9454 (+0.0018) |

- *Overall:* CERTA has the higher F1 in 11 of 12 comparisons.
- *The exception:* HG001 GATK indels, −0.0013. That is 88 false
  positives against 75, and 94 false negatives against 80, out of about
  9,900 true indels.
- *A consistent trade-off:* higher precision for SNVs in all six
  sample–caller pairs (GATK SNV false positives: HG001 508 vs 743, HG005
  537 vs 830), with slightly lower recall.
- *Certification:* lower on the NovaSeq 6000 samples, 92.9 % against
  95.4 %. These reads are 151 bp and carry more sequencing errors, so fewer
  of them fall within the certified radius.

The paired pipeline is **4.0× faster than BWA-MEM2 in paired mode**, with
5.7× less CPU time. It is not yet faster than minibwa (20M pairs, galaxy,
A100 + 64
threads, one run; `results/hg002_2026-10/galaxy_timing_pe.tsv`):

| Pipeline | Wall | CPU time |
|---|---|---|
| BWA-MEM2 | 311.22 s | 16,913 s |
| minibwa | 61.65 s | 3,598 s |
| CERTA (A100, 32 threads) + minibwa fallback (32 threads) | 76.96 s | 2,984 s |

- *Pair pass cost:* about 1,050 CPU-s. It went from 54 s to 0.35 s per
  100k pairs during development.
- *Fallback size:* a pair with one uncertified mate sends **both** mates to
  minibwa, so that minibwa can pair them. That is 12.9 % of reads, and the
  costliest ones.
- *Options:* send only the uncertified mate, with no mate rescue in the
  fallback; restrict the pair pass to the one-uncertified pairs; or move
  the pair pass to the GPU.
- *Update 2026-10-15 (`results/hg002_2026-10/galaxy_timing_pe_v2.tsv`):*
  a profile of the pair pass showed its time in the window search (an
  O(n³) span bound, a binary search per reference position) and in two
  DPs. Four changes, with byte-identical SAM and fallback output:
  - Lemma S's closed form for the span bound, used only when its
    hypotheses hold (the lemma is machine-checked);
  - a hash for the window's part keys;
  - band_score's deletion pass vectorized (16-bit lanes);
  - affine_align's per-cell bounds checks and branches hoisted.

  CERTA's own side at 32 threads drops from 56.7 s to 47.7 s, and the
  pipeline's CPU from 2,984 s to 2,680 s. End to end, the best split (CERTA
  24 threads, minibwa 40) takes 72.8 s, against 77.0 s before and 61.6 s
  for minibwa.
- *Update 2026-10-16 (`results/hg002_2026-10/galaxy_timing_pe_v3.tsv`,
  `pe_profile.txt`): level with minibwa.* A profile showed two things.
  - *The GPU and the CPU took turns.* Each batch ran its GPU passes
    (13.9 s in all), then the CPU pair pass and output (29 s at 32 threads),
    then the next batch. Now the GPU stage of batch n + 1 runs while the CPU
    stage writes batch n. Three batches are in flight; one CPU stage at a
    time, in input order, so the output does not change.
  - *The pair pass was 40 % traceback DP.* `affine_align`'s fill is now
    vectorized (4 lanes; THEORY.md "The traceback DP, vectorized"), with
    every traceback byte equal to the scalar fill's. A first version was
    slower than scalar: under Ubuntu's `_FORTIFY_SOURCE`, its `memcpy`s at a
    variable offset into stack arrays became checked library calls. With
    plain vector accesses it is 1.65× faster per cell.

  CERTA's own side at 24 threads: 56.9 s → 39.0 s (CPU 909 s → 801 s). End
  to end, CERTA 20 threads + minibwa 44: **60.6 s, against 61.4 s for
  minibwa** in the same session, with 28 % less CPU (2,603 s vs 3,595 s).
  SAM and fallback files are byte-identical to the previous build.
- *What is left is the fallback.* On 2 M pairs, 6.4 % of pairs go to minibwa
  (an earlier note here said 12.9 % of reads; it is 6.4 % of reads):
  - 2.8 % of pairs have both mates uncertified; minibwa spends 1,118 µs per
    read on them;
  - 3.7 % have one, at 461 µs per read;
  - the average over all reads is 88 µs.

  The fallback is 6.4 % of reads but 54 % of minibwa's own CPU for the whole
  sample. A smaller minibwa batch (`-K`) would overlap better, but it changes
  minibwa's per-batch insert-size estimate: with `-K 20m`, 18,658 fallback
  records change. So the fallback keeps minibwa's defaults until an
  accuracy run says otherwise.

**Variant-calling accuracy on chr20 vs GIAB v4.2.1, same caller for both
pipelines.** Both are single-end, 22×. Precision / recall / F1:

| Pipeline | GATK SNV | GATK indel | bcftools SNV | bcftools indel |
|---|---|---|---|---|
| BWA-MEM2 | 0.9898 / 0.9778 / **0.9837** | 0.9851 / 0.9603 / **0.9726** | 0.9888 / 0.9778 / 0.9833 | 0.9341 / 0.9125 / 0.9231 |
| CERTA v0.6 + minibwa (run 3) | 0.9897 / 0.9767 / 0.9832 | 0.9853 / 0.9619 / **0.9734** | 0.9884 / 0.9778 / 0.9831 | 0.9339 / 0.9144 / **0.9240** |
| **CERTA v0.6 + minibwa, MAPQ as bwa-mem (run 4, final)** | 0.9886 / **0.9810** / **0.9848** | 0.9851 / **0.9623** / **0.9736** | 0.9873 / **0.9809** / **0.9841** | 0.9335 / **0.9148** / **0.9241** |

In the final version, F1 is **higher than BWA-MEM2 for SNVs and indels
with both callers**:

| Measure | Δ F1 |
|---|---|
| GATK SNV | +0.0011 |
| GATK indel | +0.0010 |
| bcftools SNV | +0.0008 |
| bcftools indel | +0.0010 |

The gain comes from recall (GATK SNV: 228 more true positives); precision
is slightly lower (0.9886 vs 0.9898). Raw rtg vcfeval output:
`results/hg002_2026-09/giab_chr20_vcfeval.tsv`.

How the accuracy was reached (each step measured on the same protocol):
1. **First version: indels suffered.** Reporting edit-distance alignments
   gave bcftools indel F1 0.790. Reporting the affine-optimal alignment at
   the certified locus fixed it (0.924).
2. **Heuristic tier removed.** The uncertified "S2" tier was enriched at
   false positives, and on synthetic truth 50 of its 62 confident placements
   were wrong. It is off by default.
3. **MAPQ realigned with BWA-MEM2.** A remaining SNV recall gap came from
   MAPQ: reads with the same alignment as BWA-MEM2 got MAPQ 10–15 where
   BWA-MEM2 gives 22–36, because BWA-MEM's default MAPQ formula differs from
   the one first used. Run 4 uses BWA-MEM's formula.

**Speed.** 20 M reads on galaxy: an A100 plus 64 ARM cores. Every tool
gets the same 64 CPU threads; for CERTA that is 16 threads plus the GPU
for the certified tier and 48 threads for concurrent minibwa. Median of 3
runs in one session, spread < 2 %
(`results/hg002_2026-09/galaxy_timing_v06.tsv`):

| Pipeline | Wall | CPU time | vs BWA-MEM2 | vs minibwa |
|---|---|---|---|---|
| BWA-MEM2 | 99.81 s | 5,230 s | 1.00× | 0.27× |
| minibwa | 26.58 s | 1,367 s | 3.75× | 1.00× |
| **CERTA v0.6 (A100) → minibwa, concurrent** | **22.94 s** | **707 s** | **4.35×** | **1.16×** |
| *CERTA certified tier alone (94.4 % of reads)* | *11.02 s* | *69 s* | — | — |

The complete pipeline is **4.35× faster than BWA-MEM2 with better F1**. It
also uses **7.4× less CPU time**, which matters for cost per genome. It is
1.16× faster than the fastest tool (minibwa), which has more provable errors
(§1).

The pilot did not reach 10× over BWA-MEM2:
- the certified tier alone processes 94 % of reads in 11 s (9× BWA-MEM2's
  wall time);
- the 5.6 % of reads that are not certified take longer than that on 48
  CPU threads (Amdahl's law; `bench/RESULTS.md` §5b).

## 5. Contributions, stated for a paper

1. **An adaptive pigeonhole certificate for read mapping.** Choosing the
   rarest disjoint parts under a hit budget keeps the lossless guarantee and
   lets the certified radius adapt per read. With a certified-repeat rule
   and exact best/second-best verification, it certifies ~94 % of NovaSeq X
   reads, including the exact optimum, uniqueness or ties, and the second
   best within the radius.
2. **A method to measure aligner errors on real data without simulation.**
   Certificates turn "disagreement with a reference aligner" into *provable*
   errors. On HG002 every tested heuristic aligner has them, at rates that
   grow with speed.
3. **A two-pass GPU/CPU design.** The easy majority of reads is certified on
   the GPU with uniform warps; a compacted second pass handles the hard ones;
   the fallback runs concurrently on the CPU. The design is deterministic,
   with identical output on CPU, GPU, x86-64 and ARM64.
4. **Speed and accuracy together (pilot).** On GIAB HG002 chr20, F1 is
   higher than BWA-MEM2's for SNVs and indels with both GATK and bcftools,
   at 4.35× BWA-MEM2's speed and 7.4× less CPU time.

## 6. Limitations a reviewer will raise (and what would answer them)

- **One sample, single-end, chr20.** Answer: the Part-4 protocol, meaning
  several GIAB samples, paired-end, whole genome and DeepVariant as a second
  caller.
- **5.5 % of reads carry no certificate.** They get whatever minibwa gives.
  The certified local-alignment mode (Theorem L) rescues only 0.2–0.4 %:
  - 3.6 % of reads are repeat-limited;
  - 1.8 % lose ≥ 25 points, beyond what 6 exact 22-mer parts can certify.

  Answer: paired-end rescue, shorter or error-tolerant parts (§4).
- **The certificate uses edit distance, while calling uses affine
  scores.** The reported alignment is affine-optimal at the certified locus,
  but the locus itself is edit-optimal. For almost all reads these coincide,
  but not by proof.
- **MAPQ is not calibrated from first principles.** It reuses BWA-MEM's
  formula. Calibrating it on certified second-best distances is future work.
- **Speed vs the fastest tool depends on the configuration** (2026-10-18,
  64 cores, paper Table 2). Interval mode, which needs no other aligner, is
  1.60–2.00× minibwa on pairs on CPU alone and 2.11× with an A100; with
  minibwa as the fallback, pairs gain only 1.20–1.36×, because minibwa's
  work on the ~6 % of pairs with an uncertified mate exceeds CERTA's on
  all of them. Answer: fewer uncertified reads, or interval mode where its
  proven intervals suffice.
- **The GPU is part of the claim.** The comparison is A100 + 64 threads vs
  64 threads. Cost per genome (GPU hour vs CPU hours) should be reported
  alongside wall time.
- **Precision is slightly lower than BWA-MEM2's in single-end mode**
  (GATK SNV 0.9886 vs 0.9898), with F1 higher through recall. In paired
  mode the opposite holds on all three samples: higher precision, slightly
  lower recall (§4).

## 7. Open items (state as of 2026-10-03)

*Update 2026-10-08, Aim 1 of
[EXACT_ALIGNMENT_PROPOSAL.md](../EXACT_ALIGNMENT_PROPOSAL.md) §6.1.*
A second index (q = 15, 10 parts per read) breaks the q = 22 ceiling.
Single-end, the certified share rises:
- HG002: from 95.65 % to 96.29 % at about 0.27 ms per residual read, or to
  96.83 % with large budgets;
- HG005: from 93.53 % to 94.66 %, or to 95.31 %.

The new certificates show no unexplained violation against BWA-MEM2.

*Update 2026-10-11, Aim 2 (proposal §6.2).* Tier R (`certa repeat`,
`include/certa/repeat.h`) handles repeats: long rare or absent parts found
by anchor extension, copy-sampled lower bounds, and count-based tie proofs.
- **Repeat-limited reads** (|S| < P, 3.55 % of HG002, 3.99 % of HG005):
  50.5 % and 43.8 % become exact or proven ambiguous; 53.2 % and 46.0 %
  with larger limits.
- **Near-perfect copies among them:** 79–86 %.
- **The rest** are lossy reads beyond any exact-part certificate.
- **Validity:** 0 violations of tier R's bounds against BWA-MEM2 on 180,804
  residual reads.
- **Tests:** all 6 injected faults are caught.

*Update 2026-10-11, Aim 3 (proposal §6.3).* `certa map --intervals`
(`include/certa/interval.h`) writes every read itself; nothing goes to a
fallback mapper. Each record carries:
- [XL, XU], which contains OPT;
- XQ, a certified lower bound of BWA-MEM's MAPQ formula.

Findings:
- **HG005 (725.8 M reads):** 84.9 % of records are proven exact (XL = XU).
  The pre-registered go criterion is met. Variant-calling F1 against
  CERTA + minibwa, with no fallback:
  - GATK SNVs 0.9933 (+0.0007), GATK indels 0.9938 (0);
  - bcftools SNVs 0.9937 (+0.0011), bcftools indels 0.9464 (−0.0008);
  - above BWA-MEM2 in all 4.

  Cost: +13 % CPU over CERTA + minibwa on 2 M HG002 pairs.
- **XQ as MAPQ fails (F1 0.31–0.78): it is biased toward the reference.**
  Beyond the radius, Lemma E bounds SUB at one point per edit (a locus
  that matches all but d clipped bases); a SNV costs 5.
  - A one-SNV read at radius 2 gets XQ ≈ 10; reads matching the
    reference get 26–40.
  - Callers' MAPQ ≥ 20 filters therefore drop exactly the variant
    reads.
  - A SUB bound from Lemma L (exact parts, as tier SL builds its floor)
    would fix this (open).
- **Edit-distance certificates are not affine certificates.** On HG002,
  27 % of S1 reads have XL < XU: another locus, where a clipped or gapped
  alignment needs more edits, could still score higher under BWA-MEM's
  scores. Lemma E gives the bound; running tier SL on these reads would
  close most of them (open).
- **Over-confident MAPQ of certified reads at small radius** (a
  pre-existing issue, shared with the current pipeline). With R = 0 or 1
  and no second locus within R, `mapq_of` assumes sub = 19 and reports
  MAPQ 60, although a locus one edit further was never searched. On
  simulated pairs:
  - these reads cause most misplaced records with MAPQ ≥ 20: 257 of the
    539 (208 at R = 0, 49 at R = 1);
  - their XQ is 7–18, which flags them correctly.
- **XQ bounds the formula, not the error rate.** On simulated pairs:
  - all 173 misplaced records with XQ ≥ 20 are PR pairs whose score-optimal
    placement is not the true origin;
  - 150 of them are tandem shifts, which by definition count as the same
    locus;
  - BWA-MEM2 places the other 23 at the same loci, with MAPQ 0–27.
- **Bug fixed:** `Workspace::ub_magic` was never initialised. Heap memory
  that last held a Workspace kept the magic while its touched list was
  overwritten, so the PR pass wrote out of bounds. Interval mode's
  allocation pattern exposed it on HG005; it had been latent since tier PR.

*Update 2026-10-12, Aim 4 (proposal §6.4).* `certa check` (`src/check.cpp`)
is an independent checker of the interval output (Theorem V). It shares no
alignment, bound or lookup code with the mapper. It recomputes Lemma L,
enumerates every probe occurrence from an index that `certa check-index`
has validated against the FASTA, rebuilds the chains, and searches with its
own banded DP.
- **Honest records:** all accepted, 4 M each on HG002 and HG005 (2 M pairs
  each).
- **Faults false by construction:** all 23,733 injected into real SAM were
  rejected, with 0 false alarms.
- **Mapper mutants:** all 6 flagged on real data. One is a wrong constant in
  the mapper's own Lemma L.
- **Cost:** 89–134 µs CPU per record, about the mapper's. Index validation
  takes 34 s once per reference.
- **A real bug found:** records whose alignment runs off the circular chrM
  kept XL and XQ after the alignment was dropped. Fixed.
- **Gap:** the XQ of PR records (about 4 %) rests on the pair certificate,
  which the checker did not re-derive. Closed on 2026-10-14 (below).

*Update 2026-10-13, Aim 5 (proposal §6.5).* [THEORY.md](THEORY.md) holds
written proofs of every bound the implementation uses. Each one is traced to
its code, its oracle tests and the mutants that break it.
- **Exhaustive checks** of the finite lemmas (`tests/verify_lemmas.cpp`, in
  CI): Lemma L on 52.8 M and Lemma H on 502 M transcript-layout pairs, no
  violation. Their part-length conditions are necessary: 2- and 3-base parts
  give counterexamples.
- **Mutation tests** (`tests/mutants.tsv`, `tests/run_mutants.sh`): 56
  faults, 50 killed. The other 6 change no claim on any input, which
  THEORY.md §9a proves.
- **Two bugs in tier PR's S2, found and fixed:**
  - S2 did not bound a second partner inside the winning anchor's window, or
    the mate elsewhere in its own chain. This closes the S2 open item
    listed here before.
  - The window search's band pad went negative for high thresholds, missing
    an exact partner.
  - On 2 M HG002 pairs, 0.29 % of records change MAPQ (4,419 fall below
    20).
  - HG005 re-run with the fixes: variant-calling F1 stays within the
    Aim 3 margin and above BWA-MEM2 in all 4 comparisons. GATK SNVs
    0.9933, GATK indels 0.9935, bcftools SNVs 0.9938, bcftools indels
    0.9467.
- **Checker:** it now also rejects a record whose own alignment scores above
  its XU.
- **Open:** review of THEORY.md by a second person. The optional Lean 4
  mechanization has started: Lemmas L, E, H, S and W, Lemma C's core,
  Theorem G, the long-read band pruning and the vector DP row are
  machine-checked (`lean/`, 2026-10-15).

*Update 2026-10-13, Aim 6 (proposal §6.6, criterion pre-registered).*
`tools/lr_ceiling` measures whether long reads fit the exact-part
certificate. For each read it takes the Lemma L ceiling U₁ over greedy
parts that occur ≤ 16 times, and compares it with the score σ of the read's
minimap2 alignment.
- **Sample:** 242,534 real HG002 HiFi reads against GRCh38.
- **Result:** 93.29 % fit (σ > U₁), against a go threshold of 90 %. A
  typical read's loss (247 points) is a tenth of its ceiling (2,605).
- **Checks:** every part count and U₁ was recomputed independently, and
  every CIGAR checked against the reference. All agree.
- **Misses:**
  - 2.3 % clipped (mostly split alignments);
  - 2.1 % high loss (mostly an SV gap of ≥ 500 bases);
  - 1.8 % repeat-limited;
  - 0.5 % unmapped.
- **Not measured:** the verification cost, i.e. the exact alignment of the
  candidate loci.
  - *Measured on 2026-10-14 (proposal §6.6):* a known alignment of score
    lb narrows the band to L − lb − 7 (about ±250). A banded DP over the
    chains that need it (median 1) proves OPT for 226,330 reads (93.32 %),
    every read that fits.
  - *Cost:* median 3.1 ms of DP per read, 0.86× minimap2's CPU in all.
    Three changes on 2026-10-15 took it there from 20 ms and 2.8×:
    a vectorized DP row, batched index lookups, and pruning the DP to the
    cells that can beat lb. None changes a certificate.
  - *OPT vs minimap2:* OPT beats minimap2's alignment under CERTA's scores
    for 17.8 % of them.
  - *Check:* on 2,005 reads, a full DP around the known locus never beat
    the proven OPT.

*Update 2026-10-14 (proposal §6.4).* The Aim 4 gap is closed: `certa check`
verifies the XQ of PR mates too.
- **What it checks:** BWA-MEM's pairing rule, with XS proven from both
  records. Each mate's pair-certificate parts come in the new witness tag
  `XY`. Pairs with no exact part are bounded by Lemma L; the rest by
  exhaustive DP, pairing every chain of a mate's part occurrences with every
  partner placement within D.
- **Real data:** 8 M records (2 M pairs each of HG002 and HG005) accepted.
  All 303,682 XS claims were proven, and 24,306 of 24,306 false-by-
  construction faults rejected, among them three new pair faults.
- **Cost:** 181–232 µs CPU per record, up from 89–134 µs.
- **A coarse first version:** it left 12 honest HG002 pairs unproven. The
  mapper's S2 held in every case; the check was the problem.
- **Mutation tests:** 10 new mutants of the pair check (C10–C19), all
  killed. The suite now has 66 mutants: 60 killed, 6 proven equivalent.

- **Speed in paired mode:** 5.1× faster than BWA-MEM2 and level with minibwa
  since 2026-10-16 (60.6 s vs 61.4 s, 20M pairs; it was 77 s). What is left
  is the fallback's hard reads, 6.4 % of reads at about 54 % of minibwa's
  whole-sample CPU. Sending only the uncertified mate
  (`--fallback-mates single`) barely helps. Candidates: certify more of
  those reads; a smaller minibwa batch, once an accuracy run shows its
  changed pairing statistics do no harm.
  - *Update 2026-10-18.* The smaller batch is safe: on HG005, minibwa
    `-K 10m` changes 16,761 chr20 records of CERTA + minibwa but no F1 in
    the fourth decimal (GATK and bcftools, SNVs and indels;
    `bench/results/hg002_2026-10/hg005_fallback_batch_f1.tsv`).
  - *Where paired time went (a2, 64 cores, 20M pairs, build 80646c3).*
    CERTA alone took 35 s, 34 s of it in the input stage: a serial
    boundary scan, a serial interleave of the mates and the freeing of
    temporary records. minibwa on the fallback pairs alone took 1,939 CPU-s,
    more than CERTA's 1,553 for every pair, so the fallback pipeline is
    bound by minibwa's work and can gain at most ~1.3× over minibwa.
    Interval mode needs no fallback.
  - *Three changes, output byte-identical* (both hosts, CPU and A100, every
    mode): FASTQ boundaries found on several threads and records parsed in
    place; every pass of a batch run chunk by chunk (64 reads), so tier R's
    slow repeats no longer idle the other cores at each pass's end; tier R
    and `find_copies` prefetch the reference at the hit 16 ahead. a2: paired
    interval mode 51.7 → 41.0 s (minibwa 79.5 s: 1.55× → 1.94×), CERTA
    alone with fallback 34.7 → 26.2 s, fallback + minibwa (`-K 10m`)
    65.2 → 59.4 s, single-end interval mode 22.8 → 19.3 s.
  - *Variance on a2:* with the kernel's NUMA balancing, two of four
    fallback + minibwa runs took 90–92 s instead of 58–60 s, with ~40 %
    more page migrations, more kernel time and fewer CPU-seconds (cores
    idle). Under `numactl --interleave=all` four runs took 57.6–58.0 s
    (`bench/results/hg002_2026-10/a2_fallback_numa.txt`), so the final
    protocol runs every a2 pipeline that way and takes medians of three
    runs, paired too.
  - *Final protocol, build 84dfa83* (`final_{galaxy,a2}.tsv`; the
    80646c3 runs are kept as `final_*_80646c3.tsv`): paired interval mode
    39.45 s on galaxy and 39.49 s on a2 (1.60× and 2.00× minibwa; it was
    1.38× and 1.55×), 29.85 s with the A100 (2.11×, was 1.77×); the
    fallback with `-K 10m` 1.20× and 1.36×. Single-end interval mode 1.76×
    and 1.92× (was 1.63× and 1.57×), 2.10× with the A100.
  - *Tried after that, no gain (not kept):* reading the DP kernels'
    previous row as stored vectors and shifting it in registers, instead of
    the `+ b + 1` loads that straddle two stores, gave identical output but
    was 1.5–2 % slower on both hosts (paired interval mode 40.1 vs 39.5 s
    on a2, 40.2 vs 39.5 s on galaxy). GPU batch size: 0.5 M and 1 M reads
    are equal (29.7 s), 2 M and 4 M slower (30.5, 32.6 s).
  - *Where paired interval mode's CPU goes now* (galaxy, perf): tier R 26 %
    (its index lookups 9.5 %, of which `find_copies` 3.2 %), the DP fills
    for CIGARs and the pair certificate ~22 % (`certify_pair`'s tracebacks
    6.3 %, mate rescue's 10 %), passes 1–2 22 %, the pair certificate's
    window scans 7 %. Each remaining fix is worth 2–5 %; larger gains need a
    GPU pair pass or GPU tier R.
  - *Two of them, kept* (77b1df2, output identical): `find_copies` reuses
    the read's own index ranges, and `band_score_vec` keeps its window as
    16-bit lanes and skips the start mask in rows that exclude no start.
    Final protocol (`final_{galaxy,a2}.tsv`; 84dfa83's kept as
    `final_*_84dfa83.tsv`): paired interval mode 37.82 s on galaxy and
    37.97 s on a2 (1.66× and 2.09× minibwa), 28.31 s with the A100 (2.22×);
    single-end interval mode 1.88× and 2.02×, 2.32× with the A100. The same
    changes in `affine_fill_vec16` and `band_distance_vec16` were neutral
    (0.5 % faster on a2, 1 % slower single-end on galaxy) and not kept.
- **Maximum certification and the slow, reliable mappers** (2026-10-22,
  ea812d8; `bench/results/hg002_2026-10/slow_*`).
  - *Options:* `certa map --full-cert` turns on pass 3 (budget 4096) and
    pass L (t = 2); `--index2 q15.cidx` runs every pass again on a second
    index for the reads still uncertified (the Aim 1 cascade). The default
    output is unchanged (byte-identical to 77b1df2 in every mode).
  - *Certified, 2 M HG002 reads:* single-end 94.46 → 95.65 % with
    `--full-cert`, 96.83 % with `--index2` too; paired 95.42 → 96.33 →
    96.87 % of reads (both mates: 93.59 → 94.90 → 95.43 % of pairs).
  - *Paired `--full-cert` was 6.4× slower than the default:* the pair
    certificate's partner-chain lookup scanned every chain of the other
    mate, quadratic once 4096 hits per part are enumerated (80 % of the
    time). Binary search over the sorted chains (THEORY.md, Theorem P):
    22.1 → 5.9 s, output identical; mutant P17 killed.
  - *Speed* (galaxy, 64 cores, interval mode, CPU, medians of three; Yara
    and RazerS 3 one run): single-end 2.67 / 2.97 / 6.21 s (default /
    `--full-cert` / `--index2`), paired 4.95 / 6.85 / 12.28 s. Yara 1.0.5
    `-e 4 -y full` (every alignment with ≤ 4 % errors): 138.2 s and
    211.2 s, 22–46× and 17–31× CERTA's wall time and 9–17× and 7–12× its
    CPU (Yara keeps 13–19 of 64 cores busy). RazerS 3 `-i 96 -rr 100`
    (100 % recognition at ≥ 96 % identity) on the first 200,000 reads:
    239.6 s against CERTA's 3.2 s (75× wall, 510× CPU); its paired mode
    aborts with `std::bad_alloc` at the end of every run (also with
    `-le 150`). The A100 is slower than the CPU at 2 M reads (fixed
    start-up costs). CERTA with `--index2` holds two 23 GB indexes
    (~46 GB RSS; Yara 5.5 GB).
  - *Reliability* (`bench/weaknesses.py`, split by certificate): on the
    1.91 M reads certified end-to-end, Yara single-end places none
    elsewhere with more edits than the proven minimum; 89 differ only
    under BWA's affine scores (equal edits); 14 (0.001 %) are unmapped
    although ≤ 5 edits exist. Paired, Yara places 2,442 mates (0.07 %)
    elsewhere with more edits than proven, the price of its pairing.
    RazerS 3 maps every end-to-end certified read, but flags its first
    match primary, not its best: by primary it is "wrong" for 6.6 %, by its
    best record for 0.01 % (affine-score differences). Both leave ~56 % of
    local (SL) certificates unmapped: end-to-end mappers cannot clip.
  - *Long reads:* Winnowmap 2.03 has no short-read mode; on 242,534 HiFi
    reads it takes 610.4 s (88.5 s of it indexing the FASTA: its dumped
    index crashes on load) against 163.4 s for minimap2 + `lr_ceiling
    --certify` (3.7× wall, 4.1× CPU), which proves the optimum for
    93.32 %. On Winnowmap's alignments `lr_ceiling` proves the optimum for
    93.31 %, above Winnowmap's own alignment for 70,653 reads (minimap2:
    40,265), under CERTA's +1/−4 scores; Winnowmap's preset scores +2/−4,
    so this is not an error count.
  - *What is left for their patterns* (`bench/headroom.py`,
    `slow_headroom.txt`, 2 M reads single-end): Yara maps 96.83 % within 6
    edits, as many as CERTA certifies with `--full-cert --index2`
    (96.83 %), but not the same reads (CERTA's 1.09 % SL is out of Yara's
    reach). Of the 3.17 % CERTA leaves uncertified, Yara at 12 edits
    (`-e 8`, 2,405 s for these 63,368 reads) finds the best end-to-end
    alignment within 6 edits for 0.62 % of all reads (mostly with parts
    over the hit budget: 0.53 %), within 7–11 edits for 0.60 %, and none
    within 12 edits for 1.82 %. Default settings leave 1.98 %
    repeat-limited reads within 5 edits; `--full-cert --index2` already
    took three quarters of them.
  - *Patterns worth taking, by that headroom:* parts chosen by their
    rarity (Winnowmap's weighted minimizers, the optimal seed solver; the
    parts are fixed today and only their first 22-mer is looked up; Lemma
    L allows any disjoint parts, as in `lr_ceiling`) and repeats handled as
    index intervals or merged adjacent parts (Yara's bidirectional FM
    index) for the 0.62 %; parts that may match with one error (Yara's
    search schemes; radius 2|S| − 1, up to 11) for the 0.60 %; bit-parallel
    banded verification (Myers/Hyyrö, RazerS 3 and Yara) for the cost of
    `--full-cert`. Not worth taking: exhaustive search at a larger radius
    (Yara at 12 edits costs 1.3 CPU-s per read), RazerS 3's scan of the
    genome, and end-to-end-only alignment: the remaining 1.82 % have no
    end-to-end alignment within 12 edits, so only clipped (local)
    certificates or larger radii can reach them.
  - *Their payoff, measured* (`tools/payoff`, `bench/payoff_union.py`,
    `slow_payoff.txt`): a standalone certifier of the minimum end-to-end
    edit distance, one part choice per technique, on the reads CERTA leaves
    uncertified. Every certified distance agrees with Yara's (31
    comparisons, none with fewer edits; no two techniques disagree on a
    read). On top of `--full-cert` (95.65 %), on the q = 22 index alone:
    rarest disjoint 22-mers +0.85 points, one-edit parts (the six rarest
    22-mers' neighbourhoods; radius up to 11) +1.55, both in a cascade
    +2.00 (97.65 %) for 344 more CPU-s in this prototype (CERTA's whole
    run: 105 s; 4.0 ms per residual read). That beats
    `--index2` (96.83 %, two 23 GB indexes); with it too, 97.75 %; every
    technique on both indexes, 97.99 % for 20× the CPU, against a ceiling
    of 98.18 % (the 1.82 % without an end-to-end alignment within 12 edits).
    Long exact parts add little once these run, and cost most here (24–41
    ms per read through hash-index extension). CERTA's own parts with
    support 1 at budget 4096 and every hit verified on its own band
    (pass 3 uses support 2 and drops wide clusters) give +0.57 points
    alone. Next: rare parts and one-edit parts in the mapper itself.
- **Pass V: adaptive parts in the mapper** (2026-10-22; THEORY.md Lemma N,
  Theorem A; `bench/results/hg002_2026-10/pass_v_check.txt`). For reads
  still uncertified after pass 3: the rarest disjoint q-mers as parts
  (radius |S| − 1), then q-mers searched with one edit (radius 2|S| − 1, at
  most 12). On with `--full-cert` (unless `--witness`), or
  `--adaptive-parts 1|2`; opt-in, so the default output is byte-identical
  to 77b1df2 (14 files, CPU and A100), and pass V off reproduces 77b1df2's
  `--budget3 4096 --local 2`; with it on, CPU = A100 and fused = serial.
  - *2 M HG002, interval mode, 64 cores:* single-end `--full-cert` 95.65 →
    97.63 % certified, 3.0 → 5.6 s, 22.7 GiB (the old maximum, `--index2`:
    96.83 %, 6.3 s, 44.2 GiB; both: 97.73 %, 8.4 s). Paired 96.33 → 97.45 % of
    reads (both mates 94.90 → 96.18 % of pairs), 6.9 → 12.5 s. The
    prototype predicted 97.65 %.
  - *Tests:* the oracle at radius 12 (873 reads, 6 settings, 3 seeds, ARM
    and x86), Lemma N directly, the DP against brute force, the hit budget;
    11 mutants V1–V11 killed (141 in all: 133 killed, 8 equivalent).
  - *Against the slow mappers* (`slow_reliability.txt`): within their own
    error rate (4 % of the read length) Yara and RazerS 3 miss none of the
    reads CERTA certifies; every read they leave unmapped has its proven
    optimum beyond that rate (17,836 single-end reads, 0.89 %; distances
    rechecked against the reference). CERTA now certifies those too.
  - *Where pass V's time goes* (galaxy, perf, paired `--full-cert`): about
    56 % of the run, of which the exact verification of clusters is 40
    points (the banded DP itself 44 % of a single-end run), its lookups and
    part choice 17. A read in a repeat has up to 4,096 hits per strand, and a
    unique read keeps b2 at R + 1, so every cluster that could hold a second
    locus within R must be verified: that is the claim d2 makes.
  - *Faster, output unchanged* (byte-identical to c9ee1c7 in 15 outputs:
    default and `--full-cert`, single-end and paired, fallback and interval
    mode, CPU and A100; ARM = x86): batched, prefetched index lookups in
    pass V; the reference windows of clusters prefetched six ahead; and a
    support bound (THEORY.md, Theorem A) that skips a cluster whose parts
    imply more edits than the current limit. 2 M reads, `--full-cert`:
    single-end 5.58 → 4.96 s (CPU 260 → 221 s), paired 12.45 → 11.10 s
    (675 → 591 s). Two alternatives were measured and not kept: a two-phase
    verification at a smaller radius (the same speed, but a changed radius
    for some records), and a smaller pass V budget (1,024 hits: 24 % faster
    single-end, but 0.34 points fewer certified).
  - *At full scale* (`full_cert_{galaxy,a2}.tsv`, `full_cert_certified.txt`;
    20,001,073 reads or pairs, interval mode, medians of three, final
    protocol's pinning): certified single-end / paired 94.38 / 95.39 %
    (default), 95.59 / 96.29 % (`--full-cert --adaptive-parts 0`), 97.61 /
    97.44 % (`--full-cert`), the same on both hosts, CPU and A100. Wall,
    galaxy CPU: 15.8 → 20.9 → 38.7 s single-end, 37.6 → 59.8 → 100.1 s
    paired (A100: 13.0 / 17.5 / 34.7 and 28.4 / 49.2 / 88.6 s); a2: 18.2 →
    19.5 → 40.9 s and 37.8 → 58.1 → 107.1 s. Peak RSS 29.1 GiB single-end,
    36.2 GiB paired (39.5 GiB with the A100; 21.3 and 21.6 GiB since the
    memory changes below). So `--full-cert` costs about 2.5×
    the default: slower than minibwa (on galaxy 28.8 s single-end and
    63.0 s paired; also with the A100 and on a2), 2.6× and 3.3× faster than
    BWA-MEM2, and it certifies 3.2 and 2.1 points more than the default
    (pass V alone: 2.0 and 1.2).
- **Memory: the index alone** (2026-10-23; commits 0eceb42, 82b18ec;
  `bench/results/hg002_2026-10/memory_{galaxy,a2}.tsv`,
  `memory_identity.txt`). Peak memory was the whole index file (21.4 GiB,
  all of it resident) plus every page of the input read so far: 6.9 GiB
  per 20 M reads, both mates' files when paired. Two changes, every output
  byte unchanged (30 outputs on galaxy, CPU and A100; 22 on a2):
  - *Input pages are released once parsed* (`MADV_DONTNEED` on the
    private, read-only mapping of the FASTQ: a page read again would come
    back from the file). Memory no longer grows with the input.
  - *The bucket directory is held as 32-bit offsets* (1 GiB instead of 2):
    every offset is at most the number of entries, < 2^32. The file keeps
    format v3, so old indexes load and `certa check-index` (Proposition V0)
    applies unchanged; an offset ≥ 2^32 is rejected on load. The load
    converts it on 16 threads into an anonymous mapping (huge pages
    requested) while one thread faults in the other arrays: 0.7 s on
    galaxy as before; 0.2 s on a2, where the old load took 0.04 s (its
    page cache maps the file in large folios).
  - *20 M reads / pairs, interval mode, medians of three interleaved with
    the old build:* peak 29.1 → 21.3 GiB single-end and 36.2 → 21.6 GiB
    paired on both hosts (A100: 31.9 → 24.3 and 39.5 → 25.0 GiB). galaxy is
    2.7–6.5 % faster in every configuration (single-end default 15.79 →
    14.77 s, paired `--full-cert` 99.93 → 96.36 s), a2 0.3–1.2 % faster;
    at 2 M, a2 is 0.7–4 % slower (the 0.2 s load).
  - *Tried and not kept: 8-bit keys* (the keys' high 8 bits held, the low 8
    read from the reference at the entry). Byte-identical and 2.7 GiB
    less, but 21–79 % slower at 2 M: in repeats every binary-search step
    then needs two dependent random reads (the position, then the
    reference), and the lookups of tier R and `find_copies` took 37 % of a
    single-end interval run.
  - *Tried and not kept: a 2-bit reference* (`--packed-ref`; N runs kept
    apart, windows unpacked by table and lazily). Byte-identical and
    2.2 GiB less, but pass V and the pair passes compare many short-lived
    bands, which the byte reference serves by vector loads in place: at
    2 M, `--full-cert` 9 % slower on galaxy and 17–19 % on a2, paired
    default 8 % / 13 % more CPU; even switched off it cost 0.5–0.8 % CPU.
  - What is left is the index itself: positions 11.0 GiB, keys 5.5 GiB,
    reference 2.9 GiB (1 byte per base), directory 1 GiB. Below that would
    take a compressed index (FM-index or sampled positions), which changes
    what each lookup costs rather than how the arrays are stored.
- **HG001 GATK indels** are the one comparison where CERTA's F1 is lower
  (0.9909 vs 0.9922). Look at the 13 extra FP and 14 extra FN.
- **Long reads:** the mapper's limits (320 bp, 8 parts, radius 5) exclude
  them. `tools/lr_ceiling --certify` certifies them (Aim 6, below): OPT is
  proven for 93.3 % of real HiFi reads, at 0.86× minimap2's CPU.
  - *Reads that span SVs:* a split (jump) objective is defined, and Lemma
    L-J bounds it (`certa/longread.h`). With it, 94.4 % of reads fit at
    P = 30, about a fifth of the clipped and high-loss misses.
  - *The exact certificate of the split optimum* (`certify_jump`) is
    oracle-tested. On the reads that fit only the split objective it proves
    OPT_J for 1,551 (0.64 % of all reads), so 93.96 % are now certified
    exact.
  - *Its pruning* is now iterative: exact prefix and suffix values of the
    chains evaluated first, read windows, lb raised to the best alignment
    found, and a static step when rounds cost more. That cut the median from
    45 s to 3.4 s per read and the CPU in half. The first version had used a
    single loose bound. The vectorized DP row then took the median to 1.4 s
    (39 % of the first version's CPU).
  - *Still open:* reads with a long lossy stretch (for example an insertion
    of repeat sequence) stay costly and keep wide intervals; and building
    all this into the mapper with a faster DP.

  A full long-read CERTA would likely be slower than shmap-rs, which does
  no base-level alignment: about 11.7 s for 2,000 HG002 HiFi reads
  single-threaded (a2:`~/shmap-rs/profiling/other-mappers/`). A more
  promising idea is to certify only the reads shmap-rs flags as ambiguous.
- **Data on the hosts:**
  - galaxy:
    - HG005 FASTQs (52 GB, `~/certa-data/SRR14724528_*`); delete them if
      the space is needed;
    - 20M/2M HG002 subsamples with R2 mates (`~/certa-bench/reads`);
    - per-sample results in `~/certa-bench/acc_HG00{1,5}` and `acc2`.
  - a2: the full HG002 R2 (`~/certa-data/SRR37356339_2.fastq.gz`).
