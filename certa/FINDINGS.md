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
  BWA-MEM-style pair MAPQ. S2 is not a uniqueness proof: a second, equal
  partner inside the winning anchor's own window is bounded only by R1.
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
    the floor, band pad or window threshold, survive. They need pairs
    exactly at those thresholds, for example a ~46-bp deletion against a
    decoy one point below.
- **Real data:**
  - 1,672,120 certified HG002 reads were compared with BWA-MEM2's unclipped
    alignments. BWA-MEM2 never has fewer edits than the certified minimum,
    except 3 reads whose windows contain a reference `N` (W4).
  - Output is byte-identical on ARM64, x86-64 and an A100.
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
- **Speed vs the fastest tool is modest.** The certified tier is fast, but
  the uncertified ~5.6 % dominate the remaining cost (`bench/RESULTS.md`
  §5b). The speedup is 4.35× against BWA-MEM2 and only 1.16× against
  minibwa. Answer: fewer uncertified reads (paired-end, a certified
  local-alignment mode) or a faster engine for them.
- **The GPU is part of the claim.** The comparison is A100 + 64 threads vs
  64 threads. Cost per genome (GPU hour vs CPU hours) should be reported
  alongside wall time.
- **Precision is slightly lower than BWA-MEM2's** (GATK SNV 0.9886 vs
  0.9898), with F1 higher through recall. A multi-sample study is needed to
  show this is systematic and not sample noise.
