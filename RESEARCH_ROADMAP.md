# Research Roadmap: Audit, Future Directions, a Proposed Adaptive Aligner, and a Clinical Benchmark Protocol

*Prepared September 2026 as a companion to [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md).
It has four parts. Part 1 audits the literature review. Part 2 sets out
the upgraded future directions. Part 3 proposes a new adaptive aligner
for short and long reads. Part 4 specifies a pre-registrable benchmark
that can support (or refute) the claim that the proposed aligner is
faster than current tools in common clinical scenarios without losing
accuracy.*

> **Status of the claims in this document.** Part 3 is a *design*
> proposal, and every performance statement about it is a *hypothesis*
> to be tested in Part 4, not a result. No benchmark has been run yet.

## Contents

1. [Audit of the existing review](#1-audit-of-the-existing-review)
2. [Future directions for performance and correctness](#2-future-directions-for-performance-and-correctness)
3. [Proposed algorithm: CERTA, a certified, error-tiered adaptive aligner](#3-proposed-algorithm-certa-a-certified-error-tiered-adaptive-aligner)
4. [Benchmark protocol](#4-benchmark-protocol)
5. [Realistic expectations and risks](#5-realistic-expectations-and-risks)
6. [References cited only in this document](#6-references-cited-only-in-this-document)

---

## 1. Audit of the existing review

### 1.1 Verdict

**Was the review complete?** No. The first version described itself as
"cross-checked … to confirm no major algorithm family was missing", but
it omitted several methods that a 2026 PhD committee would expect,
including the two strongest current baselines for a "faster aligner"
claim:

| Omission | Why it matters |
|---|---|
| **DRAGEN** (Behera et al., *Nat. Biotechnol.* 2024) | De facto clinical standard for Illumina data (FPGA; ~30 min FASTQ→VCF for a 30x genome; multigenome graph reference). Any clinical speed claim has to address it. |
| **minibwa** (Li & Homer, arXiv June 2026) | Reported >2x faster than BWA-MEM2 at comparable accuracy. It is now the strongest open CPU short-read baseline. |
| **Strobealign with multi-context seeds** (Tolstoganov et al., *Genome Biol.* 2026) | Fixed strobealign's short-read accuracy gap; very fast. |
| **mapquik** (Ekim et al., *Genome Res.* 2023) | Reported ~30x faster than minimap2 for accurate long reads, using minimizer-space seeds. |
| **A\*PA / A\*PA2** (Groot Koerkamp 2024) | Exact alignment competitive with *approximate* methods; the main exact competitor to WFA. |
| **Block Aligner** (Liu & Steinegger 2023) | Adaptive-effort SIMD DP, directly relevant to an adaptive aligner. |
| **Parabricks, Sentieon** | GPU and CPU clinical production baselines. |
| **mm2-fast, mm2-gb, Winnowmap2, ERT, BLEND, search schemes, syncmers/mod-minimizers, move-structure/Movi, RawHash2/UNCALLED, Minigraph-Cactus, long-read Giraffe** | Each is a significant family or a state-of-the-art point that was missing. |

All of these are now covered in [SWOT_ANALYSIS.md §K](SWOT_ANALYSIS.md#k-recent-and-previously-omitted-methods-20112026).
The review is still deliberately *selective*: it covers DNA
read-to-reference alignment and excludes RNA-seq splice alignment,
protein search and metagenomic classification.

**Was it scientifically correct?** Mostly, but it contained factual
errors that a committee would catch, plus a recurring overclaim. All of
them are now corrected in the main documents:

| # | Original claim | Correction |
|---|---|---|
| 1 | WFA is "provably near-linear" / "near-linear exact" | O(ns+s²) is linear in n only for *bounded* s. At a fixed error rate e, s = Θ(e·n), so time is Θ(e·n²): a large constant-factor gain, still quadratic. |
| 2 | WFA is "adopted inside minimap2" | minimap2 uses KSW2 (Suzuki-Kasahara difference recurrence), not WFA. WFA is used in wfmash, PGGB and vg. |
| 3 | BiWFA paper first author "Eizenga" | Marco-Sola, Eizenga, Guarracino, Paten, Garrison & Moreto (2023). |
| 4 | BiWFA costs "roughly 2x runtime" | Same asymptotic time, and often competitive in practice. |
| 5 | WFA-GPU first author "Marco-Sola" | Aguado-Puig et al. (2023). |
| 6 | Myers' bit-vector is "the extension kernel inside GenASM" | GenASM uses modified Bitap (Wu-Manber); BWA-MEM uses affine-gap SW. |
| 7 | "Graph alignment is NP-hard in general" (twice, including the conclusions) | Sequence-to-graph alignment is O(\|E\|·m) (Navarro 2000). It is NP-hard only when the graph may also be edited (Jain et al. 2020); there is a conditional quadratic lower bound (Equi et al. 2019). |
| 8 | minimap2 chaining is O(N log N) "with a Fenwick tree" | minimap2's default is an O(N·h) bounded-look-back heuristic; O(N log N) needs RMQ-based exact chaining. |
| 9 | HISAT2 uses the minimizer + chaining architecture | HISAT2 is hierarchical-FM-index based. |
| 10 | GATK GPU paper (BMC Genomics 2019) accelerates PairHMM | It accelerates semi-global **Smith-Waterman with traceback**. Authors were missing (Ren, Ahmed, Bertels, Al-Ars). |
| 11 | GPU BWA-MEM (Pham et al.) is from 2024 and accelerates only extension | ICS '23 (2023), and it ports the whole pipeline. |
| 12 | BWA-MEME "provably" preserves output; memory not mentioned | Output is identical *by construction* (bounded last-mile search) and checked empirically, not formally proven. Memory is ~38–118 GB, far above BWA-MEM2. |
| 13 | AnySeq/GPU targets tensor cores; cited as arXiv only | Uses warp shuffles, not tensor cores; published at ICS '22. |
| 14 | Winnowmap uses "learned" minimizers | *Weighted* (frequency-aware) minimizers. |
| 15 | Darwin ">39,000x energy efficiency"; GACT's trade-off understated | Unverified figure removed. GACT is a tiled *heuristic* with no global optimality guarantee. |
| 16 | SneakySnake "up to 979x"; GraphAligner "13x faster, 3x less memory" | Unverified headline numbers replaced with baseline-specific wording. |
| 17 | SETH section: real data "essentially never" adversarial | Repeats, segmental duplications and STRs behave adversarially, and many clinically important genes sit in them. Four-Russians, approximation and LCS/DTW/graph hardness results were added. |
| 18 | Conclusion: WFA changes exact alignment "for the first time in ~40 years" | WFA generalizes 1985–89 diagonal-transition algorithms to affine gaps. |
| 19 | Minor bibliographic errors | SOAP3-dp author list, GASAL2 author order, MaxSSmap availability, FM-index count vs. locate cost. |

**Were the conclusions sound?** The central thesis still holds:
SETH-hardness pushes the field toward parameterized, filtered,
heuristic and hardware-accelerated methods. Three conclusions were
wrong or incomplete, and were rewritten:

1. The NP-hardness claim was false (row 7 above).
2. The "near-linear" framing overstated the case (row 1).
3. The conclusions left out the two facts that matter most for the
   thesis goal. First, **Amdahl's law**: alignment is only part of
   FASTQ→VCF time. Second, **clinical incumbents are hardware-accelerated
   and validated**, so a new method must show accuracy non-inferiority
   and not just speed.

**Verification performed.** All 79 DOIs in `SWOT_ANALYSIS.md` were
resolved against Crossref/doi.org (Sept. 2026), and title, first author
and year were checked. The Python test suite passes. (It used to
crash on Windows cp1252 consoles when printing ✓/✗; it now switches
stdout to UTF-8.)

---

## 2. Future directions for performance and correctness

This replaces the short "Opportunities" list with directions ordered by
their expected value for a clinical-speed thesis. Each direction pairs a
**performance** lever with the **correctness** property that must be
preserved.

### 2.1 Spend effort only where the data needs it (adaptive alignment)

- *Observation.* Modern reads are much more accurate than the reads
  most tools were tuned for (XLEAP-SBS on NovaSeq X; ONT R10.4.1 with
  v5 basecalling at Q20+; HiFi). In the error and variant model below,
  most 150 bp reads differ from the reference by only a few edits. This
  fraction must be **measured** on real NovaSeq X data (Aim 1 in §3.8).
- *Existing partial solutions:* Block Aligner (adaptive band), minibwa
  (less effort in repeats), strobealign (read-length-dependent seed
  parameters), DRAGEN (tiered hardware pipeline), BWA-MEM's heuristic
  early exits.
- *Open problem:* treating engine choice as a formal **algorithm-selection
  problem** (Rice, 1976). That means an explicit cost model, a
  guaranteed fall-back, and a proof that the fast path cannot lower
  accuracy beyond a stated bound.

### 2.2 Lossless filtering for the high-accuracy regime (correctness)

- For reads with ≤k differences, the pigeonhole principle (split into
  k+1 parts, so at least one part matches exactly) and optimum *search
  schemes* over a bidirectional FM-index (Kucherov et al. 2016; Kianfar
  et al. 2018) **provably enumerate every locus within k edits**.
- Most production mappers use lossy seeding for every read. For a
  high-accuracy read, a certified search can be both faster (little
  chaining or DP needed) *and* more rigorous than lossy seeding. This is
  the main theoretical lever in §3.

### 2.3 Make the chaining layer rigorous

- Chaining sits between seeds and base-level DP. It is the least formally
  guaranteed layer and a major source of mis-mappings in repeats.
  Exact chaining with gap costs (Jain et al. 2022) and RMQ-based
  O(N log N) chaining can replace bounded-look-back heuristics. The
  research question is whether they can do so at equal speed.

### 2.4 Calibrated mapping quality (probabilistic correctness)

- MAPQ is used as a probability by every downstream caller, but it is
  rarely validated for *calibration*. Future work: reliability diagrams
  and Brier scores per platform; recalibration models trained on
  simulated reads from HPRC assemblies; and, as a novel direction,
  **conformal prediction** to bound the mis-mapping rate among reads
  above a MAPQ threshold, with a distribution-free guarantee.

### 2.5 Memory latency, not arithmetic, bounds seeding

- FM-index and hash lookups are cache-miss bound. Promising levers are
  prefetching and batching (minibwa, BWA-MEM2), radix-tree seeding
  (ERT), learned indexes with smaller memory than BWA-MEME,
  cache-efficient r-index variants (Movi), and huge pages/NUMA-aware
  index placement. Correctness is easy to keep here, because these are
  exact data-structure substitutions that can be verified by differential
  testing.

### 2.6 Beat Amdahl's law: fuse the pipeline

- Sorting, duplicate marking, CRAM/BAM compression and I/O take a large
  share of the end-to-end time of a CPU pipeline. Streaming the mapper
  directly into coordinate-bucketed sorting, duplicate marking and
  CRAM output is how DRAGEN and Parabricks gain much of their speed. An
  open-source CPU implementation of a fully fused path is a practical,
  measurable contribution.

### 2.7 Reference choice and reference bias

- GRCh38 with ALT masking stays the clinical coordinate system. Pangenome
  references (DRAGEN multigenome, Giraffe with haplotype sampling,
  Minigraph-Cactus/HPRC) reduce reference bias, especially in
  medically relevant difficult regions (CMRG). The open problem is how to
  get pangenome accuracy while still reporting in linear coordinates that
  clinicians and ClinVar use.

### 2.8 Hardware portability

- GPU for irregular long-read chaining (mm2-gb); WFA-GPU; AVX-512/AVX10
  and ARM SVE2 (Graviton/Grace-class CPUs); FPGA (DRAGEN,
  minimap2-fpga); processing-in-memory (SAFARI group). DSL/compiler
  approaches (AnySeq, FILTR 2026) reduce the cost of re-targeting
  kernels. That matters because hardware-specific tools are abandoned
  quickly (NVBIO).

### 2.9 Real-time and raw-signal mapping (ONT)

- For adaptive sampling (RawHash2, UNCALLED, Movi), the relevant metric
  is *decision latency per read*, not throughput. This is a distinct
  clinical use case (targeted repeat-expansion or pharmacogene panels on
  a single flow cell).

### 2.10 Difficult, medically relevant regions

- Repeat expansions (e.g., *HTT*, *FMR1*, *C9orf72*, *RFC1*),
  *SMN1/SMN2*, *CYP2D6*, HLA, *PMS2*/*PMS2CL*, *GBA1*/*GBAP1*: this is
  where heuristic seeding fails and where clinical value concentrates.
  Any "faster" aligner must be reported *stratified* by these regions
  (GIAB stratifications, CMRG).

### 2.11 Verification and reproducibility as first-class goals

- Differential testing of every fast kernel against exact DP (as this
  repository already does for its samples), property-based testing,
  and bit-identical output across thread counts and runs. Deterministic
  tie-breaking is a regulatory requirement in practice (ISO 15189
  laboratories; EU IVDR for in-house software). Formal verification of
  small SIMD kernels is feasible and largely unexplored.

### 2.12 Better benchmarking science

- Most published speedups are component-level, measured on the authors'
  own hardware and data. Standardized, pre-registered, end-to-end,
  hardware-matched benchmarks, reported with confidence intervals,
  would benefit the field and are a publishable contribution in their
  own right (Part 4).

---

## 3. Proposed algorithm: CERTA, a certified, error-tiered adaptive aligner

*CERTA is a working name.*

### 3.1 Novelty statement

Individual ingredients exist: lossless filters (Yara, RazerS3, GEM,
Columba), minimizer-space seeds (mapquik), adaptive DP (Block Aligner),
WFA, repeat-aware effort reduction (minibwa), and tiered pipelines
(DRAGEN). A defensible PhD contribution is therefore **not** "a new seed
type". It is three things:

1. A **per-read (short reads) and per-gap (long reads) algorithm-selection
   framework** with an explicit cost model and a *correctness
   certificate* for every fast-path decision.
2. A proof that the fast path returns the same primary alignment as an
   exhaustive ≤k-error search (short reads), or the optimal gap-affine
   alignment between anchors (long reads).
3. An **end-to-end, fused** implementation benchmarked under the
   protocol in Part 4.

### 3.2 Architecture overview

```
FASTQ/uBAM ──► Profiler (first ~10^5 reads, seconds)
                 │  platform, read-length dist., est. error rate,
                 │  PE/SE, insert size, WGS/WES/panel, target BED
                 ▼
             Policy (deterministic, trained offline, version-locked)
        ┌────────┴──────────────┐
        ▼                       ▼
  Short-read engine        Long-read engine
  (Tier S0–S3)             (Tier L0–L3)
        └────────┬──────────────┘
                 ▼
  Fused output: MAPQ calibration → coordinate-bucketed sort
                → streaming duplicate marking → CRAM/BAM
```

### 3.3 Short-read engine (Illumina and other SBS platforms)

| Tier | Trigger | Method | Guarantee |
|---|---|---|---|
| **S0 exact** | Whole read (or both mates) found by a direct hash of its k-mers or by an FM-index count; ≤ c_max loci | Emit directly (CIGAR `150M`), no DP | Exact; all loci within 0 edits enumerated |
| **S1 certified ≤k** | S0 fails; k = 2–4, chosen by the profiler from the estimated error rate | Pigeonhole partition into k+1 parts (or an optimum search scheme on a bidirectional FM-index) → candidate loci → bit-parallel Hamming check, then Myers/banded verification; mate-pair consistency | **Lossless**: every locus within k edits is found, so the best alignment and the second-best gap are exact. That makes MAPQ well-founded. |
| **S2 standard** | No locus within k edits, or ambiguity (best − second-best < Δ) | minibwa/BWA-MEM-style SMEM seeding + exact RMQ chaining + banded affine DP (KSW2-style, or WFA with a small predicted s) | Heuristic (the same class as current tools) |
| **S3 rescue / difficult** | Read or mate in a difficult region (precomputed mask: segdups, STRs, CMRG loci), or S2 fails | Higher-sensitivity seeding, mate rescue, optional pangenome (haplotype-sampled) alignment projected to GRCh38 | Heuristic; reported separately |

Key points:

- **The certificate.** In S1 the read is split into k+1 non-overlapping
  parts, so any alignment with ≤k edits contains at least one part
  exactly. Enumerating all exact hits of all parts therefore finds every
  locus within k edits, and verification is exact. The emitted primary
  alignment is then *identical* to what an exhaustive ≤k search would
  return. A full ≤k-edit guarantee needs indel-aware verification
  windows of ±k around each candidate diagonal; this is standard.
- **Cost control.** Parts that hit highly repetitive sequence are
  capped, and a read whose cap is reached *loses its certificate*
  and goes to S2/S3. The certificate is never silently weakened.
- **Expected share of S0/S1.** For a 150 bp read, with per-base
  sequencing error around 0.1–0.3 % and heterozygous SNV density around
  1/1,000 bp, the expected number of differences per read is below one.
  A large share of reads should therefore be S0/S1 eligible. **The exact
  figure is Aim 1's first measurement**; it determines the achievable
  speedup, because by Amdahl's law the S2/S3 residue bounds the gain.
- **Determinism.** Ties are broken by a hash of read name and locus,
  never by thread scheduling.

### 3.4 Long-read engine (ONT R10.4.1, PacBio HiFi)

| Tier | Trigger | Method | Guarantee |
|---|---|---|---|
| **L0 anchor** | Q20+ read; unique minimizer-space anchors (k consecutive minimizers occurring once in the reference) | mapquik-style coarse placement | Anchors are exact and unique |
| **L1 gap fill** | Between consecutive anchors | Per-gap engine selection: identical gap → copy; small predicted cost → **WFA/BiWFA** (exact); larger or divergent → Block-Aligner/KSW2 with Z-drop; long indel → two-piece affine (SV-aware) | Exact per gap when WFA is chosen; the router's lower bound is provable (below) |
| **L2 standard** | Too few unique anchors (repeats, centromeres, segdups) | minimap2-style chaining with weighted minimizers (Winnowmap-like) + exact RMQ chaining + KSW2 | Heuristic (the same class as minimap2) |
| **L3 difficult** | STR/repeat-expansion loci, highly divergent reads | Repeat-aware scoring; hand-off to specialized genotypers (Straglr, TRGT-class tools) | Reported separately |

Key points:

- **Provable routing bound.** For a gap between anchors with query
  length Δq and reference length Δr, any alignment of the gap costs at
  least `o + e·|Δq − Δr|` (gap open plus extend for the net length
  difference). That is a *lower bound on WFA's score s* for the gap.
  Together with a cheap upper bound (e.g., a Hamming or diagonal
  estimate), it bounds WFA's work, O(n·s + s²), before any work is done.
  The router can then pick the cheaper engine with a bounded worst case.
- Gaps between anchors in Q20+ reads are short and highly similar,
  which is WFA's best regime. That keeps the "WFA is quadratic at a
  fixed error rate" problem (§1, row 1) local to each short gap instead
  of the whole read.
- Base-modification tags (MM/ML) and read groups pass through
  unchanged, which is required for ONT clinical methylation workflows.

### 3.5 Adaptive policy (the "adapts to data" part)

- **Features** (cheap, computed before alignment): read length, mean and
  minimum base quality, fraction of k-mers in the reference, number of
  seed hits per part, whether the read falls in a difficult-region mask,
  platform and profile from the profiler.
- **Cost model:** predicted CPU time per tier (from micro-benchmarks) and
  predicted probability that the tier succeeds or certifies.
- **Policy:** minimize expected time subject to (i) never skipping a
  fall-back when a certificate fails, and (ii) a validation-set
  constraint that mapping accuracy is non-inferior to S2/L2 alone.
- **Offline training only.** Thresholds are fitted on training samples
  disjoint from all test samples (§4.4), then frozen and versioned.
  **There is no online learning in clinical mode**; reproducibility
  requires deterministic behaviour for a given version and input.

### 3.6 Correctness framework

1. **Differential testing** of every kernel against exact DP (this
   repository's samples provide the oracle pattern).
2. **Certificate audits:** on simulated reads with known truth, every
   S0/S1 decision is re-checked by exhaustive search.
3. **MAPQ calibration:** reliability curves per tier and per platform
   (§2.4).
4. **Invariance tests:** identical output across thread counts,
   operating systems and repeated runs.
5. **Error budget:** the paper reports the mis-mapping rate by tier, so
   the fast path's contribution to errors is visible, not averaged away.

### 3.7 Hardware mapping

- CPU first (AVX2/AVX-512/AVX10 on x86, SVE2 on ARM). Portable SIMD via
  a thin abstraction such as Highway or SIMDe.
- Optional GPU offload of S2 DP batches and L1 gap batches (batched
  kernels in the GASAL2/WFA-GPU style).
- Compare *within hardware class* (§4.5). A CPU tool is not claimed to
  beat an FPGA on wall time. Cross-class comparison uses cost and energy
  per genome.

### 3.8 Work plan (suggested thesis aims)

| Aim | Deliverable | Go/no-go criterion |
|---|---|---|
| **1. Feasibility** | Measure the share of S0/S1-eligible reads and the per-tier time on NovaSeq X and ONT R10.4.1 GIAB data | If S0+S1 is below ~60 % of reads, the short-read speedup ceiling is low and the thesis should focus on L1 and fusion |
| **2. Short-read engine** | S0–S2 + fused sort/markdup; certificate proofs; differential tests | Non-inferior accuracy on GIAB HG002 (pilot) |
| **3. Long-read engine** | L0–L2 with a WFA gap router; routing lower-bound proof | ≥1.5x faster than minimap2 `lr:hq` on the pilot at non-inferior accuracy |
| **4. Benchmark** | Pre-registered protocol (Part 4), full run, open data and containers | Pre-registered hypotheses tested as stated, whatever the outcome |

### 3.9 Prototype status (September 2026)

Tiers S0/S1 are prototyped in [`certa/`](certa/README.md) for CPU (x86-64,
ARM64) and CUDA (A100 target), with one shared per-read core.

- **Correctness:** a brute-force oracle test (Sellers DP over the whole
  reference) passes for soundness, optimality, losslessness, uniqueness
  and determinism. The GPU output is byte-identical to the CPU output.
- **Synthetic data:** on a 20 Mbp genome with repeats and 150 bp reads,
  95.7 % of reads are certified at k = 2.
- **Lab hosts:** the prototype runs on galaxy (ARM64 CPU and A100) and
  a2 (x86-64). Outputs are byte-identical across all three back-ends; on
  synthetic data the map step reaches 8.85 M reads/s on the A100.
- **Next step:** run the same measurement on real NovaSeq X GIAB data to
  produce Aim 1's go/no-go number.

---

## 4. Benchmark protocol

Designed so that the result is publishable whatever the outcome. The
claim being tested is: *"CERTA is faster than current tools in the most
common clinical scenarios, without loss of accuracy."*

### 4.1 Pre-registered hypotheses

Register on OSF before running the main benchmark.

- **H1 (speed, per scenario × competitor within the same hardware
  class).** The geometric-mean wall-clock ratio (competitor / CERTA) for
  FASTQ → sorted, duplicate-marked CRAM is ≥ 1.10, *and* CERTA is faster
  on a majority of samples.
- **H2 (accuracy non-inferiority).** Per-sample differences, CERTA minus
  competitor, with the variant caller held fixed, have a 97.5 %
  one-sided lower confidence bound above −Δ. Margins:
  - SNV F1: Δ = 0.001
  - indel F1: Δ = 0.005
  - SV F1 (long reads): Δ = 0.01
  - CMRG-gene SNV+indel F1: Δ = 0.005
  - Mendelian-violation rate on trios: +0.05 % absolute
- **Confirmatory claim:** both H1 and H2 must hold. A speed win with an
  accuracy loss is reported as a trade-off, not a win.
- **Exploratory:** mapping-only time, peak memory, cost and energy per
  genome, MAPQ calibration, stratified accuracy.

### 4.2 Clinical scenarios

These are the "most common clinical cases".

| ID | Scenario | Platform | Typical data |
|---|---|---|---|
| S1 | Germline WGS, rare disease | Illumina NovaSeq X (XLEAP-SBS), PE150 | 30–40x |
| S2 | Exome / large panel | NovaSeq X or NextSeq 2000, PE150 | 100–200x on target |
| S3 | Small targeted panel (e.g., hereditary cancer) | MiSeq i100 / NextSeq | ≥500x on target |
| S4 | Somatic tumor/normal | NovaSeq X | e.g., 60–100x tumor / 30–40x normal |
| S5 | Long-read germline WGS (SVs, repeats, phasing, methylation) | ONT PromethION, R10.4.1, Dorado v5 SUP/HAC | ~30x |
| S6 | Long-read targeted / adaptive sampling | ONT R10.4.1 | Panel; latency metric |
| S7 | Rapid (STAT/NICU) WGS | NovaSeq X or ONT | Time to first VCF |

Report every scenario, including those where CERTA loses. "Most" means
a pre-declared majority of S1–S5, with S1, S2 and S5 as primary.

### 4.3 Datasets

**With truth sets (accuracy).**

- **GIAB HG001–HG007:** v4.2.1 small-variant truth (GRCh38) and the
  **HG002 Q100 v5.0q** benchmark (T2T-HG002-based, covering more
  difficult regions). Also **CMRG v1.00** (Wagner et al. 2022), GIAB
  **stratifications** (difficult regions), and HG002 SV truth.
- **GIAB HG008** (tumor/normal) and SEQC2 HCC1395/HCC1395BL (Fang et al.
  2021) for S4.
- **ONT:** `s3://ont-open-data/giab_2025.01/`, with all seven GIAB
  samples on R10.4.1, Dorado v5 basecalling and two PromethION flow
  cells each (≈140 Gb per flow cell, per ONT).
- **Illumina NovaSeq X:** public GIAB data from XLEAP-SBS chemistry is
  scarce; at the time of writing, only HG002 runs appear to be public.
  **Recommendation: sequence the NIST reference materials (RM 8391
  HG002 trio, RM 8392 HG005 trio, RM 8398 HG001) in-house on NovaSeq X.**
  Use ≥3 independent libraries per sample, including PCR-free WGS and one
  exome kit. This is the cheapest way to get latest-chemistry data with
  truth, and it provides the replicate structure needed in §4.6.
- **Platinum Pedigree (CEPH 1463)**: 28 members over four generations,
  with Illumina, ONT and HiFi data (AWS Open Data). It supports
  Mendelian-consistency evaluation on many related genomes (Porubsky et
  al. 2025; Kronenberg et al. 2025). Check chemistry and basecaller
  version per sample and keep only R10.4.1 in the primary ONT analysis.

**Without truth (speed, scalability, Mendelian checks).**

- **1000 Genomes high-coverage** (NYGC; Byrska-Bishop et al. 2022),
  3,202 samples including 602 trios. The data are NovaSeq 6000, older
  chemistry, so use them for ancestry diversity and runtime scaling, not
  as "latest hardware".
- **Clinical cohort (strongest evidence).** Use ≥30–50 de-identified
  clinical WGS/WES samples from a partner laboratory, with orthogonally
  confirmed reportable variants (Sanger or MLPA) as truth for clinically
  relevant calls. This requires ethics-committee approval and a GDPR
  data-processing agreement. Run on-premises; the data do not leave
  the laboratory.

**Simulated reads (read-level truth).** Simulate from HPRC diploid
assemblies (not from the reference, to avoid optimistic bias) using
platform-matched error models: a NovaSeq X profile for Illumina, and a
simulator validated on R10.4.1 for ONT (see the 2026 comparison of ONT
simulators). This measures the correct-placement rate against MAPQ, and
certificate audits.

**Data leakage rule.** Samples used to fit the policy (§3.5), for example
HG001, HG005 and a few 1000G samples, are **excluded** from confirmatory
testing. Declare the split in the pre-registration.

### 4.4 Competitors

Pin versions and container digests.

| Class | Short reads | Long reads |
|---|---|---|
| CPU, open source | BWA-MEM2 (clinical reference point), **minibwa**, **strobealign (MCS)**, Bowtie2, BWA-MEME (if memory allows) | **minimap2** (`lr:hq`, `map-ont`, `map-hifi`), mm2-fast, Winnowmap2, **mapquik** (HiFi/Q20+), minibwa (HiFi) |
| CPU, commercial | Sentieon BWA (license) | — |
| GPU | Parabricks `fq2bam` | mm2-gb; Parabricks minimap2 (if available) |
| FPGA | DRAGEN (on-premises server or cloud FPGA instance) | — |
| Graph | vg Giraffe (HPRC, haplotype-sampled) | vg Giraffe long-read |

Every tool uses its **recommended preset** for the data type and the
**same reference build** (GRCh38 analysis set, ALT-masked or ALT-aware
consistently). Every tool gets the **equivalent downstream steps** to
produce sorted, duplicate-marked CRAM (e.g., `bwa-mem2 | samtools sort |
samtools markdup`, compared with Parabricks/DRAGEN's fused equivalents).
Ask the authors of the leading competitors to review your
configurations. This is the best defence against "you ran our tool
wrong".

### 4.5 Hardware tiers

Comparisons of wall time are made **only within a tier**. Across tiers,
report cost per genome (on-demand cloud price, stated date) and energy
per genome.

- **CPU tier:** one current x86 server (e.g., AVX-512-capable, 64–128
  cores, ≥512 GB RAM, NVMe) and one ARM server (SVE2). Run the same
  thread counts (8, 32, all cores) for every tool.
- **GPU tier:** one current data-centre GPU node, specified exactly.
- **FPGA tier:** a DRAGEN server or cloud FPGA instance.
- **Measurement:** `/usr/bin/time -v` (wall time, CPU time, peak RSS);
  RAPL/IPMI/`nvidia-smi` for energy. Drop the OS page cache before cold
  runs and report both cold and warm runs. Index load time is reported
  separately *and* included in end-to-end time. Isolate the machine
  (no other jobs; turbo policy recorded).

### 4.6 Statistical design and sample size

- **Unit of analysis = sample (library).** Repeated runs of the same
  sample (≥5 per configuration, interleaved in random order) measure
  machine noise. Use their median; they are *not* extra samples.
- **Speed (H1).** For each sample, compute the log wall-time ratio.
  Test with a one-sided paired t-test on the logs (or Wilcoxon
  signed-rank as a sensitivity analysis) plus a sign test for "faster on
  the majority". Apply Holm correction across all
  (scenario × competitor) tests and report geometric-mean speedups with
  95 % CIs.
- **Accuracy (H2).** Use a paired, one-sided non-inferiority test on the
  per-sample F1 differences. When GIAB replicates of the same individual
  are pooled, use a **mixed-effects model with individual as a random
  effect** (replicates are not independent in genome content). Report
  within-sample uncertainty with a block bootstrap over genomic regions.
- **Sample-size calculation** (normal approximation with t-correction,
  power 0.9; to be re-run with the SDs from the pilot):

| Quantity | Assumed SD of paired difference | α (one-sided) | n needed |
|---|---|---|---|
| ≥10 % speedup (log ratio) | 0.05 / 0.10 / 0.15 / 0.20 | 0.05 | 4 / 11 / 23 / 40 |
| same, Bonferroni over 24 tests | 0.05 / 0.10 / 0.15 / 0.20 | 0.0021 | 9 / 24 / 47 / 80 |
| SNV F1 non-inferiority, Δ = 0.001 | 0.0010 / 0.0015 | 0.025 | 13 / 26 |
| indel F1 non-inferiority, Δ = 0.005 | 0.004 / 0.006 | 0.025 | 9 / 18 |
| "faster on majority" sign test (Bonferroni 24), true share 0.9 / 0.8 / 0.7 | — | 0.0021 | 23 / 44 / 107 |

- **Recommendation.**
  1. Run a **pilot** (n ≈ 8 per primary scenario) to estimate SDs.
  2. Plan for **≥30 samples per primary scenario** for speed, raised to
     ≥50 if the pilot SD is ≥0.15.
  3. For truth-based accuracy, use **≥24 truth-set libraries per
     platform** (7 GIAB individuals × 3–4 independent libraries or runs),
     plus Mendelian consistency on trios and the Platinum Pedigree.
  4. The **headline "faster in most cases" claim** needs ≥44 samples if
     the true win rate is ~0.8.

### 4.7 Metrics

- **Speed:** end-to-end wall time; mapping-only wall time; CPU-hours;
  throughput (reads/s or Gbp/h); scaling with thread count; peak memory;
  index size and build time; cost and energy per genome; for S6/S7,
  time to first result and per-read decision latency.
- **Accuracy:**
  - SNV/indel precision, recall and F1 with `hap.py` / `vcfeval` (GA4GH
    best practice; Krusche et al. 2019), overall and by GIAB
    stratification.
  - CMRG genes.
  - SVs with Truvari; tandem repeats and repeat expansions against
    benchmark STR sets.
  - Mendelian violation rate.
  - Concordance on clinically reportable (ClinVar P/LP) sites.
  - Read level: correct-placement rate vs. MAPQ (simulated reads);
    calibration (Brier score, reliability curve).
- **Variant callers held fixed** to isolate the aligner's effect:
  DeepVariant for Illumina, Clair3 for ONT, Sniffles2 for SVs. A
  separate, clearly labelled end-to-end comparison uses each vendor's
  own caller (DRAGEN, Parabricks).

### 4.8 Threats to validity and mitigations

| Threat | Mitigation |
|---|---|
| Tuning on test data | Frozen policy; disjoint training split; pre-registration |
| Mis-configured competitors | Recommended presets; author review; configurations published |
| Hardware or noise artefacts | Isolated nodes; randomized run order; ≥5 repeats; cold and warm runs |
| Truth-set bias toward easy regions | Q100 v5.0q and CMRG; stratified reporting |
| Simulation optimism | Simulate from HPRC assemblies with platform-validated error models; the real-data results are the primary ones |
| Selective reporting | Report all scenarios, competitors and failures; open containers, scripts and raw timing logs (Zenodo DOI) |
| Version drift (tools evolve during the PhD) | Freeze a version snapshot at pre-registration; add a later "current versions" sensitivity run |

### 4.9 Reproducibility package

Include:

- Snakemake or Nextflow workflow
- Apptainer/Docker images pinned by digest
- Reference files with checksums
- A sample manifest with accessions
- Raw `time -v` logs
- `hap.py` and Truvari outputs
- An analysis notebook that regenerates every table and figure

---

## 5. Realistic expectations and risks

- **Against DRAGEN or Parabricks on wall time:** a CPU-only aligner is
  unlikely to win within their hardware tier. Compete on
  cost and energy per genome, openness and portability, and certified
  correctness. Do not claim a wall-time win across tiers.
- **Against minibwa and strobealign (CPU, short reads):** both are highly
  optimized, so any gain must come from the S0/S1 fast path share and
  from pipeline fusion. Aim 1's measurement decides whether this is
  realistic. If the S2/S3 residue is large, Amdahl's law caps the gain.
- **Long reads (ONT R10.4.1 Q20+ and HiFi):** this is the most promising
  target. mapquik shows that seeding can be around 30x faster, and
  base-level DP dominates the remaining time. WFA gap filling between
  unique anchors attacks exactly that DP time.
- **End-to-end CPU pipelines:** fusing sort, duplicate marking and CRAM
  output is a low-risk, measurable win against the standard
  `aligner | samtools` chains that many laboratories still run.
- **Main scientific risk:** a speed win bought with accuracy loss in
  difficult regions. The stratified non-inferiority design in §4 is
  built to detect this and report it honestly.

---

## 6. References cited only in this document

(References for all methods are in [SWOT_ANALYSIS.md](SWOT_ANALYSIS.md).)

- Rice, J. R. (1976). The algorithm selection problem. *Advances in
  Computers*, 15, 65-118. https://doi.org/10.1016/S0065-2458(08)60520-3
- Krusche, P., et al. (2019). Best practices for benchmarking germline
  small-variant calls in human genomes. *Nature Biotechnology*, 37,
  555-560. https://doi.org/10.1038/s41587-019-0054-x
- Wagner, J., et al. (2022). Curated variation benchmarks for challenging
  medically relevant autosomal genes. *Nature Biotechnology*, 40,
  672-680. https://doi.org/10.1038/s41587-021-01158-1
- Byrska-Bishop, M., et al. (2022). High-coverage whole-genome sequencing
  of the expanded 1000 Genomes Project cohort including 602 trios.
  *Cell*, 185(18), 3426-3440. https://doi.org/10.1016/j.cell.2022.08.004
- Porubsky, D., et al. (2025). Human de novo mutation rates from a
  four-generation pedigree reference. *Nature*, 643, 427-436.
  https://doi.org/10.1038/s41586-025-08922-2
- Kronenberg, Z., et al. (2025). The Platinum Pedigree: a long-read
  benchmark for genetic variants. *Nature Methods*, 22, 1669-1676.
  https://doi.org/10.1038/s41592-025-02750-y
- Fang, L. T., et al. (2021). Establishing community reference samples,
  data and call sets for benchmarking cancer mutations using whole-genome
  sequencing. *Nature Biotechnology*, 39, 1151-1160.
  https://doi.org/10.1038/s41587-021-00993-6
- English, A. C., Menon, V. K., Gibbs, R. A., Metcalf, G. A., & Sedlazeck,
  F. J. (2022). Truvari: refined structural variant comparison preserves
  allelic diversity. *Genome Biology*, 23, 271.
  https://doi.org/10.1186/s13059-022-02840-6
- ONT Open Data, GIAB 2025.01 release (R10.4.1, Dorado v5):
  https://epi2me.nanoporetech.com/giab-2025.01/
- NIST Genome in a Bottle (v4.2.1, HG002 Q100 v5.0q, stratifications,
  HG008): https://www.nist.gov/programs-projects/genome-bottle
