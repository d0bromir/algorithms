# Research guideline: certified, biologically-correct, fast read placement

**Status:** living guideline for the PhD R&D programme on CERTA.
**Date:** 2026-11-21. **Companion documents:** [`THEORY.md`](../THEORY.md),
[`FINDINGS.md`](../FINDINGS.md), [`Biologically correct read placement ideas.md`](Biologically%20correct%20read%20placement%20ideas.md),
[`Adaptive certified read alignment.md`](Adaptive%20certified%20read%20alignment.md),
[`Theoretical limits of certified read placement.md`](Theoretical%20limits%20of%20certified%20read%20placement.md),
[`../../RESEARCH_ROADMAP.md`](../../RESEARCH_ROADMAP.md), [`../../SWOT_ANALYSIS.md`](../../SWOT_ANALYSIS.md)
(and its correctness review [`../../SWOT_REVIEW.md`](../../SWOT_REVIEW.md)).

This document consolidates the repo analysis of 2026-11-21 into a single
forward plan: what to change, where it lands in the code, what it should
buy, how to validate it, and how to run it on the lab hosts without
filling their disks. It is the entry point for implementation; each step
below is meant to be executed one at a time, with a test and a host
confirmation, after the guideline PR is merged.

---

## 1. Goal and the correct framing of "the gap"

**Goal (stated 2026-11-21):** a FASTQ-to-human-reference aligner that is
biologically correct to the highest possible extent, fast, and with better
read-placement accuracy than any known aligner.

**The repo already meets the comparative part of this goal.** On the HG005
holdout (the latest refit, PR #43), misplaced *confident* mates per million:

| run | N=99.75% | N=99.29% | N=98% | N=95% |
|---|---|---|---|---|
| CERTA `--mapq pair --t2t` | **95.1** | 24.8 | 7.1 | 2.6 |
| CERTA `--mapq pair` | 185.2 | 60.8 | 21.3 | 10.0 |
| minibwa | 208.0 | 79.9 | 37.4 | 37.4 |
| BWA-MEM2 | 202.2 | 107.9 | 107.9 | 107.9 |
| minimap2 | 252.2 | … | … | … |
| strobealign | 322.6 | … | … | … |

Overall at-locus placement (test B): CERTA 99.913% vs minibwa 99.909%,
BWA-MEM2 99.908%. So CERTA is already best-in-class on confident
misplacement (~2.2× fewer than the best competitor at N=99.75%, ~14× fewer
at N=95%).

**The residual error is a specific, hard class — not a search failure.**
The oracle decomposition (`tools/truth_oracle.cpp`, `bench/truth_oracle.py`)
splits every misplacement into *missed* (truth scores higher → search
failure), *tie*, and *reference* (truth scores lower → no reference-only
score-maximiser can place it there). For minibwa, **88.0 of its 88.7
confident errors/M are reference-class**; for CERTA only ~9 of 122/M are
avoidable by another GRCh38 run. These are reads whose true source copy is
**missing from or collapsed in GRCh38** (segmental duplications). For them
the biologically correct output is to **abstain** (low MAPQ), or to name the
true copy against an augmented reference — not to place confidently in
GRCh38. This is a **MAPQ-calibration + reference-content problem**, and it
is the frontier.

**Number reconciliation (do not blend these):**
- `122/M` = interval `--iv-mapq cert`, **HG002** test B, N=99.44% (2026-10-31).
- `95.1/M` = refit `--mapq pair --t2t`, **HG005**, N=99.75% (2026-11-20).
- `--bio` = 138/M, **HG002**, N=95.58%; **`--bio` has never been run on
  HG005** (open "Next" item).
- The "~2 Phred per competitor mismatch" slope and "85% of confident errors
  in the 3% of mates with BWA-MEM2 gap < 50" are **inferences from a marginal
  stratification**, not per-read likelihood ratios; the hs1→hg38 collapse
  prior still needs checking in the chain file. Experiment A1 converts these
  inferences into a measured ceiling.

---

## 2. Decision: evolve CERTA, do not rewrite

Keep and extend CERTA. A from-scratch aligner would discard the certificate
infrastructure, the Lean 4 proofs, the independent `certa check` auditor,
byte-identical x86/ARM/A100 output, and the validated GIAB benchmarks — and
would still hit the same wall, because the dominant residual error is not a
search problem. The only "from scratch" component worth building is the
**copy-level joint assigner** (§3, B5), a downstream tool, not a replacement
aligner.

---

## 3. The plan: three pillars + enabling work

Gains are the placement report's own estimates (overlapping rows do not
add). Baselines: 122/M default, 138/M `--bio` (HG002). Code sites verified
against the tree at commit `7e6de1b`.

### Pillar A — Calibrated confidence (cheapest, largest lever; do first)

| Step | Change | Code site | Gain | Cost | Status |
|---|---|---|---|---|---|
| **A1** | Offline ceiling test: re-rank existing mates by `min(MAPQ, isotonic(AS−XS))` fit on HG002; re-score at matched-N. No mapper code. | `bench/` scripts | sets ceiling (~104/M); go/no-go for B | hours | proposed |
| **A2** | Gap-graded competitor search: relative horizon (best + 10–12 edits) instead of absolute 8; count **clipped local** competitors (BWA's XS is local); or `sub` from seed chains. Certificate untouched. | `src/map_config.cpp:54` (`--mapq-limit`) → `core.h:1911` `finish_certified` → `other_candidates` (`core.h:1336–1353`, `S2_TOP=8`/`S2_MAX=10`, end-to-end `band_distance`); `sub=19` floor `interval.h:66–72` `mapq_se`. SL already uses a local chain sub. | ≥33/M (parity with minibwa's graded 48–56), up to ~104/M | days | proposed |
| **A3a** | `--bio` near-competitor **floor**: `MAPQ = min(posterior, gap-graded MAPQ)`. | `bio.h:1198` (`o.mapq=`); pair `1305–1306`; rescue `1462–1463` | 138/M → ≤ default from floor alone | hours | proposed |
| **A3b** | `--bio` **phantom-copy term** `π_∅(locus)·ℓ_phantom(r)` in the posterior denominator (paralog-shaped alternative at divergence δ≈1e-2, replacing the i.i.d. null). Fit δ, π_∅ on HG002 to reproduce ~2 Phred/mismatch. | `bio.h` `bio_place:1188–1198` (`others`); pair `bio_pair:1291,1305`; rescue `bio_rescue:1435–1447`; or `posterior_mapq:993–999` | closes most of the 138-vs-46 gap | days | proposed |
| **A3c** | MAQ-style **divergence cap** (Q20) in copy-deficit loci, **MAPQ-side only** (do not touch `BioModel::mis_`, which changes the certificate). | `bio.h:1198`. **Caveat:** `BioLocus` stores no per-locus high-Q mismatch count — a count-keyed cap needs new traceback/DP plumbing; a `p1`-keyed cap works with existing fields. | part of A4 | days | proposed |
| **A4** | Divergence × locus prior: penalise high-Q mismatches only where CHM13 has extra copies (chain) or GIAB SD strata, so het reads in unique sequence are not penalised. | `bio.h` + per-locus prior track | ~18/M (gap≥50 tail) | days | proposed |
| **A5** | Leakage-free calibrator: Qtip-style model trained on **CHM13-simulated** reads mapped to GRCh38, CHM13-only reads labelled *unplaceable*. Optionally **conformal prediction** for a distribution-free mis-mapping-rate bound (roadmap-level, `RESEARCH_ROADMAP.md §2.4`; not yet in `certa/`). | `bench/truth_mapq_*.py`, `mapq_table.h`/`mapq_model.h` | calibrates A2–A4 jointly; transfer to HG005 | days | proposed |

**Guardrail:** fit on HG002, **freeze**, score HG005 once. Report the MAPQ
tie-group size at the matched-N boundary (a saturating MAPQ cannot be
operated at that N). Pre-fixed adoption rule: ≥5% lower errors at 2 of 3
operating points and ≤5% rise at any.

### Pillar B — Reference content (the missing copies)

A CHM13 second opinion removes **42–59%** of every aligner's confident
errors; a second GRCh38 aligner removes 0–3%. The information is in the
sequence, not the algorithm.

| Step | Change | Gain | Cost | Status |
|---|---|---|---|---|
| **B1** | Upgrade hs1 → **T2T-CHM13 v2.0**; re-lift with the T2T **nf-LO chain** (`chm13v2-grch38.chain`). Move from features-only/veto to **levioSAM2-style commit/defer/suppress**; measure how often the lifted CHM13 placement lands within 150 bp of truth and commit those. | keeps 42–59% (~51–72/M), fewer false flags (veto precision only 14–23%), some rescues | hours (alignments exist); two indexes | proposed. **Today `--t2t` is features-only** (27-feat trees, `mapq_model_t2t.h`); the veto lives only in `bench/truth_t2t.py`, `lift_sam.py`; lift class in `src/chain.cpp:156` `t2t_view` (cls1 = within 150 bp). Commit/defer/suppress is new code. |
| **B2** | **One-pass CHM13 sink contigs**: CHM13 segments unaligned/collapsed in the chain, **satellites excluded**, **alt-masked**, appended as **ordinary contigs** (ALT-listed sinks would not lower the GRCh38 hit's MAPQ). Sink-winning mate → MAPQ 0 + sink-coordinate tag; certificate becomes a statement over GRCh38 ∪ sinks. **Validate every sink per-locus** (FixItFelix's FANCD2 recall collapsed 0.941→0.318 when a decoy stole real reads). | up to the veto's 42–59% in one pass; **also kills the `--t2t` 2× cost** (no second index) | weeks | proposed |
| **B3** | Reference hygiene: mask GRC false duplications + chrY PARs (`maskedGRC_exclusions_v2`). | small on confident errors; adds confident correct mates in CBS/CRYAA/KCNE1; removes PAR artefacts | hours | proposed |
| **B4** | **HPRC2 / pangenome as copy catalog** (after sinks): CHM13-based graph, diploid haplotype sampling, surjected to GRCh38; used as confidence signal + copy catalog, not a coordinate replacement. First test the ~21,000 unmapped mates as a veto. | uncertain per-read; strong at variant level | weeks + large downloads | proposed |
| **B5** | **Copy-level joint assignment (research frontier):** second-pass **EM over reads × copies** of each SD family (GRCh38 copies + sinks + phantom), so mismatches recurring across reads at the same offsets vote for an unseen copy. Prototype on FixItFelix genes (KMT2C, KCNJ18, MAP2K3, DUSP22). | the only route to naming the true copy in SDs | months | proposed |

**Leakage rule:** HG002 Q100 haplotypes are truth/calibration only, never
sink sequence. **CHM13 is the only donor that keeps both HG002 and HG005
clean.** HG005 haplotypes may donate sinks for HG002 only in a diagnostic
that sacrifices HG005 as a holdout. HPRC eval graphs exclude HG002/HG005/NA19240.

### Pillar C — Speed (fast *and* correct)

Measured costs (galaxy, 64 cores): `--bio` **16–17× default** (211/379 s vs
13.5/22.0 s; ≈1.0× BWA-MEM2 wall, ≈5× minibwa); `--t2t` **1.9–2.2×** CPU/wall,
+21 GiB; default interval mode already ~1.8–2.4× faster than minibwa wall.
Profile (paired interval): tier R **26%**, passes 1–2 22%, CIGAR/pair DP
fills ~22%, window scans 7%. **GPU runs passes 1–2 only.**

| Step | Change | Property / gain |
|---|---|---|
| **C1** | Adaptive per-read routing, classified: (i) *execution-only* (scheduling, CPU/GPU split, cache reordering, seed caching, memoisation) → byte-identical by construction; (ii) *search-effort* (pass budgets, skip/reorder for predicted-uncertifiable reads) → **per-read rules frozen after calibration**, invariant to threads/batches/device/order/sharding; (iii) *certified-problem parameters* (scoring, the fixed 850-bp pair distance) → **never adapt online**. | preserves the determinism moat (BWA-MEM `-K` is the cautionary tale) |
| **C2** | **CERTA-X exact completion** (replace heuristic fallback): **tier F** second index at **q≈15** (10 parts → loss ceiling 45–50 vs 25–30 at q=22; lift `PMAX=8`/`KMAX=5` in `core.h`); **Theorem G** greedy optimal parts; **Theorem R′** window-collapse (DP once per distinct window string → up to 100× budget in near-identical repeats); **Theorem R** count-based tie proofs; budgeted branch-and-bound (R5) as the terminal *bounded* fallback. | raises certification past the q=22 barrier; removes the fallback cost (6.4% of reads = 54% of minibwa's CPU). Break-even: residual cost ≲1.2 ms/read keeps CPU ≤ minibwa; ≲2.4 ms stays 3× below BWA-MEM2. |
| **C3** | Kernel + scheduling: **widen SIMD to AVX2/AVX-512** (every kernel today is **128-bit only** — NEON/SSE, 8×int16 or 16×int8; the repo's own "Next" lists batched DP, 4 clusters per AVX-512 register, and Myers/Hyyrö bit-parallel verification); **GPU tier R + GPU pair pass**; port `--bio` band DP to GPU (its edit certificate is now 48% of `--bio` time; 11% of reads take 87%); the `d1=2` shortcut (explicitly not done). | largest untapped kernel lever; output-preserving (differential-testable) |
| **C4** | Fuse the pipeline: stream mapper → coordinate-bucketed sort → streaming markdup → CRAM/BAM (DRAGEN/Parabricks style). | end-to-end wall-clock (Amdahl) |

**Joint pair seeds** cost +45–47% CPU (+68% on repeat-enriched subsets);
`--iv-pair-votes 0` gives it back — gate per sample.

### Enabling work — evaluation rigour

- **Multiple correctness criteria:** 150 bp locus *and* strict 20 bp
  clip-aware endpoint *and* lenient interval-overlap >0.1; chrX/chrY PAR
  equivalent unless masked.
- **Projection circularity:** truth is projected via minimap2/dipcall and
  minibwa shares minimap2 chaining — restrict to the v5.0q benchmark BED and
  stratify by GIAB v3.6 SD / low-mappability strata.
- **Downstream:** DeepVariant + `hap.py --engine vcfeval --gender male` vs
  v5.0q smvar + CMRG. Meaningful = several-percent overall or concentrated in
  SD strata.
- **Proofs + audit:** mechanise new lemmas (phantom-copy bound, gap-graded
  horizon soundness) in Lean; extend `certa check` for new tags.
- **Adaptive modes:** report static / adaptive / calibrated-static / oracle
  with nonparametric CIs from ≥6 runs, plus input-order/sharding/controller-fuzz
  checks.

---

## 4. New human-genome data that pushes accuracy further

Confirmed by literature search 2026-11-21 (see Sources). These target exactly
the missing-copy class that dominates the residual.

1. **T2T-CHM13 v2.0** (`GCA_009914755.4 / GCF_009914755.1`, marbl/CHM13): gapless genome **plus
   a finished chrY from HG002/NA24385**, Cen/Sat v2.1 satellite annotation,
   resolved rDNA arrays on acrocentric short arms, updated SD annotation, a
   mapping analysis set on AWS. A direct **upgrade of the hs1 (v1.1)** CERTA
   uses today → better sinks, better lift, chrY handling. *Use in B1/B2.*
2. **HPRC Release 2 / HPRC2** (Lucas et al., bioRxiv **2026**,
   `10.64898/2026.07.21.739710`): **460 haplotypes / 232 individuals** (~5×
   Release 1), capturing **>99% of common variation in All of Us v8**; **~50%
   fewer** structurally unreliable regions and per-haplotype base errors vs
   HPRC1; deep-learning polishing (Google/UCSC) at **<1 error per 500 kb**;
   **thousands of T2T chromosomes**; **>100,000 variants per individual** not
   describable against existing references; formal pangenome coordinates,
   pantranscriptome/panepigenome. The **copy catalog** for B4/B5.
3. **Pangenome-aware DeepVariant** (Asri et al., bioRxiv **2025**,
   `10.1101/2025.06.05.657102`, Google+UCSC): **>20% fewer variant-calling
   errors**, **greatest in segmental duplications**; best short-read call set
   **99.65% precision / 99.1% recall** (AVITI + vg Giraffe, T2T-Q100-v1.1
   truth). Confirms the gain is concentrated where CERTA's residual sits.
   Forward-looking: HPRC graphs up to **1,000 haplotypes** + personalized
   pangenomes.
4. **Population resources:** All of Us v8, 1KGP ONT pangenome — useful as
   priors (π_∅, locus copy-number) and as a cross-sample sink diagnostic,
   with the leakage caveat.

**The honest ceiling.** Short reads cannot fully resolve SD ambiguity
per-read (the SDrecall caution). "A level not known so far" is therefore not
a better per-read score; it is the combination no published aligner offers:

> certify over GRCh38 ∪ CHM13/HPRC2 sinks → report the best GRCh38 locus
> **and** a certified statement that a better copy exists outside GRCh38
> (sink-copy tag) → abstain with a calibrated phantom-copy posterior →
> resolve the true copy by copy-level joint EM (B5) using the near-complete
> HPRC2 haplotypes as the copy catalog.

The prior-art check (companion limits doc) confirms this combination is
novel, but the claim must be **qualified**: per-read placement benchmarks
*do* exist — VerityMap (real long reads vs assembly truth, 2022), minibwa
(diploid-assembly-*simulated* reads mapped to source, 2026), Q100/GQC
`readbench` (base-level error rates, not placement). What is novel is the
specific design: **real short reads + diploid-assembly truth projected to
GRCh38 + a missing-copy/unplaceable stratum where abstention is correct**,
together with the locus-specific phantom-copy MAPQ term (Qtip 2017 only
*proposed* global unseen-source scaling; Locityper's null is locus-genotyping,
not paralog-shaped). So the **measurement and the calibrated posterior are
publishable**, not only the fix — positioned against Qtip/Locityper/VerityMap.

---

## 5. Theoretical limits

Formal, publication-grade statements (with prior-art status) are in
[`Theoretical limits of certified read placement.md`](Theoretical%20limits%20of%20certified%20read%20placement.md).
Summary: the seed-completeness loss ceiling (Thm 1/A), the MAPQ posterior
saturation proposition, the reference-content ceiling + missed/tie/reference
oracle decomposition, the confidence ceiling on the sample's own genome
(~3/M at 96.43%), the phantom-copy posterior bound, the short-read SD
identifiability limit, the adaptive-determinism invariant, and the
speed/certification break-even cost model. SETH-hardness (Backurs–Indyk 2015)
is the published worst-case backdrop.

---

## 6. Host execution & disk-safety protocol (galaxy / a2)

**The disks must not fill.** Every host run follows this protocol. The bench
scripts already refuse the wrong host; add a disk guard alongside.

**Hosts:** `galaxy` = ARM64 Ampere Altra, 128 cores (benchmarks use 64) +
2×A100 (one used); `a2` = 4×Xeon Gold 5218, 64 cores, 4 NUMA nodes, no GPU.
Use `numactl --interleave=all` on a2 (NUMA-balancing variance is 58–60 s vs
90–92 s without it); take medians of 3.

**Before any run:**
1. `df -h` on the working filesystem **and** `/tmp` (and the scratch mount).
   Require free space ≥ index size + output size + 20 GiB safety margin.
2. `du -sh` of the benchmark working dir and the index dir; record a baseline.
3. Index sizes (GRCh38, s=1): v3 **20.4 GiB**, `--key-bits 8` 17.6 GiB,
   `0` 14.9 GiB; hs1 **21.3 GiB**. `--t2t` holds **two** indexes (peak RSS
   42.6–46.1 GiB) — prefer **B2 sinks** (one index) to avoid this. `--bio`
   run table adds **1.45 GiB** RAM (not disk).
4. Never build two full indexes concurrently; build once, reuse, and verify
   with `certa check-index`.

**During the run:**
5. **Stream, don't spill:** use the concurrent pipe
   `certa … -u /dev/stdout | minibwa … /dev/stdin` (already in
   `bench/run_bench.sh` `certa-pipe`) so fallback reads never hit disk as
   intermediate FASTQ. AsyncWriter queues are RAM (2 GiB), not disk.
6. Write outputs to a designated scratch with a quota; prefer **CRAM** over
   BAM where the downstream tool accepts it (smaller). Keep intermediate BAMs
   only as long as needed for scoring.
7. Poll `df` during long runs (e.g. every batch); **abort** if free space
   drops below the safety margin. Log disk usage with the timings.

**After the run:**
8. Delete intermediate BAMs/FASTQs once the per-read score/truth TSVs are
   extracted. **Keep the small TSVs, not the BAMs**, in the results tree
   (`bench/results/…`); raw TSVs are what the prose tables cite.
9. `du -sh` again; confirm the working dir returned to ~baseline. Record peak
   disk in the run log.
10. Indexes and reference data are shared, read-only, mmap'd — never duplicate
    them per run.

**Go/no-go per step:** each implementation step (§7) is confirmed on **both**
hosts (CPU paths on a2; GPU paths on galaxy) with (a) the unit/oracle test
green, (b) byte-identity or the pre-registered accuracy rule met, (c) the
disk protocol satisfied (peak free space never below margin). A step is not
"done" until all three hold on both hosts.

---

## 7. Implementation sequencing (step-by-step, after merge)

Execute one step at a time; each gets a test + a galaxy/a2 confirmation +
the disk check. Do not start the next until the current is confirmed.

1. **A1** offline ceiling test (hours, no mapper code, no downloads). *Gate:*
   does the ceiling leave enough behind to justify B?
2. **A2** gap-graded MAPQ + relative horizon → **A3a** `--bio` floor.
3. **B1** v2.0 + lift-and-commit, **B3** hygiene, **A3b** phantom term +
   **A4** prior, **A5** calibrator. Freeze on HG002, validate on HG005.
4. **B2** one-pass sinks (per-locus validated); **C1** adaptive routing +
   **C3** AVX-512 widening / GPU tier R / `--bio` GPU port.
5. **C2** CERTA-X exact completion (tier F + window collapse); **C4** fused
   output.
6. **B5** copy-level EM on FixItFelix genes; **B4** HPRC2 graph second opinion
   only if sinks leave much behind.

Nine of the twelve placement-report experiments need **no downloads**; the
reference data already in hand (Ensembl GRCh38, hs1→upgrade to CHM13 v2.0,
HG002 Q100 v1.1, HG005 HPRC2 haplotypes, the existing BAMs) covers most of A
and B1/B3.

---

## 8. Corrections to earlier analysis (2026-11-21)

- **`--bio` is ~16–17× default mode, not ~6×.** The "~5–6×" figure is the
  ratio to *minibwa*, not to CERTA's default. This makes C3 more important.
- **122/M and 95.1/M are different runs** (HG002 vs HG005, different N) — see
  §1. `--bio` has never been run on HG005.
- **Conformal calibration** is a roadmap-level direction (`RESEARCH_ROADMAP.md
  §2.4`), not present in `certa/` code or the placement report.
- **`--t2t` is features-only today**; the veto is bench-only; lift-and-commit
  is new code (B1).
- The merged PRs #44/#45 (2026-11-21) are **WES-tertiary only** (couple carrier
  screening, gnomAD coverage/pext); nothing in `certa/`, `THEORY.md`,
  `FINDINGS.md`, the benchmarks, or `SWOT_ANALYSIS.md` changed, so the aligner
  findings stand.

---

## 9. Skills & tooling

No installed skill encodes bioinformatics knowledge; the domain reasoning is
the assistant's + the repo docs + literature. Useful bundled skills:
`agent-delegation` (parallel deep-dives), `review` (per-change correctness
review), `dataviz` (reliability/calibration diagrams, matched-N curves),
`new-app` (only if prototyping B5 standalone), `workflow-creator` (encode the
benchmark protocol). The `bailian-*` skills are not relevant. Tooling: C++17/
CUDA on galaxy/a2, Python for benchmark/calibration, `gh`/`git` (authenticated
as the owner), Lean 4 for proofs.

---

## Sources (literature search 2026-11-21)

- https://humanpangenome.org/hprc-data-release-2/
- https://www.barthel-lab.com/publications/hprc2-a-human-pangenome-reference-with-near-complete-coverag/
- https://www.biorxiv.org/content/10.64898/2026.07.21.739710v1.full-text
- https://github.com/marbl/CHM13
- https://ccb.jhu.edu/T2T.shtml
- https://www.elementbiosciences.com/blog/pangenome-aware-variant-calling-sets-a-new-standard
- https://www.biorxiv.org/content/10.1101/2025.06.05.657102v1.full-text
- https://github.com/wwliao/hprc_release2_graph_variant_benchmarking
