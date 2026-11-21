# Theoretical limits of certified short-read placement

**Status:** publication-grade statements for the PhD thesis. Prior-art
verdicts below were confirmed by a skeptical literature check (2026-11-21,
~30 searches + full-text reads of MAQ, Qtip, Locityper, minibwa, Q100,
BWA-MEM source `bwamem.c`, bwa-mem3 docs). Verdicts: **PUBLISHED** (cite, do
not claim), **PARTIALLY** (content known, formalisation/packaging new),
**LIKELY NOVEL** (no direct prior art found). Re-run a bioRxiv/arXiv sweep
before submission — late-2026 preprints cannot be ruled out.
**Companion:** [`Research guideline for certified biologically-correct placement.md`](Research%20guideline%20for%20certified%20biologically-correct%20placement.md),
[`../THEORY.md`](../THEORY.md) (proofs), [`../lean/`](../lean/) (mechanised lemmas).

This document states the limits that bound what *any* certified short-read
aligner can achieve, separates restatements of known results from new
formalisations and genuine contributions, and records the closest prior art
for each. Notation follows `THEORY.md`.

## Notation

- `G` — the reference (GRCh38 unless stated); `r` — a read of length `L`;
  `S` — the sample's true (diploid) genome; `π(u)` — the projection of a
  sample locus `u` onto `G` (may be absent).
- `d(r,ℓ)` — semi-global unit-cost edit distance; `score(r,ℓ)` — the affine
  objective (BWA-MEM scores: +1/−4, gap −6−k, clip −5); `pen(r,ℓ)` — the
  `--bio` integer-Phred penalty (lower is better).
- `m` — number of disjoint read *parts*; `q` — index k-mer length; `s` —
  sampling stride; parts have length ≥ `q+s−1` (Lemma I).
- `w(x) = 10^(−x/10)`; MAPQ in Phred.

---

## Limit 1 — Seed-completeness loss ceiling

**Statement (edit radius; Theorem 1 / Theorem A, `THEORY.md:130,181`).**
Let `S` be `m` disjoint parts, each ≥ `q+s−1`, all occurrences enumerated.
Every locus with `d ≤ m−t` has ≥ `t` parts occurring exactly there, on
diagonals within ±`d` of the start diagonal; certified radius `R = m−t`. With
one-edit neighbourhoods (Lemma N, Theorem A2), `R = 2m−1`. Sound for *any*
choice of disjoint parts, including parts chosen adaptively after hit counts.

**Statement (affine-score loss ceiling; Lemma L/S, `THEORY.md:229,270`).**
For disjoint parts ≥ 5, an alignment with fewer than `t` exact parts scores
≤ `U_t = L − 5·max(0, m−(t−1)) − 5·[clips]`; the full span attains the maximum
(Lemma S). With `t=1` the certifiable score loss is ≤ `5(m−1)`, and
`m ≤ ⌊L/(q+s−1)⌋`:

| q (s=1, L=150) | max parts m | edit radius R=m−1 / 2m−1 | score-loss ceiling ≈5(m−1) |
|---|---|---|---|
| 22 | 6 | 5 / 11 | ~25–30 |
| 15 | 10 | 9 / 19 | ~45–50 |

This is the **q=22 completeness barrier**; `PMAX=8`/`KMAX=5` in `core.h` cap
`m`. Lemma L is sharp (parts of length 2 fail) and Lemma S is sharp (parts of
length 4 give counterexamples), both mechanised in Lean.

**Prior art — PARTIALLY.** The pigeonhole lemma and its per-part error-budget
generalisation are the **search-scheme** formalism: Kucherov, Salikhov, Tsur,
*Theor. Comput. Sci.* 638:145–158 (2016), DOI 10.1016/j.tcs.2015.10.043
(CPM 2014, arXiv:1310.1440); the one-edit-neighbourhood parts are the
**01\*0 lossless seeds** of Vroland, Salson, Bini, Touzet, *J. Discrete
Algorithms* 37:3–16 (2016), DOI 10.1016/j.jda.2016.03.002; optimal part
orders/budgets in Kianfar, Pockrandt, Torkamandi, Luo, Reinert,
arXiv:1711.02035 (2017/2018), bioRxiv 10.1101/301085. Adaptive
frequency-driven seed selection is published: LAST (Kiełbasa et al., *Genome
Res.* 21:487–493, 2011, DOI 10.1101/gr.113985.110), GEM's Adaptive Seeds
Filter (Marco-Sola et al., *Nat. Methods* 9:1185–1188, 2012,
DOI 10.1038/nmeth.2221), and the **Optimal Seed Solver** (Xin et al.,
*Bioinformatics* 32:1632–1642, 2016, DOI 10.1093/bioinformatics/btv670), which
provably returns the least-frequent set of seeds per read from known index
frequencies; see also Xin, *Algorithms Mol. Biol.* 15:10 (2020),
DOI 10.1186/s13015-020-00172-3, and the unification in Depuydt et al.,
*OASIcs* 131 (2025), DOI 10.4230/OASIcs.Manzini.9, and Gottlieb et al.
(SeARChA), *NAR Genom. Bioinform.* 7:lqaf025 (2025),
DOI 10.1093/nargab/lqaf025. **Do not claim the lemma or adaptive
frequency-based seeding as novel.** What is not in print: the *certificate*
framing (a soundness proof that rarest-first adaptive part selection after
observing hit counts preserves losslessness) and the named "loss ceiling"
packaging — new packaging of a corollary, not a new mechanism.

---

## Limit 2 — MAPQ posterior saturation

**Proposition (saturation).** CERTA's `--bio` posterior is
`MAPQ = min(60, round(10·log10((1+Ω)/Ω)))`,
`Ω = w(YH−p₁) + Σ_{i≠pick} w(pen_i−p₁)` (`bio.h:993–999,1184–1198`), `YH` the
certificate horizon, `p₁` the best penalty. If the best locus is unique and
`YH−p₁ ≥ 60`, then `Ω ≤ 10⁻⁶` and **MAPQ = 60 regardless of `p₁`** (the read's
absolute divergence from `G`).

*Proof.* `w(YH−p₁) ≤ 10⁻⁶`; the other-loci sum is empty; `(1+Ω)/Ω ≥ 10⁶` ⇒
≥ 60 Phred, capped at 60. ∎

**Attribution (corrected).** This `(1+Ω)/Ω` posterior is **MAQ's** mapping
quality (Li, Ruan, Durbin, *Genome Res.* 18:1851–1858, 2008,
DOI 10.1101/gr.078212.108, Eq. 1: `p_s = p(z|x,u)/Σ_v p(z|x,v)`,
`Q_s = −10log10(1−p_s)`), with weights `10^(−gap/10)`. **It is not BWA-MEM's
formula.** BWA-MEM's actual MAPQ (verified in `lh3/bwa` `bwamem.c`,
`mem_approx_mapq_se`) is score-gap-based, `≈ 6.02·(score−sub)/a·tmp²` with
`tmp` carrying an **identity² damping**, minus `4.343·log(sub_n+1)`, capped at
60; with no suboptimal hit `sub` defaults to `min_seed_len·a = 19`. So
"regardless of absolute divergence" is *strictly* true of CERTA's gap-only
posterior and only *approximately* true of BWA-MEM (it still saturates at 60
for unique hits across a wide divergence range, but the identity² term damps
very divergent reads).

**Corollary (overconfidence in copy-deficit loci).** A read with many
high-quality mismatches (four Q40 mismatches, `p₁ ≈ 124–144`) still gets MAPQ
60 when its true competing copy is absent from `G`: the denominator is
conditioned on `G`'s competitor set. The certificate ("no better locus *in G*
below YH") is true; the posterior ("this locus is correct") is not. The model
prices a Q40 mismatch at ~31–36 Phred vs the measured ~2 Phred per competitor
mismatch — a ~9–20× overstatement.

**Prior art — PARTIALLY.** The formula is MAQ (2008). The overconfidence
mechanism is documented repeatedly: MAQ's Discussion already warns
sample–reference differences "can invalidate the calculation of mapping
qualities" (with an extra-copy *S. paratyphi* example); Qtip (Langmead,
*Genome Biol.* 18:152, **2017**, DOI 10.1186/s13059-017-1290-3) measures MAPQ
miscalibration and *assumes absent-source errors rare enough to ignore*;
Aganezov et al., *Science* 376:eabl3533 (2022), DOI 10.1126/science.abl3533
(false-PSV pileups from collapsed dups); Aqil & Gokcumen, bioRxiv
10.1101/2025.10.08.681258 (2025); **minibwa** (Li & Homer, arXiv:2606.15357,
2026) states verbatim that extra-copy reads "will be mismapped at high mapping
quality"; bwa-mem3 docs (`--supp-rep-hard-cap`) document competitors filtered
before scoring ⇒ MAPQ 60. **New packaging only:** a formal named statement
that horizon-limited unique hits saturate at the cap independent of absolute
divergence, with the unseen-competitor-at-horizon derivation. No verbatim
theorem found, but it is a corollary of published theory.

---

## Limit 3 — Reference-content ceiling and the oracle decomposition

**Proposition (reference-content ceiling).** Let
`ℓ* = argmax_{ℓ∈G} score(r,ℓ)` be any reference-only score-maximising
placement and the true source copy project to `π(u) ∈ G` (or be absent). If
`score(r,π(u)) < score(r,ℓ*)`, then **no** reference-only score-maximiser
places `r` at its truth. This class is irreducible for the whole family of
reference-only aligners, independent of search quality.

**Definition (oracle decomposition).** For reported `ℓ_rep` and truth `π(u)`,
re-align `r` at both under a fixed objective:
- **missed** — `score(r,π(u)) > score(r,ℓ_rep)`: search failure (avoidable);
- **tie** — equal: representation/ambiguity;
- **reference** — `score(r,π(u)) < score(r,ℓ_rep)`: irreducible for
  reference-only maximisers (truth scores *worse*).

**Measurement (repo, `bench/truth_oracle.py`).** minibwa test B: 382 missed /
9,435 tie / 19,155 reference; **88.0 of 88.7 confident errors/M are
reference-class** (test A: 425.1 of 429.1). 16,639 real mates (0.052%) placed
at truth by *no* GRCh38 run; BWA-MEM2 on CHM13 rescues 38.3%, Giraffe 27.9%.
For CERTA only ~9 of 122/M are avoidable by another GRCh38 run.

**Prior art — PARTIALLY.** The *ceiling* is a restatement of known
reference-bias facts — cite as a lemma, not a discovery: MAQ (2008); minibwa
(2026); Aganezov et al. (2022, 203 loci / 8.04 Mbp missing copies, 12
false-dup loci); Vollger et al., *Science* 376:**eabj6965** (2022),
DOI 10.1126/science.abj6965 (SDs absent from GRCh38; note abj6987 is Nurk et
al., the T2T sequence); levioSAM2 (Chen et al., *Nat. Methods*, 2024,
DOI 10.1038/s41592-023-02069-6, map-to-CHM13-then-lift, 11–52% fewer variant
errors); FixItFelix (Behera et al., *Genome Biol.* 24:31, **2023**,
DOI 10.1186/s13059-023-02863-7, masks ~1.2 Mbp false dups); SDrecall (Yang,
She et al., *Genome Biol.*, 2026, DOI 10.1186/s13059-025-03928-5). The closest
taxonomy is Qtip's MAQ-derived categories (1 = origin not in reference ≈
"reference"; 2 = origin in reference but no alignment ≈ "missed"; 3 = wrong
locus). **The per-read missed/tie/reference *truth-score* oracle decomposition
— scoring the true locus and comparing to the best found score — is LIKELY
NOVEL** (a score-based refinement of the outcome-based MAQ/Qtip taxonomy; no
direct prior art found). Monitor the unpublished `kira-ls-aligner` README
(ARyaskov) as conceptually adjacent concurrent work.

---

## Limit 4 — Confidence ceiling on the sample's own genome

**Statement (measured floor).** Placing the same reads on the sample's *own*
diploid assembly `S` (which contains the true copies) rather than `G` collapses
confident misplacement to **~3 per million at 96.43% of mates**
(`FINDINGS.md`, "how far from the limit"). This is the floor imposed by read
length and residual within-genome ambiguity once reference content is no
longer the bottleneck.

**Interpretation.** The gap between the GRCh38 confident-error rate
(~95–122/M) and this ~3/M floor is attributable to reference content (Limit 3)
— the headroom sinks/pangenome augmentation can recover. The ~3/M itself is
the short-read identifiability limit (Limit 6).

**Prior art — LIKELY NOVEL as an empirical measurement** (not a theorem). The
*principle* (map to the sample's own genome to remove reference-content
errors) is published since AlleleSeq (Rozowsky et al., *Mol. Syst. Biol.*
7:522, 2011, DOI 10.1038/msb.2011.54; vcf2diploid; ~0.5% of reads map to the
personal genome but not the reference), MAQ's iterative-reference suggestion,
WASP (van de Geijn et al., *Nat. Methods*, 2015), and most recently Vaddadi
et al., *Genome Res.* 36:740–753 (2026), DOI 10.1101/gr.280989.125 (imputed
personalized reference). VerityMap (2022) counts misplacements on complete
assemblies (long reads); Q100/GQC `readbench` reports base-error rates, not
placement. **No published equivalent of the high-MAPQ ~3/M confidence-ceiling
measurement was found.** Publishable as the headline empirical result with the
Limit-6 benchmark; cite AlleleSeq/Vaddadi for the principle.

---

## Limit 5 — Phantom-copy posterior bound

**Statement.** A calibrated placement posterior in a copy-deficit locus must
include the unseen-copy hypothesis:
`P(ℓ*|r) ∝ w(p₁) / [w(p₁) + Σ_{i≠*} w(pen_i) + w(YH−p₁) + π_∅(locus)·ℓ_phantom(r)]`,
`ℓ_phantom/ℓ(u*) ≈ Π_{mismatches} δ/(θ+ε) × (1−δ)^{matches}`,
`θ` the SNV rate, `ε` the sequencing error rate, `δ ≈ 10⁻²` for recent SDs,
`π_∅(locus)` the prior probability of an unseen copy (CHM13 copy excess, HPRC
duplications, GIAB false-dup lists, or sample depth/k-mer excess).

**Bound.** Dropping the phantom term (Limit 2) overstates MAPQ up to the cap;
adding it reduces MAPQ in copy-deficit loci by an amount governed by `δ`,
`π_∅`. Fitting `δ, π_∅` to the measured ~2 Phred/mismatch slope calibrates the
posterior. This is the BWA-PSSM structure with a *paralog-shaped* (not
i.i.d.-random) alternative.

**Prior art — PARTIALLY; the specific term is LIKELY NOVEL (strongest theory
candidate).** The general idea is published **as a suggestion**: Qtip (2017)
states verbatim that "category 1 errors could be included in our model, e.g.
by assuming a global prevalence of category 1 errors and scaling predictions
accordingly" — global, not locus-specific, not paralog-shaped. Locityper
(Prodanov et al., *Nat. Genet.*, 2025, DOI 10.1038/s41588-025-02362-4;
bioRxiv 2024.05.03.592358) has a penalised null window
`p(r,w○) = Λ·max_h p(r,h)`, `Λ = 10⁻⁵`, plus an "unexplained reads" filter —
but it is a *per-locus genotyping* likelihood, the null is an escape sink for
unalignable reads (not a paralog-shaped copy at divergence δ), and it cannot
predict haplotypes missing from its database. BWA-PSSM is **Kerpedjiev,
Frellsen, Lindgreen, Krogh**, *BMC Bioinformatics* 15:100 (2014),
DOI 10.1186/1471-2105-15-100 (position-specific error/damage priors, not
missing paralogs — *not* Li). Paralog-aware tools recover copies present in
the reference or by assembly, with no probabilistic unseen-copy term:
DuploMap (Prodanov & Bansal, *Nucleic Acids Res.* 48:e114, 2020,
DOI 10.1093/nar/gkaa829), Parascopy (*Nat. Commun.*, 2022,
DOI 10.1038/s41467-022-30930-3), Paraphase (*Nat. Commun.* 16:2340, 2025,
DOI 10.1038/s41467-025-57505-2). **A locus-specific phantom-copy term in a
genome-wide placement posterior, with a paralog-divergence-parameterised
likelihood replacing the i.i.d. null, was not found.** It *must* be positioned
against Qtip's explicit proposal and Locityper's null window, or a reviewer
will call it incremental.

---

## Limit 6 — Short-read segmental-duplication identifiability

**Statement.** A read of length `L` uniquely identifies its source copy among
`k` paralogs only if it spans a *copy-diagnostic* site (where the copies
differ). If the mean spacing `d` of diagnostic sites exceeds `L`, a fraction
of reads span no diagnostic site and are **fundamentally ambiguous per-read**,
regardless of alignment quality. Per-read posteriors cannot resolve this;
aggregation across reads (copy-level joint EM, guideline B5) can, because
mismatches recurring at the same offsets across reads vote for an unseen copy.

**Interpretation.** The information-theoretic reason "better than any aligner"
cannot mean a better per-read score in SDs; it means certified abstention +
cross-read copy resolution. The SDrecall authors' caution is the empirical
statement of this limit.

**Prior art — PARTIALLY (formalisation new).** SDrecall (2026,
DOI 10.1186/s13059-025-03928-5) and SD genotyping (Locityper, Parascopy,
Paraphase) embody the difficulty empirically; no formal identifiability
statement was found. State it as a lemma motivated by these, not a discovery.

---

## Limit 7 — Adaptive-determinism invariant (engineering theorem)

**Statement.** Output `O` is byte-identical across thread counts, batch sizes,
CPU/GPU back-ends, input order and sharding **iff** every per-read
search-effort decision is a function of `(read, index, frozen parameters)`
alone, not of scheduling state. *Execution-only* levers (scheduling,
cache-locality reordering, memoisation that does not change which alignments
are evaluated) preserve `O` by construction; *search-effort* levers must be
frozen per read after calibration; *certified-problem* parameters (scoring,
the fixed 850-bp pair distance) must never adapt online.

**Prior art — PARTIALLY.** The problem is documented: BWA `-K`; bwa-mem3
threading docs ("`-t` changes the output … each batch's `mem_pestat`
insert-size estimate is drawn from the reads in that batch … pass `-K INT`");
bwa-mem2/bwa-mem3 equivalence docs; and Icer Baykal et al., bioRxiv
10.1101/2025.05.08.652934 (2025), measuring up to 13.53% inconsistent mappings
among 8 aligners under read shuffling/reverse-complementing (input-perturbation
reproducibility, not thread/batch/hardware determinism). **A formal
cross-platform invariant via frozen per-read effort rules: no direct prior art
found** — new packaging of an engineering property; system-paper material.

---

## Limit 8 — Speed / certification break-even cost model

**Statement.** With residual fraction `f`, exact-completion cost `c` per
residual read, certified-tier CPU `C_cert`, read count `N`, baseline CPU
`C_base`: total `= C_cert + f·N·c ≤ C_base` iff `c ≤ (C_base − C_cert)/(f·N)`.
With `f ≈ 0.046`: **`c ≲ 1.2 ms/read`** to stay ≤ minibwa CPU, **`c ≲ 2.4
ms/read`** to stay 3× below BWA-MEM2 (`EXACT_ALIGNMENT_PROPOSAL.md`).

**Prior art — restatement** of Amdahl-style cost accounting; a go/no-go budget
for CERTA-X (guideline C2), not a novel theoretical claim.

---

## Limit 9 — SETH-hardness of edit distance (published backdrop)

**Theorem (Backurs–Indyk, STOC 2015, DOI 10.1145/2746539.2746612).** Unless
SETH is false, edit distance cannot be computed in `O(n^{2−δ})` for any
constant `δ>0`. Conditional, worst-case, exact. Provable guarantees here are
necessarily *parameterised* (WFA `O(ns+s²)`, Ukkonen `O(nd)`, r-index `O(r)`),
*lossless-filter* (Limit 1), or *approximation* guarantees. Already in
`SWOT_ANALYSIS.md §J`. **PUBLISHED; restated.**

---

## Limit 10 — Window-collapse and count-based tie proofs

**Theorem R′ (window collapse, `THEORY.md:989`).** Identical reference windows
have identical band maxima; one DP per *distinct* window string suffices,
reducing the enumeration budget up to ~100× in near-identical repeats at the
same DP cost. (Not yet built; proposal §6.2.)

**Theorem R (ties, `THEORY.md:975`).** If the reference string under the best
alignment `A*` occurs at another non-overlapping position (or the other
strand), the same transcript gives an equal-score alignment elsewhere — a tie
proven by occurrence count, without a second DP.

**Prior art.** Window-collapse memoisation: **LIKELY NOVEL** as a specific
optimisation but engineering-only (DP memoisation is textbook; GenCache,
MICRO 2019, DOI 10.1145/3352460.3358308, is in-cache acceleration, not reuse
across identical windows; FPGA/ASIC caching literature not exhaustively
checked). Count-based tie proofs: **PUBLISHED in essence** as genome
mappability — GEM (Derrien et al., *PLoS ONE* 7:e30377, 2012,
DOI 10.1371/journal.pone.0030377) and GenMap (Pockrandt et al.,
*Bioinformatics* 36:3687–3692, 2020, DOI 10.1093/bioinformatics/btaa222)
compute exact (k,e)-frequency counts per position (count-based uniqueness
certificates via optimum search schemes); the per-read certificate framing
differs in presentation, not concept.

---

## Novelty summary (verified)

| Limit | Kind | Verdict |
|---|---|---|
| 1 seed-completeness loss ceiling | restatement + certificate packaging | PARTIALLY (cite Kucherov/Vroland/Kianfar/Xin) |
| 2 MAPQ posterior saturation | formalisation of a known mechanism | PARTIALLY (formula = MAQ 2008; mechanism known) |
| 3a reference-content ceiling | restatement | PARTIALLY (cite as lemma) |
| 3b missed/tie/reference truth-score oracle | new diagnostic | **LIKELY NOVEL** |
| 4 own-genome confidence ceiling (~3/M) | new measurement | **LIKELY NOVEL** (empirical) |
| 5 locus-specific phantom-copy posterior | new model term | **LIKELY NOVEL** (strongest theory; position vs Qtip/Locityper) |
| 6 SD identifiability | formalisation | PARTIALLY |
| 7 adaptive-determinism invariant | engineering theorem | PARTIALLY (system paper) |
| 8 break-even cost model | restatement | restatement |
| 9 SETH-hardness | published | PUBLISHED |
| 10a window-collapse memoisation | engineering | LIKELY NOVEL (weak theory) |
| 10b count-based tie proofs | restatement | PUBLISHED in essence (GEM/GenMap) |

**Publishable contributions (with positioning):** the locus-specific
phantom-copy posterior (Limit 5, against Qtip's category-1 scaling proposal
and Locityper's null), conformal MAPQ thresholding (guideline A5; application
novelty — address exchangeability honestly, likely needs Mondrian/stratified
CP keyed on mappability/SD strata, cf. Papangelou et al., *Front. Bioinform.*
5:1507448, 2025, DOI 10.3389/fbinf.2025.1507448; ProtCOVER, bioRxiv
10.64898/2026.11.13.766406, shows marginal guarantees fail under group shift),
the truth-score oracle decomposition (Limit 3b), and the real-read
GRCh38-projected diploid-truth placement benchmark with a missing-copy
abstention stratum + the ~3/M confidence ceiling (Limits 4/6-benchmark).

**The benchmark claim must be qualified.** The placement report's unqualified
"no per-read placement benchmark exists" is **false as stated**: VerityMap
(Bzikadze, Mikheenko, Pevzner, *Genome Res.* 32:2107–2118, 2022,
DOI 10.1101/gr.276871.122) benchmarks real long-read placement against complete
assembly truth (wrong-haplotype counts); minibwa (2026) benchmarks
diploid-assembly-*simulated* reads mapped to source (IoU ≥ 10%); Q100/GQC
`readbench` (Hansen et al., bioRxiv 10.1101/2025.09.21.677443, 2025; Cell 2026)
reports base-level error rates, not placement; GIAB v5.0q (Olson et al.,
bioRxiv 10.64898/2026.09.23.752440, 2026) and synDip (Li et al., *Nat. Methods*
15:595–597, 2018, DOI 10.1038/s41592-018-0054-7) are variant-centric. **What is
novel:** *real short reads*, truth = *diploid assembly projected to GRCh38*
(dipcall-style), scored *per-read for placement*, with an explicit
*missing-copy/unplaceable stratum where abstention is the correct answer*.

---

## Mandatory citation corrections (apply before the thesis cites these)

- Posterior `(1+Ω)/Ω` formula → **MAQ (Li, Ruan, Durbin 2008)**, *not* BWA-MEM.
  BWA-MEM MAPQ is score-gap-based with identity² damping (`bwamem.c`).
- "Kucherov et al. 2018, A unifying framework" → **2006**, *J. Bioinform.
  Comput. Biol.* 4:553–569, DOI 10.1142/S0219720006001977. Search schemes =
  CPM 2014 / TCS 2016.
- "Optimal search schemes (Rahn et al.)" → first author **Kianfar**;
  arXiv:1711.02035 / bioRxiv 10.1101/301085 (journal version unverified).
- "Qtip (2024)" → **2017**, *Genome Biol.* 18:152, DOI 10.1186/s13059-017-1290-3.
- "BWA-PSSM (Li 2014)" → **Kerpedjiev et al. 2014**, *BMC Bioinformatics*
  15:100, DOI 10.1186/1471-2105-15-100.
- "FixItFelix (2024)" → **2023**, *Genome Biol.* 24:31,
  DOI 10.1186/s13059-023-02863-7.
- Vollger DOI → 10.1126/science.**abj6965** (abj6987 = Nurk et al. T2T).
- Sirén/Giraffe DOI → 10.1126/science.**abg8871**.
- SDrecall → *Genome Biol.* 2026, DOI 10.1186/s13059-025-03928-5 (author list
  partial: Yang X, She CH, et al.).
- Q100 → cite bioRxiv 10.1101/2025.09.21.677443 (Cell 2026 DOI unverified).
- Unverified, do not cite without checking: Ruffalo 2012 (*Bioinformatics*
  28:i349) DOI; the 2020 MAPQ-recalibration paper (PMC7727374) venue.

**Gray literature used (not peer-reviewed):** bwa/bwa-mem2/bwa-mem3 source +
docs (2026); qcfail.com MAPQ blog (2016); `kira-ls-aligner` README (ARyaskov,
unpublished — monitor as concurrent work).
