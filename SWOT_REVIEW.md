# SWOT_ANALYSIS.md — correctness review (2026-11-21)

**Scope:** full read of all 1864 lines of [`SWOT_ANALYSIS.md`](SWOT_ANALYSIS.md)
plus ~50 targeted web verifications (publisher / PubMed / arXiv / LIPIcs pages
and the Ukkonen 1985 PDF). Prepared for the PhD thesis literature review.
**Verdict: the document is unusually accurate.** Nearly every complexity
claim, attribution and DOI spot-checked is correct, including several commonly
gotten wrong (GenASM ≠ Myers; Backurs–Indyk "*Strongly* Subquadratic"; r-index
JACM Article 2; strobealign/syncmers sole-authorship; WFA-adaptive is a
heuristic). The defects are a handful of small factual slips, a few
overstatements, and one genuine completeness gap. The corrections below have
been **applied in place** in this PR; this file records what changed and why.

Legend: **[A]** definite error (fixed) · **[B]** overstatement/nuance
(softened) · **[C]** completeness gap (entries added) · **[D]** verified
correct (left unchanged).

---

## [A] Definite factual errors — corrected

- **A1. K.3 strobealign-MCS author initial.** Third author is **Nicolas
  Buchin** (Stockholm Univ.), not "K. Buchin". Fixed `Buchin, K.` →
  `Buchin, N.`; added locator *Genome Biology* **27**, 118 (2026). Verified
  author list: Tolstoganov I, Martin M, Buchin N, Sahlin K (PubMed 41764549).
  DOI `10.1186/s13059-026-04017-x` was already correct.
- **A2. K.4 ONT `lr:hq` misdated.** The preset was added in **minimap2 v2.27
  (r1193), 12 March 2024**, not "since 2025" (minimap2 NEWS.md/releases).
  Reworded to "added in minimap2 v2.27, 2024".
- **A3. K.3 mapquik speedup.** The cited *Genome Research* 2023 version
  reports **37×** on human (410× maize); 30×/350× are the bioRxiv preprint
  numbers. Updated to ~37× (noting 30× preprint). Citation
  `10.1101/gr.277679.123`, 33(7):1188–1197 otherwise correct.

## [B] Overstatements / missing nuance — softened

- **B1. B.6 Ukkonen "provably optimal".** O(nd) is an *upper bound*; no
  matching Ω(nd) lower bound is known. Replaced with "output-sensitive /
  parameterized". Added the score-only refinement: Ukkonen's Theorem 2 gives
  O(s·min(m,n)) time/space with the edit sequence, reducible to
  O(min(s,m,n)) score-only (verified in the 1985 PDF, p.107).
- **B2. K.1 A\*PA2 "on all tested datasets".** The paper says "often faster
  and never much slower", with an explicit exception: on ~10–11 kbp ONT reads
  Block Aligner is ~1.6–2× faster (but reports only ~53% of alignments
  optimally); comparators are WFA-adaptive and Block Aligner, not minimap2.
  Reworded accordingly.
- **B3. F.19 SA-IS "first … O(n)".** Linear-time SA construction predates
  SA-IS (Kärkkäinen–Sanders 2003; Ko–Aluru 2003/2005). Softened to "the first
  simple, compact (<100 line) almost-pure induced-sorting linear-time
  construction".
- **B4. F.17 / inventory #17 BWT build bound.** With SA-IS the BWT is
  buildable in **O(n)**; O(n log n) is an older bound. Updated inventory and
  body to "O(n) via SA-IS; naive rotation sort O(n² log n)".
- **B5. K.3 minibwa "comparable accuracy".** Added the centromere caveat:
  marginally lower accuracy than BWA-MEM on simulated WGS short reads, driven
  by centromeric/acrocentric regions, no variant-calling impact; noted the
  repo README's ">3×/2×" vs the paper's ~4×.
- **B6. K.3 Winnowmap2.** The description covered only Winnowmap (2020).
  Added Winnowmap2's key addition: **MCAS** (minimal confidently alignable
  substrings) + paralog-specific-variant sensitivity.
- **B7. K.4 DRAGEN coverage.** The paper benchmarks **35×** 1kGP genomes at
  ~30 min. Changed "30x" → "~35x".
- **B8. Inventory #62 year range.** Cited entries span 2021–2024 (Movi =
  iScience 2024). Changed "2021–2025" → "2021–2024".

## [C] Completeness gap — added K.7 "Reference content, coordinates and calling"

The SWOT covered repetitive/low-complexity regions and clinically important
genes but had **no entry** for T2T/HPRC reference content,
false-duplication/collapse errors, lift-over, pangenome-aware or SD-specific
calling, or MAPQ calibration — precisely the class that bounds CERTA's
residual accuracy (see [`certa/reports/Theoretical limits of certified read placement.md`](certa/reports/Theoretical%20limits%20of%20certified%20read%20placement.md),
Limits 3–6). Added subsection **K.7** + inventory rows + a cross-cutting
bullet ("reference-content errors, not just algorithmics, bound achievable
accuracy"). Entries added (all web-verified):

- **T2T-CHM13 / v2.0** — Nurk et al. 2022, *Science* 376:44–53,
  `10.1126/science.abj6987`; v2.0 = **GCA_009914755.4 / GCF_009914755.1**,
  adds complete chrY from HG002 (T2T-Y): Rhie et al. 2023, *Nature*
  **621**:344–354, `10.1038/s41586-023-06457-y`.
- **HPRC draft pangenome** — Liao et al. 2023, *Nature* 617:312–324,
  `10.1038/s41586-023-05896-x` (47 individuals / 94 haplotypes).
- **HPRC2** — Lucas et al. 2026, bioRxiv `10.64898/2026.07.21.739710`
  (460 haplotypes / 232 individuals; >99% of common variation, All of Us v8;
  halves unreliable regions and base errors vs HPRC1).
- **FixItFelix** — Behera et al. **2023**, *Genome Biol.* 24:31,
  `10.1186/s13059-023-02863-7` (~1.2 Mbp false dups + ~8.04 Mbp collapsed SDs,
  33 protein-coding genes; CMRG SNV F 0.012→0.980 in false-dup genes). *Not 2024.*
- **levioSAM2** — Chen et al. 2024, *Nat. Methods* 21:41–49,
  `10.1038/s41592-023-02069-6` (map to CHM13 then lift; commit/defer/suppress).
- **Pangenome-aware DeepVariant** — Asri et al. 2025, bioRxiv
  `10.1101/2025.06.05.657102` (up to ~25.5% fewer variant errors, greatest in
  SDs; **caveat:** induces ~half as many new errors as it fixes, 5,497 induced
  vs 9,272 fixed SNPs, concentrated in high-identity SDs). DeepVariant itself:
  Poplin et al. 2018, *Nat. Biotechnol.* 36:983–987, `10.1038/nbt.4235`.
- **SDrecall** — Yang, She et al. 2026, *Genome Biol.* 27:23,
  `10.1186/s13059-025-03928-5` (~95% SD short-variant sensitivity from short
  reads alone).
- **Qtip (MAPQ calibration)** — Langmead **2017**, *Genome Biol.* 18:152,
  `10.1186/s13059-017-1290-3`. **There is no 2024 Qtip paper** — the 2024
  "QTIP" (arXiv:2406.11235) is LLM weight quantization; do not cite it.
- **Reference-content impact / truth** — Aganezov et al. 2022, *Science*
  376:eabl3533, `10.1126/science.abl3533`; GIAB Zook et al. 2020,
  *Nat. Biotechnol.* 38:1347–1355, `10.1038/s41587-020-0538-8`.

**Not added (unverified DOIs — confirm before citing):** double-decycling
(Pellow et al. 2023, PMC10538483) and GreedyMini (Golan et al. 2025,
*Bioinformatics* 41(Suppl 1):i275) — DOIs not individually resolved. FILTR now
has a published version (ACM `10.1145/3839504`) alongside arXiv:2607.06225 —
optional update.

## [D] Verified correct — left unchanged (representative)

- **Complexity:** NW/Gotoh/SW O(mn); Hirschberg O(mn) time / O(min(m,n)) space;
  Myers O(n·⌈m/w⌉); WFA O(ns+s²), BiWFA O(s) space; FM-index O(m) count +
  O(occ·s_SA) locate; SA-IS O(n); r-index O(r) space, JACM 67(1) Article 2,
  `10.1145/3375890`; minimap2 bounded-lookback h≈50 + `--rmq`; Backurs–Indyk
  conditional/worst-case/exact; Navarro O(|E|·m); Jain NP-hard only with graph
  edits; Equi OV/SETH lower bound.
- **Attributions:** GenASM = modified Bitap/Wu-Manber, NOT Myers;
  Suzuki–Kasahara = KSW2 core; WFA btaa777 / BiWFA btad074; A\*PA btae032 /
  A\*PA2 WABI 2024 LIPIcs 312:17 (sole author); minibwa arXiv:2606.15357;
  DRAGEN `10.1038/s41587-024-02382-1`; strobealign 2022 sole-author Sahlin;
  syncmers sole-author Edgar; BWA-MEME btac137 (bit-identical output);
  giraffe Science abg8871; HISAT2 s41587-019-0201-4.
- **WFA-adaptive IS a heuristic** (drops lagging diagonals; WFA2-lib default
  `wf-adaptive` is approximate). Entry 9's "fast variants give up the exactness
  guarantee" is correct — do not "fix" it.
- **Internal consistency:** inventory ↔ body ↔ code-sample table match; every
  linked repo file exists; TOC anchors valid; no dangling cross-references;
  "over 60" tools = 64 rows. The Sept-2026 "all DOIs checked" audit claim holds
  for every DOI resolved here.

## Unverified / flagged (consistent with the canonical record; not re-fetched)

- DOIs not individually re-fetched (high confidence): SneakySnake btaa1015;
  GASAL2 s12859-019-3086-9; MaxSSmap 1471-2164-15-969; SOAP3-dp pone.0065632;
  BarraCUDA 1756-0500-5-27; Parasail s12859-016-0930-z; SSW pone.0082138;
  Farrar btl582; BLAST; Roberts bth408; BWT TR-124; FM-index
  10.1109/SFCS.2000.892127; SA-IS 10.1109/DCC.2009.42; Hirschberg
  10.1145/360825.360861; NW/SW/Gotoh JMB DOIs; Myers 1999/1986;
  Landau–Vishkin; Ukkonen DOI; Bowtie/Bowtie2; Mash/MashMap; Block Aligner
  btad487; mm2-fast s43588-022-00201-8; Sentieon 115717; Kianfar
  arXiv:1711.02035; Liao 2023; Zook 2020.
- TALCO DOI `10.1109/HPCA57654.2024.00044`: paper/venue/authors confirmed
  (IEEE doc 10476438) but the exact DOI string not directly resolved.
- Author-reported speedups not independently re-measured: GenASM "22–12,501×"
  (116× confirmed); Parasail "136 GCUPS"; BWA-MEM2 "1.3–3.1×"; mm2-fast
  "1.6–1.8×"; Bowtie "~1.3 GB". MaxSSmap "no maintained public repository"
  not verified.

**Bottom line:** A1–A3 fixed, B1–B8 softened, K.7 reference-content /
pangenome-calling / MAPQ entries added (with corrected FixItFelix=2023,
Qtip=2017, T2T accession GCA_…4/GCF_…1, Rhie=Nature 621). Everything in [D] is
sound and unchanged.
