# MAPQ models and calibration in short-read aligners, and features that flag reads from paralogs or copies missing from the reference

Scope: human Illumina short reads, 2010–2026. Framing from our own measurements (HG002 2x150 NovaSeq X, GRCh38 primary without decoys or ALTs, Q100 truth, locus-level scoring at matched sensitivity): minibwa misplaces 89 per million confident mates, strobealign 111, BWA-MEM2 114, CERTA 122, minimap2 130. Nearly all errors are on mates whose true copy is missing from GRCh38 or diverged there, so what separates aligners is MAPQ ranking. Unless marked "Inference", every statement below is sourced.

## Q1. How each major aligner computes MAPQ (source-level)

### Takeaway
Every aligner examined (BWA-MEM/BWA-MEM2, minimap2 -x sr, strobealign, DRAGMAP/DRAGEN, LAST, Stampy, vg Giraffe) computes MAPQ mainly from the score gap between the best hit and the best competing hit it found, with smaller terms for the count of near-best hits, seed repetitiveness, and pairing. When no competitor is found, they all go to about 60, whatever the absolute divergence of the best hit. Bowtie2 end-to-end is the one exception that matters here: its table scales MAPQ by the best hit's absolute score relative to the minimum valid score, so a unique hit with many high-quality mismatches gets a low MAPQ.

### Cited Findings

**BWA-MEM (bwamem.c, `mem_approx_mapq_se`)**
- Verbatim logic: `sub = a->sub ? a->sub : opt->min_seed_len * opt->a; sub = a->csub > sub ? a->csub : sub; if (sub >= a->score) return 0;`. Then `l = max(qe-qb, re-rb)` and `identity = 1 - (l*a - score)/(a+b)/l` — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- The default branch runs when `mapQ_coef_len > 0`: `tmp = l < mapQ_coef_len ? 1 : mapQ_coef_fac/log(l); tmp *= identity*identity; mapq = (int)(6.02 * (score - sub)/a * tmp*tmp + .499)` — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- The defaults are `mapQ_coef_len = 50; mapQ_coef_fac = log(mapQ_coef_len)` (in `mem_opt_init`), so the default path is the 6.02·Δscore formula with a length factor and an identity⁴ factor — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- The legacy branch runs only when `mapQ_coef_len <= 0`: `mapq = MEM_MAPQ_COEF * (1 - sub/score) * log(seedcov)`, then `if identity < 0.95: mapq *= identity²`. Here `MEM_MAPQ_COEF = 30.0` and `MEM_MAPQ_MAX = 60`. So the "identity < 0.95" penalty and seed coverage are **not** in the default path — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c); [bwamem.h](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.h)
- Post-terms in all branches: `if (sub_n > 0) mapq -= (int)(4.343*log(sub_n+1) + .499)`, then a clamp to [0,60], then `mapq = (int)(mapq * (1 - frac_rep) + .499)` — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- `frac_rep` is the fraction of read bases covered by seeds with more than `max_occ` occurrences (`chain.a[i].frac_rep = (float)l_rep / len`). It is the k-mer-multiplicity term: a read whose seeds are mostly repetitive has its MAPQ scaled down even with a good score — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- `sub` is the score of the first lower hit overlapping a kept hit by ≥ `mask_level`. `sub_n` counts hits with `score_best − score_i ≤ max(a+b, o_del+e_del, o_ins+e_ins)` — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- Struct comments: `sub` is the "2nd best SW score", `csub` the "SW score of a tandem hit", `sub_n` the "approximate number of suboptimal hits", `seedcov` the "length of regions coverged [sic] by seeds" — [bwamem.h](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.h)
- Pair MAPQ (bwamem_pair.c): `raw_mapq(diff,a) = (int)(6.02*diff/a + .499)`. `subo = max(subo, score_un)` where `score_un = score_mate1 + score_mate2 − pen_unpaired`. `q_pe = raw_mapq(o − subo)`, then `−4.343·ln(n_sub+1)`, a clamp to [0,60], then `q_pe *= 1 − 0.5·(frac_rep1 + frac_rep2)` — [bwamem_pair.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem_pair.c)
- Per-mate combination: `q_se = q_se > q_pe ? q_se : (q_pe < q_se + 40 ? q_pe : q_se + 40)`, so pairing can raise a mate's MAPQ by up to 40. It is then capped by `raw_mapq(score − csub)` (the tandem-hit cap) — [bwamem_pair.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem_pair.c)
- Sub-optimal pairs: `n_sub` counts candidate pairs whose score is within `tmp = max(a+b, o_del+e_del, o_ins+e_ins)` of the runner-up pair — [bwamem_pair.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem_pair.c)

**BWA-MEM2**
- Its README says BWA-MEM2 "produces alignment identical to bwa", and that output should match bwa-mem 0.7.17. It does not mention MAPQ explicitly — [bwa-mem2 README](https://raw.githubusercontent.com/bwa-mem2/bwa-mem2/master/README.md)

**Bowtie2 (unique.h, `BowtieMapq2::mapq`)**
- Variables: `scPer` is the perfect score (plus the mate's if paired), `scMin` the minimum valid score, `diff = max(1, scPer − scMin)`, `bestOver = best − scMin`, and `bestdiff = |best − secbest|` — [unique.h](https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h)
- End-to-end with no second-best: `bestOver ≥ 0.8·diff → 42; ≥0.7 → 40; ≥0.6 → 24; ≥0.5 → 23; ≥0.4 → 8; ≥0.3 → 3` (else 0) — [unique.h](https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h)
- End-to-end with a second-best, by `bestdiff/diff` band (perfect best / imperfect best): ≥0.9: 39/33; ≥0.8: 38/27; ≥0.7: 37/26; ≥0.6: 36/22; ≥0.5: 35/25 (bestOver ≥ 0.84·diff) /16 (≥0.68) /5; ≥0.4: 34/21/14/4; ≥0.3: 32/18 (≥0.88)/15 (≥0.67)/3; ≥0.2: 31/17/11/0; ≥0.1: 30/12/7/0; >0: 6 (bestOver ≥ 0.67·diff)/2; =0: 1/0. Local mode uses a different table (44, 42, 41, 36, 28, 24, 22 with no second-best) — [unique.h](https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h)

**minimap2 (hit.c, `mm_set_mapq2`; the sr preset sets `is_sr`)**
- `q_coef = 40`. `uniq_ratio = sum_sc/(sum_sc + rep_len)` (rep_len is the repetitive-seed length). `pen_s1 = (score > 100 ? 1 : 0.01·score) · uniq_ratio`. `pen_cm = cnt > 10 ? 1 : 0.1·cnt`, then `pen_cm = min(pen_s1, pen_cm)`. `subsc = max(r->subsc, min_chain_sc)` — [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c)
- With base alignment, non-splice: `x = dp_max2·subsc/(dp_max·score0)`, then `mapq = identity · pen_cm · 40 · (1 − x²) · ln(dp_max/match_sc)`. `identity = mlen/blen`. The BWA-style `mapq_alt = 6.02·identity²·(dp_max − dp_max2)/match_sc` applies only when `!is_sr`, so not to short reads — [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c)
- Final terms: `mapq −= 4.343·ln(n_sub+1)`, then a cap at 60. `n_sub` is incremented for a secondary hit if `ri->cnt ≥ rp->cnt` or `dp_max(primary) − dp_max(secondary) ≤ sub_diff` — [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c)

**strobealign**
- v0.13.0 C++, single-end: `mapq = (60·(best_score − second_best_score) + best_score − 1)/best_score`, i.e. ⌈60·(1 − s2/s1)⌉ — [strobealign v0.13.0 aln.cpp](https://raw.githubusercontent.com/ksahlin/strobealign/v0.13.0/src/aln.cpp)
- Paired, `joint_mapq_from_high_scores`: 60 if there is only one candidate pair; 0 if `score1 == score2`; otherwise `min(60, score1 − score2)` when both scores are positive (a raw SW-score difference, not log-scaled) — [aln.cpp](https://raw.githubusercontent.com/ksahlin/strobealign/v0.13.0/src/aln.cpp)
- `proper_pair_mapq` (NAM-level): copied from minimap2 v1 as `40·(1 − s2/s1)·min(1, n_hits/10)·ln(s1)`, capped at 60. Secondary alignments get MAPQ 0 — [aln.cpp](https://raw.githubusercontent.com/ksahlin/strobealign/v0.13.0/src/aln.cpp)
- Current strobealign has been rewritten in Rust (`src/aligner.rs`, `chainer.rs`, `modes/` ...). Whether the MAPQ formulas carried over unchanged was not verified — [strobealign src listing](https://api.github.com/repos/ksahlin/strobealign/contents/src)

**DRAGEN / DRAGMAP**
- DRAGEN docs: "The primary contributor to MAPQ calculation is the difference between the best and second-best alignment scores." AS holds the best score and XS the second-best. When nothing competes, `aln-min-score` serves as the suboptimal score, so raising it lowers MAPQ for low-scoring alignments. `pe-max-penalty` caps how much a mate can raise MAPQ — [DRAGEN v4.2 DNA Aligning](https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/GPipelineAlign_fDG.htm) (via search summary); [DRAGEN v4.2 MAPQ Scoring](https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/TPipelineMAPQ_fDG.htm)
- For RNA multimappers, DRAGEN uses MAPQ = int(−10·log10(1 − 1/NH)) — [DRAGEN v4.2 MAPQ Scoring](https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/TPipelineMAPQ_fDG.htm)
- DRAGMAP `PairBuilder::updateEndMapq`: suboptimal-count penalty `sub_mapq_pen_v = sub_count ? 3·sub_count_log2 : 0` (fixed-point, `>>7`). XS-heuristic cap `xs_heur_mapq = max(0, ((xs_score_diff·a2m_scale) >> 13) + aln_cfg_xs_pair_penalty_)`. Final value = `min(xs_heur_mapq, pe_mapq − pen)`. Near-best pairs and singles are counted within a window one SNP cost wide (`list_pe_min = list_pe_max − similarity_.getSnpCost()`). Pair penalties use an insert-size lookup table (`petab_rom_c`) and `aln_cfg_unpaired_pen_`. The core `computeMapq`/`aln2mapq` are defined in another file not retrieved — [DRAGMAP PairBuilder.cpp](https://raw.githubusercontent.com/Illumina/DRAGMAP/master/src/lib/align/PairBuilder.cpp)

**Illumina Isaac**
- Final MAPQ is the alignment score truncated to 60 — [Isaac docs](https://support.illumina.com/help/BS_App_ENR_OLH_15050961/Content/Source/Informatics/IsaacAlignerDetermScore.htm)

**LAST (probabilistic)**
- `last-map-probs` gives each alignment a mismap probability, the probability that it is not the query's true origin. It assumes one reported alignment is the origin: P(first is origin) = exp(s1/t)/Σ_j exp(s_j/t). The default filter discards mismap > 0.01 (`-m`). The method is Frith, Wan & Horton 2010 NAR (quality-aware alignment) — [last-map-probs doc](https://home.cc.umanitoba.ca/%7Epsgendb/doc/last/last-map-probs.rst); [Debian source](https://sources.debian.org/src/last-align/1609-1/bin/last-map-probs)

**Stampy (Lunter & Goodson 2011)**
- MAPQ is an approximate Bayesian posterior that the inferred locus is wrong, computed from the likelihoods of candidate loci with base qualities. A separate P-value, estimated from nucleotide qualities, that the true candidate was never considered caps the posterior; if it is too high the mapping is not reported. In simulations Stampy, BWA and MAQ were well calibrated, while ELAND and Novoalign were over- or under-confident — [Stampy, Genome Research 21:936](https://genome.cshlp.org/content/21/6/936) (from a search summary; full text blocked)

**vg Giraffe (graph)**
- Developer statement: MAPQ is a phred-scaled probability of wrong placement, estimated from multiple candidate locations. Mismatches count "only to see how good the final alignment is relative to other potential alignments", so a read with several mismatches can still get MAPQ 60 if no rival exists — [Biostars 9581513](https://biostars.org/p/9581513)

**minibwa (Li & Homer, arXiv June 2026)**
- A hybrid of BWA-MEM variable-length seeding (ropebwt3-style SMEM finding) and minimap2 chaining with SIMD alignment. About 4× faster than BWA-MEM and over 2× faster than BWA-MEM2. The preprint gives **no MAPQ formula** — [arXiv 2606.15357](https://arxiv.org/html/2606.15357)

### Inferences
- **BWA-MEM default, worked numbers** (assuming defaults a=1, b=4, o=6, e=1, min_seed_len=19). For l=150 the length factor is ln50/ln150 = 0.781, so MAPQ ≈ 6.02·0.610·id⁴·Δ ≈ 3.67·id⁴·Δ, with Δ = score − sub in match units. One extra mismatch at the competitor (Δ=5) gives about 18, minus 3 for sub_n=1, so about 15. Two extra (Δ=10) give about 36. Three (Δ=15) give about 55. Four or more give 60. If the competing copy is not found or is missing from the reference, sub falls back to 19. Even a best hit with 8 mismatches (score 110, id ≈ 0.947, id⁴ ≈ 0.80) then gives 3.67·0.80·91 ≫ 60, so MAPQ is 60. The identity term therefore never bites for a unique 150 bp hit, and absolute divergence is effectively invisible to BWA-MEM/BWA-MEM2 MAPQ.
- **BWA-MEM pairing makes this worse.** With no alternative pair, o − subo ≥ pen_unpaired (17 by default), so q_pe ≥ 102, capped to 60. Each mate is then raised by up to 40 toward q_pe. A missing-copy pair that is unique in GRCh38 therefore gets 60 on both mates.
- **minimap2 sr.** Assuming match_sc=2, a perfect 150-mer gives 40·ln(150) ≈ 200, capped to 60. The identity factor is linear (5 mismatches give about ×0.97) and multiplies a value far above the cap. The (1 − x²) form is lenient: a competitor at 90% of the best score still gives about 40·0.19·5.0 ≈ 38. That is consistent with minimap2 being the most overconfident in our data (130/M), but it is not proven by it.
- **strobealign paired.** MAPQ = min(60, ΔSW) in raw score units. With its scoring (match 2, mismatch 8, assumed), a 6-mismatch gap reaches 60. A single candidate pair also gives 60.
- **Bowtie2 is the only aligner here with an absolute-divergence term.** Take end-to-end, L=150, default `--score-min L,-0.6,-0.6` so diff ≈ 90.6, and MX ≈ 5–6 per high-quality mismatch. A **unique** hit then gets 42 with ≤3 mismatches, 40 with about 4–5, 24 with about 6–7, 23 with about 8–9, 8 with about 10, 3 with about 11–12. With a distant second-best, any imperfect hit is capped at 33 instead of 39. Li 2014's observation that Bowtie2 gives low MAPQ to reads from paralogs missing from the reference (Q4) follows from this design.
- **CERTA's certified "unique by a margin of edits" is the same statistic all these aligners use** (Δ to the best competitor in the reference). Its proof makes CERTA confident exactly where the reference lacks the competitor. None of BWA, minimap2, strobealign or DRAGMAP has a term that would save it.

### Gaps
- Novoalign and GSNAP MAPQ formulas: no primary documentation found (Novocraft docs not retrieved). Novoalign is described as posterior-based only second-hand.
- minibwa's MAPQ code (lh3/minibwa source) was not retrieved. The preprint is silent on it. It likely inherits minimap2-style `mm_set_mapq2` given the shared chaining, but that is unverified.
- Current Rust strobealign MAPQ, minimap2's pair-level MAPQ (pe.c), and DRAGMAP's `computeMapq` body were not retrieved.
- Stampy's exact likelihood model and divergence prior: the PMC full text was blocked by a CAPTCHA. Only the abstract-level description is confirmed.
- Giraffe's MAPQ caps: I recall a cap based on the probability that distinguishing minimizers are explained by sequencing errors, but I could not verify it in a source.

## Q2. Which features best predict mismapping in segmental duplications and paralogs

### Takeaway
No source I found ranks per-read features specifically for segdup or missing-copy errors on real data. The evidence points to two features. The first is absolute divergence of the best hit, especially high-quality mismatches not explained by known variation: Bowtie2's absolute-score MAPQ resists missing-paralog artifacts, and reads carry 20–25% fewer mismatches on T2T-CHM13 than on GRCh38. The second is locus priors (known collapsed or missing-copy regions, excess heterozygosity). The ML recalibrators use AS, AS−XS, base-quality sums, soft-clip qualities and fragment length, and were trained on reference-drawn simulations that cannot contain missing-copy errors.

### Cited Findings
- Qtip's model features: best alignment score, best − second-best score difference, aligned length, sum of base qualities over aligned bases, sum of base qualities over soft-clipped bases, and fragment length for concordant pairs. Aligner-specific features are in Supplementary Note 6 and feature importance in Supplementary Note 12 (not retrieved) — [Qtip, Genome Biology 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
- LoQuM's logistic-regression features: base quality scores, alignment statistics such as match and mismatch counts, and the aligner's own MAPQ — [Ruffalo et al., Bioinformatics 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3436835)
- k-mer multiplicity is already in BWA-MEM (`frac_rep`, the fraction of the read covered by seeds with occurrence > max_occ, scales MAPQ by 1 − frac_rep) and in minimap2 (`uniq_ratio = sum_sc/(sum_sc + rep_len)` inside pen_s1) — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c); [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c)
- Number of near-optimal hits: BWA-MEM, minimap2 and BWA pair mode all subtract 4.343·ln(n_sub+1), and DRAGMAP subtracts about 3·log2(sub_count) in fixed point — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c); [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c); [PairBuilder.cpp](https://raw.githubusercontent.com/Illumina/DRAGMAP/master/src/lib/align/PairBuilder.cpp)
- Mismatch burden as a signal of missing copies: with 3,202 1KGP samples mapped by BWA-MEM, the per-read mismatch rate is 20–25% lower on T2T-CHM13 than on GRCh38, and 7.4×10⁶ (0.97%) more reads map properly paired — [Aganezov et al., Science 2022](https://pmc.ncbi.nlm.nih.gov/articles/9336181)
- Heterozygosity as a signal: 23.6% (87,005/368,574) of CHM13 heterozygous calls (HiFi/simulated) fall at gnomAD excess-heterozygosity sites, which the authors read as false positives from collapsed copies (permutation P = 1×10⁻⁴). Li 2014's Hardy–Weinberg / negative-inbreeding filter (P < 1e−10, clustered within 10 kb) was almost as effective as max-depth on CHM1 — [Aganezov 2022](https://pmc.ncbi.nlm.nih.gov/articles/9336181); [Li 2014](https://ar5iv.labs.arxiv.org/html/1404.0929)
- Depth as a signal: Li 2014 found that copy-number variants and paralogs missing from the reference cause many false heterozygotes, and the maximum-depth filter (about d+3√d to d+4√d) removes many of them — [Li 2014, Bioinformatics 30:2843](https://ar5iv.labs.arxiv.org/html/1404.0929)
- Absolute score versus relative gap: Li 2014 gives a pair near chr1:13.5 Mb with two mismatches where BWA-MEM gave MAPQ 27 and Bowtie2 gave 6, and only BWA-MEM produced the SNP call. Li reports that Bowtie2 is more likely to give low MAPQ to reads from paralogs missing from the reference — [Li 2014](https://ar5iv.labs.arxiv.org/html/1404.0929)
- Large SVs: the minibwa preprint states that "most false read alignment at high mapping quality is caused by large structural variants between individuals" — [minibwa, arXiv 2606.15357](https://arxiv.org/html/2606.15357)

### Inferences
- For our error class, the predictive features are those not conditioned on the reference's competitor set: (i) the number of mismatches and indels at the best hit, weighted by base quality (high-Q mismatches are unlikely to be sequencing errors); (ii) whether those mismatches sit at known population variant sites (gnomAD/HPRC) or are novel, since novel high-Q mismatches in a segdup look like paralogous sequence variants; (iii) a locus prior: segdup membership, the GRCh38 collapsed-duplication set (203 loci / 8.04 Mbp with missing copies per Aganezov), and gnomAD excess-het sites; (iv) whether the mate also lies in a segdup, or the pair's both-mate divergence. The gap to the second-best hit (CERTA's margin) is near-useless for this class by construction: the competitor is absent.
- Base qualities at mismatches matter more on NovaSeq X than they used to. With binned qualities, a mismatch at a high-Q bin is strong evidence of a real sequence difference. This should be checked empirically on our HG002 errors versus correct confident mates.

### Gaps
- No real-data per-feature ranking (AUC or importance) for segdup or missing-copy mismaps was found. Qtip's feature-importance supplement and LoQuM coefficient tables were not retrieved.
- No paper found that uses k-mer counts from the sample itself (sample-specific k-mer multiplicity above the reference count) as a per-read MAPQ feature for short reads.

## Q3. Machine-learned or recalibrated MAPQ (2012–2026): methods and gains

### Takeaway
Recalibration (LoQuM 2012, logistic regression; Qtip 2017, random forest on tandem simulations; a PeerJ 2020 recalibration study) improves ranking and calibration on simulated data substantially: Qtip cuts the area metrics by 7–69%. Real-data gains are small (Qtip: ΔF-β ≈ +0.001 in variant calling; PeerJ 2020: "modest improvements"). Qtip explicitly ignored reads from sequence missing from the reference, which is the error class that dominates our measurements. No DeepVariant-era ML-MAPQ method for short reads was found.

### Cited Findings
- **LoQuM (Ruffalo, Koyutürk, Ray, LaFramboise; Bioinformatics 28(18):i349, 2012):** logistic regression trained on simulated reads. Motivation: aligners report many accurate mappings at MAPQ 0 and many inaccurate ones at high MAPQ. LoQuM's estimates "strongly correlate" with accuracy for most aligners, recover many MAPQ-0 mappings, and improve SNP-call precision — [PMC3436835](https://pmc.ncbi.nlm.nih.gov/articles/PMC3436835)
- **Qtip (Langmead, Genome Biology 18:152, 2017):** trains a scikit-learn random forest on "tandem reads" simulated from the reference to mimic the input reads, with separate models for unpaired, concordant, discordant and bad-end alignments. Metrics are RCA (relative change in CID, a ranking metric) and RCE (relative change in CSED, a calibration metric); negative values favor Qtip — [PMC5557537](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
  - Bowtie2, GRCh38, simulated Mason reads, end-to-end: unpaired 150 nt RCA −9.11%, RCE −16.22%; paired 150 nt RCA −22.39%, RCE −25.94%; paired 250 nt −37.65%/−33.48%; paired 500 nt −54.59%/−44.91%. Local paired 150 nt: −46.87%/−52.76% — [Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
  - BWA-MEM, GRCh38: unpaired 100 nt −15.49%/−47.42%, 250 nt −15.78%/−51.14%; paired 100 nt −14.19%/−41.35%, 250 nt −12.36%/−42.78% (SD < 0.6). The large RCE gains mean BWA-MEM's MAPQ is poorly calibrated in absolute terms even on simulated data — [Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
  - SNAP, GRCh38 paired: RCA −51.36%/−51.32%, but RCE −11.16% at 100 nt and **+4.34% (worse)** at 250 nt — [Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
  - Real data: only ERR194147 (NA12878, Bowtie2, 100×100, Freebayes). F-β rose by about 0.0009–0.0011 for β = 0.5–4. For example, at β=1, F went from 0.9835 to 0.9845 (+3,999 TP, −2,745 FP), and at β=4, FP fell by 13,600. Qtip MAPQs spanned 0–36 versus 0–42 for Bowtie2 — [Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
  - Limits: "Category 1" errors (reads from sequence missing from the reference) were assumed rare and ignored. Adding CHM1-specific sequence had "little effect" (Table S2). Category 3 errors (wrong locus within the reference, typically repeats) made up 95.8–99.7% of simulated errors. Qtip was worse at many thresholds for BWA-MEM unpaired reads, especially in Z. mays — [Qtip](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
- **PeerJ 2020, "Recalibration of mapping quality scores in Illumina short-read alignments improves SNP detection results in low-coverage sequencing data":** applies a logistic-regression recalibration and reports that it "delivers modest improvements" in single-sample variant calling — [PeerJ 10501](https://peerj.com/articles/10501) (from a search snippet; the full text returned 403); [PMC7727374](https://pmc.ncbi.nlm.nih.gov/articles/PMC7727374/)

### Inferences
- A tandem-simulation recalibrator trained on reads drawn from GRCh38 cannot learn our error class: by construction every training read has its true copy in the reference. To learn "missing copy", training reads must be simulated from a genome that contains the copies GRCh38 lacks and then mapped to GRCh38. T2T-CHM13 (8.04 Mbp of missing-copy loci per Aganezov) or the HG002 Q100 haplotypes would serve, though for calibration the latter would leak into our held-out evaluation. This matches the direction of our latest commit ("a MAPQ informed by T2T-CHM13"). Per our tune-for-recent-Illumina rule, calibrate on HG002 and hold out HG001/HG005.
- Qtip's finding that adding CHM1-specific sequence had little effect was measured on simulated reads, where such reads are rare. On our real HG002 data the missing or diverged-copy class is nearly all of the confident errors, so Qtip's assumption does not transfer.

### Gaps
- No 2018–2026 ML-MAPQ method for Illumina short reads found (e.g., deep-learning MAPQ, or a DeepVariant-integrated MAPQ). DeepVariant consumes MAPQ as a channel but, per my sources, does not recalibrate it.
- LoQuM's and PeerJ 2020's quantitative calibration and SNP results were not retrieved (CAPTCHA/403).

## Q4. Evaluations of MAPQ calibration on real data (assembly/GIAB truth) and overconfidence in segdups

### Takeaway
I found no published per-read MAPQ calibration curve on real Illumina data against an assembly truth, stratified by segdups. Published calibration evidence is almost all simulated (Stampy, Qtip, minibwa, Giraffe). The best real-data proxies are variant-level. GRCh38 collapsed duplications push Illumina SNV precision to 64–67%, versus 98–99.7% on T2T-CHM13. False duplications push recall to 1–2%. CHM1 on GRCh37 without decoys gives about twice the false heterozygotes of hs37d5. Our locus-level, matched-sensitivity measurement on HG002 Q100 seems to fill a real gap.

### Cited Findings
- Simulated calibration only: minibwa plots the "fraction of correctly mapped at a mapping quality threshold" on mason 2×150 reads simulated from GRCh38, T2T-CHM13 and the HG002 haplotypes. It reports no numeric error rates and no real-data MAPQ calibration, and notes that simulated accuracy "may not reflect accuracy in downstream data processing" — [minibwa](https://arxiv.org/html/2606.15357)
- minibwa real-data downstream (HG002, GRCh38, DeepVariant 1.10, GIAB Q100 v1.1): SNP FN/FP were minibwa 46,367/7,544, BWA-MEM 46,895/7,585, strobealign 57,545/8,108. Indel FN/FP were 36,321/5,308, 37,425/5,218 and 37,828/5,655 — [minibwa](https://arxiv.org/html/2606.15357)
- Giraffe vs BWA-MEM MAPQ calibration is in the Giraffe paper's Supplementary Fig. 5, on simulated reads — [Biostars 9581513](https://biostars.org/p/9581513)
- Real-data, variant-level effect of GRCh38's reference gaps (Aganezov 2022). In a ~910 kbp subset of collapsed duplications, Illumina precision was 64.3–67.3% on GRCh38 versus 98.1–99.7% on T2T-CHM13, with recall 21.5–25.4% versus 20.1–28.3%. In corrected false duplications, recall was 1.1–1.8% (GRCh38) versus 57.4–68.3% (T2T) and precision 76.5–95.8% versus 98.5–99.3%. Variant density in 8 medically relevant genes in collapsed regions was 58/kbp (GRCh38) versus 40/kbp (T2T). HiFi false positives in 269 challenging medically relevant genes in HG002 fell by up to 12× on T2T. Mendelian-discordant variants (median per trio) were 8,879 versus 7,484 — [Aganezov et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/9336181)
- KCNE1: on GRCh38 many reads mismap to KCNE1B (a false duplication on chr21p), lowering KCNE1 coverage. This is the opposite failure, where an extra reference copy splits reads to MAPQ 0 — [Aganezov 2022](https://pmc.ncbi.nlm.nih.gov/articles/9336181)
- GRCh38 problem regions: 908 putative problematic regions cover 20.8 Mbp, 73.6% overlap segdups, and 203 loci (8.04 Mbp, 308 gene features) have copies missing in GRCh38 — [Aganezov 2022](https://pmc.ncbi.nlm.nih.gov/articles/9336181)
- Li 2014 (CHM1, haploid, so any heterozygous call is an error): mapping to the GRCh37 primary assembly gave about twice as many heterozygotes as hs37d5, which adds 35.4 Mb of decoys. Raw error was about 1 per 10–15 kb, and about 1 per 100–200 kb after filtering, leaving about 15–30k false heterozygotes per genome — [Li 2014](https://ar5iv.labs.arxiv.org/html/1404.0929)

### Inferences
- Our setup (GRCh38 without decoys or ALTs) is the configuration Li 2014 showed doubles reference-gap artifacts. Adding hs38DH-style decoys or ALTs would remove part of the missing-copy class at the reference level. But the HG002 Q100 truth projected to GRCh38 would then need decoy-aware scoring.
- Because aligners compute MAPQ from the reference's competitor set, overconfidence in segdups should concentrate in the 8 Mbp of missing-copy loci and in highly diverged paralogs. A stratified plot of our HG002 errors against Aganezov's collapsed-duplication BED would test this directly.

### Gaps
- No published real-data, assembly-truth per-read MAPQ calibration (error rate per MAPQ bin, stratified by GIAB segdup stratification) was found for BWA-MEM2, minimap2, strobealign, DRAGEN or minibwa.
- The paftools mapeval MAPQ-60 error rates surfaced in search (about 1.85e-3 for minimap2) are for long reads and were not used.

## Q5. Methods that lower MAPQ when the best hit has excess divergence (a "reference-gap" or "paralog" term)

### Takeaway
No widely used short-read aligner has an explicit term for "this read may come from a copy not in the reference". The closest existing mechanisms are four. Bowtie2's absolute-score MAPQ buckets. BWA-MEM's identity factors, which are inert for 150 bp unique hits in the default path. Stampy's P-value cap that the true locus was never considered. LAST/Stampy-style posteriors, which could carry a pseudo-candidate for an unseen copy. Building such a term is open territory and matches CERTA's need.

### Cited Findings
- Bowtie2 lowers MAPQ as `bestOver/diff` falls, i.e. as the best hit's score approaches the minimum valid score, even with no second-best (42 → 40 → 24 → 23 → 8 → 3) — [unique.h](https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h). Li 2014 reports that Bowtie2 is less susceptible to missing-paralog artifacts for this reason — [Li 2014](https://ar5iv.labs.arxiv.org/html/1404.0929)
- BWA-MEM: the default path multiplies by identity⁴ (via tmp²). The legacy path multiplies by identity² only when identity < 0.95 — [bwamem.c](https://raw.githubusercontent.com/lh3/bwa/master/bwamem.c)
- minimap2 multiplies MAPQ linearly by identity = mlen/blen — [hit.c](https://raw.githubusercontent.com/lh3/minimap2/master/hit.c)
- Stampy caps the posterior by an estimated P-value that the true candidate was not considered, computed from nucleotide qualities. It addresses missed candidates within the reference, not absent sequence — [Stampy](https://genome.cshlp.org/content/21/6/936)
- LAST's mismap probability is a softmax over alignment scores, exp(s_i/t)/Σ exp(s_j/t), so one more term in the denominator lowers confidence — [last-map-probs](https://home.cc.umanitoba.ca/%7Epsgendb/doc/last/last-map-probs.rst)
- Decoys as a reference-level "missing copy" sink: hs37d5's 35.4 Mb of decoys halved false heterozygotes for CHM1 — [Li 2014](https://ar5iv.labs.arxiv.org/html/1404.0929)

### Inferences (candidate designs for CERTA; not from sources)
- **Unseen-copy pseudo-hit.** Use a posterior in the LAST/Stampy style: P(wrong) = [π_miss(locus)·L_miss(read)] / [L_best + Σ L_alt + π_miss·L_miss].
  - L_best is the likelihood of the read at the certified hit under sequencing error and population heterozygosity, using base qualities and known-variant sites.
  - L_miss is the likelihood under a "diverged paralog" model: the same hit plus an extra substitution rate d (about 0.5–3% for recent segdups), integrated over d.
  - π_miss is a locus prior: high in GIAB segdups and in Aganezov's 8.04 Mbp missing-copy loci, about 0 elsewhere.
  - The effect is that MAPQ falls with each high-Q, non-population mismatch at the best hit, but only where a missing copy is plausible. This is the term BWA/minimap2/strobealign/DRAGMAP lack and Bowtie2 approximates crudely, without a locus prior.
- **Margin × divergence interaction.** CERTA's certified edit margin m still matters, but conditioned on k, the number of edits at the best hit. A read with k=0 and m=5 is safe. A read with k=4 high-Q novel mismatches and m=5 in a segdup is the typical missing-copy error. A simple feature is k_highQ_novel − expected(k | error rate, heterozygosity), thresholded or fed to a calibrator.
- **Pair-level.** Use the sum of both mates' excess divergence when both lie in segdups. Do not let a unique mate lift a diverged mate's MAPQ the way BWA's "+40" rule does.
- **Calibration data.** Simulate from T2T-CHM13 or the HG002 haplotypes and map to GRCh38, so training includes reads whose copy is absent. Fit on HG002 real data at matched sensitivity and validate on HG001/HG005.

### Gaps
- No published per-read "reference-gap" or excess-divergence MAPQ term for short reads was found. Paralogous-sequence-variant methods such as DuploMap and Parascopy (Prodanov & Bansal) appeared only as long-read or copy-number tools in search results and were not verified for short-read MAPQ.
- No source quantifies the expected-divergence baseline (sequencing error at NovaSeq X binned qualities plus HG002 heterozygosity) per 150 bp read. It must be measured from our own data.
