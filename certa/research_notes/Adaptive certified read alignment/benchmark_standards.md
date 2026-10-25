# Methodological standards for benchmarking short-read aligners, and for validating tuned or adaptive parameters without overfitting

Scope: what reviewers and the community expect from the evaluation in an aligner paper (truth sets, accuracy metrics, speed and memory measurement, replicates and variance, determinism, artifact availability, calibration vs held-out validation, validation of certified claims), applied to CERTA. CERTA is an exact, certified short-read aligner. It is calibrated on GIAB HG002, held out on HG001 and HG005, run on ARM64 + A100 and on x86-64, reported as medians of 3 runs, and compared with minibwa, BWA-MEM2, strobealign, minimap2, Bowtie 2, Yara and RazerS 3. Run-time adaptive parameters are planned.

Method note: Hoefler & Belli's rules and survey numbers come from text extracted from the authors' own PDF. Some sources could be read only in part, because of PMC reCAPTCHA pages and Springer redirects. Where that happened, the finding says so.

---

## Q1. What general and aligner-specific benchmarking guidelines require: neutrality, fair parameter tuning, datasets, metrics

### Takeaway
The core guidelines agree on four points:
- Use one tuning strategy for every tool, either all defaults or equal tuning effort. Tuning your own method more is a known source of bias.
- Include all relevant methods, or use inclusion criteria that favour none, and document every configuration.
- Use several datasets and metrics, because no single winner usually exists.
- Publish code, parameters and versions.

Papers that introduce a method are measurably over-optimistic compared with neutral benchmarks. A non-neutral paper like CERTA's has to pre-empt that criticism explicitly.

### Cited Findings
**Weber et al. 2019, "Essential guidelines for computational method benchmarking", *Genome Biology* 20:125**
- Neutrality: the benchmarking group should know all included methods about equally well, or involve the original authors. Methods whose authors decline to take part should still be reported. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Method selection: a neutral benchmark should include all available methods. Otherwise the inclusion criteria must favour no method, and leaving out a widely used method must be justified. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Three tuning strategies are acceptable: defaults, published or expert values, or systematic tuning such as a grid search. Comparing one tool at defaults with another that was tuned extensively is unfair, so the same strategy should be used for every tool. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Weber et al. name tuning a new method more extensively than its competitors as a source of bias, and recommend reporting the tuning strategy transparently. Final results should generally use the latest software versions. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Datasets:
  - Simulated data give ground truth, but must be checked against the properties of real data. Real data often lack ground truth.
  - Spike-ins, validation experiments and semi-simulated data are alternatives.
  - Diversity of evaluations is what makes a benchmark excellent. Datasets should be neither so easy nor so hard that they cannot separate methods.
  - [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Metrics and interpretation:
  - Use rankings from several metrics to find a set of consistently strong methods.
  - Characterise stability and robustness.
  - Expect no single winner; present the trade-offs instead.
  - Discuss the benchmark's limitations openly.
  - [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Secondary measures (runtime, memory, scalability, usability) vary widely and should be measured systematically. Performance often trades off against computational cost. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Reproducibility:
  - Report parameter values, including random seeds, and software versions.
  - Consider containers.
  - "Available on request" is often not enough; use public, archival repositories.
  - [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)

**Mangul et al. 2019, "Systematic benchmarking of omics computational tools", *Nature Communications* 10 (doi 10.1038/s41467-019-09406-4)**
- Compile a comprehensive list of tools, since choosing the best tool in advance is risky. Tools dropped because they would not install or ran too long should be logged. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Parameter optimization can substantially improve accuracy. In one example, tuning RNA-seq aligners raised the share of correctly mapped reads by about 10% on average. For RNA-seq alignment, the most influential parameters were the number of allowed mismatches and the seed length. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Benchmarks should report server specifications and the number of processors used. Only 28% of surveyed studies reported runtime and memory. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Only 40% of surveyed studies shared their data fully. Mangul et al. recommend sharing all data and code, the raw tool outputs (so readers can apply their own metrics), and scripts that reproduce the figures. Only one surveyed study (CAMI) packaged its tools as containers. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Every read-alignment benchmarking study in their survey used simulated data. They argue that simulated data can only complement real gold-standard data, because it is less complex than real data. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)

**Buchka et al. 2021, "On the optimistic performance evaluation of newly introduced bioinformatic methods", *Genome Biology* (doi 10.1186/s13059-021-02365-4)**
- In their case study (normalization methods for 450K methylation arrays), the new method beat its predecessors much more often in the papers that introduced it than in neutral comparisons:

  | Level | Introducing papers | Neutral comparisons |
  |---|---|---|
  | Paper | 94.7% (18 of 19) | 64.3% (49.5 of 77) |
  | Substudy | 83.2% | 61.2% |

  [Buchka et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8111726/)
- Sources of over-optimism they list:
  - selective reporting of datasets, metrics and the best-performing method variant
  - debugging bad results of one's own method more carefully than competitors'
  - choosing one's own parameters more attentively
  - testing on familiar datasets

  [Buchka et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8111726/)
- Recommendations: more neutral benchmarks, reviewers who accept a method that is not best on every measure, and open reporting of weaknesses. Preregistration is not explicitly recommended in the paper. — [Buchka et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8111726/)
- Nießl et al. (*WIREs Data Mining and Knowledge Discovery* 12(2) e1441) show that many design and analysis options, combined with questionable research practices, bias benchmarks toward over-optimistic conclusions. The options include the choice of datasets and performance measures, the handling of missing performance values, and how results are aggregated across datasets. — [Nießl et al. 2021/2022 (arXiv 2106.02447)](https://arxiv.org/abs/2106.02447)

**Heiser, "Systems Benchmarking Crimes" (web essay, systems community)**
- Heiser lists these practices as benchmarking crimes:
  - Selective benchmarking: showing only the cases that help, subsetting a suite without justification, or picking a flattering input range.
  - Unfair treatment of competitors: running them sub-optimally (for example with debugging enabled) or not documenting their configuration. Heiser says this can amount to misconduct.
  - Evaluating only against yourself.
  - Reporting only overall scores, which hides individual degradations.
  - Reporting only relative numbers.
  - Ignoring the current state of the art.

  [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)
- Baselines should be state-of-the-art, optimal or theoretical, or an unperturbed system. — [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)

**Aligner-specific benchmarks**
- Rabema (Holtgrewe et al. 2011, *BMC Bioinformatics* 12:210) defines the read-mapping problem formally. A feasible match has distance ≤ k, edit or Hamming. Three categories are defined over feasible matches:
  - *all*: every feasible match
  - *all-best*: every match with the smallest distance for the read
  - *any-best*: one best match

  Matches are merged into equivalence classes ("intervals") by neighbour and trace equivalence. — [Holtgrewe et al. 2011 (Europe PMC full text)](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)
- Rabema gold standard and metric:
  - A fully sensitive mapper (RazerS) supplies at least one match per class, and the intervals are extended with Myers' bit-vector algorithm.
  - Each read counts at most one point: if a read has n intervals, each found interval contributes 1/n ("normalized found intervals").
  - This scheme avoids over-weighting reads in long tandem repeats.

  [Holtgrewe et al. 2011](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)
- Rabema found that different parameter settings of the same tool can shift sensitivity by more than 20 percentage points. Soap2 and Bowtie at defaults were tuned for few errors, and their sensitivity dropped sharply beyond 2 errors on long Illumina reads. — [Holtgrewe et al. 2011](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)
- In the `rabema_evaluate` tool, the default category is *all* and the default metric is edit distance. An extra valid hit that the gold standard lacks indicates a flaw in the gold standard, not an error by the mapper. — [rabema_evaluate man page](https://dyn.manpages.debian.org/testing/seqan-apps/rabema_evaluate.1)
- Hatem et al. 2013 (*BMC Bioinformatics* 14:184) ran each mapper with its own defaults, using 1M wgsim-simulated 125-bp reads, among other tests. They also ran BWA with Bowtie's defaults ("BWA-ND") to separate the effect of the algorithm from that of the defaults. They concluded that no tool beat all others in every test, so users must state their needs. — [Hatem et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3694458/)
- Alser et al. 2021 ("Technology dictates algorithms", *Genome Biology* 22:249):
  - Benchmarked 11 Bioconda-installable DNA aligners on 10 randomly chosen 1000 Genomes WGS samples, all at default parameters.
  - Recorded CPU time and RAM.
  - Compared algorithm classes with likelihood-ratio and Wald tests, adjusting for publication year and algorithm features. After adjustment, BWT-FM tools needed on average 3.8× fewer computational resources than hashing-based tools.
  - Put all data and code on GitHub, and had tool authors verify the information about their tools.

  [Alser et al. 2021 (Europe PMC full text)](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8390189/fullTextXML)
- Strobealign (Sahlin 2022, *Genome Biology*, doi 10.1186/s13059-022-02831-7) is a direct precedent for CERTA's comparators:
  - All aligners, strobealign included, ran at defaults. The paper acknowledges that defaults may not be optimal, and reports one non-default strobealign setting (`-M 40`) separately.
  - Alignment used 16 threads, with thread scaling tested at 4, 8 and 16. Indexing used 1 thread.
  - Peak memory was reported for both indexing and alignment, and exact tool versions were listed.
  - Accuracy was defined as overlap of an alignment with the true genomic location on simulated reads.
  - Variant calling used bcftools on simulated data and on GIAB HG004 reads.

  [Sahlin 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264/)

### Inferences
- CERTA's paper is an introducing, non-neutral comparison in Buchka's sense, so reviewers will look for the bias sources Buchka lists. The strongest defence is a single, documented tuning policy:
  - every tool at its defaults, plus
  - one stated speed or sensitivity preset per tool, for example minimap2 `-x sr` and Bowtie 2 `--very-fast/--sensitive`,
  - with CERTA's default and adaptive modes reported side by side.
- Yara and RazerS 3 are fully sensitive or exact mappers. They are natural "optimal baselines" in Heiser's sense, and can serve as Rabema-style gold-standard generators that CERTA's exactness claim can be checked against.
- Rabema's all-best category matches "exact/certified" claims well (all co-optimal hits within k). CERTA can report all-best and any-best sensitivity on simulated reads, which most heuristic aligners cannot reach by construction.

**CERTA checklist for Q1**
- [ ] State one tuning policy for all tools. Give the exact command line, version, commit or conda build for each, and put them in the supplement.
- [ ] Justify any tool left out (for example SNAP, Accel-Align, BWA-MEM) and log failures to install or run.
- [ ] Give minibwa a full citation, version and a description. It is the least-known comparator, and reviewers will ask what it is.
- [ ] Report several metrics per tool (speed, memory, mapping accuracy, variant-calling F1 per stratum). Show which tool wins where, rather than declaring one winner.
- [ ] Show where CERTA is slower or worse, to avoid the "selective benchmarking" crime: short or long reads, high-divergence reads, small thread counts.
- [ ] Ask the authors of the competing tools to check your configurations, as Alser et al. did.

### Gaps
- I did not find a "Sarwal et al." paper benchmarking short-read aligners. The only Sarwal benchmarking paper I know of covers deletion SV callers (Brief Bioinform 2022), and I did not verify it. The coordinator should confirm which Sarwal paper is meant.
- Alser et al.'s hardware, thread counts and accuracy metrics are in Supplementary Note 3, which I could not read.
- I could not open Mangul et al.'s own numbered list of principles in full. The page was truncated and only paraphrases were extracted.

---

## Q2. Truth sets: GIAB benchmarks, hap.py/vcfeval with stratifications, held-out samples, simulated reads with known origin

### Takeaway
The accepted evaluation for aligners is two-pronged:
- Simulated reads with known origin, to measure mapping correctness directly. MAPQ-stratified error curves are the de-facto format.
- Downstream germline small-variant calling, measured against GIAB v4.2.1 (or newer T2T/Q100-based) benchmarks with GA4GH hap.py + vcfeval.

Results should be stratified by GIAB difficult regions and carry confidence intervals. Calibrating on HG002 and evaluating on other GIAB genomes mirrors how the community itself separates training data from test data. One case is precisionFDA Truth Challenge V2, where HG002 truth was public and HG003/HG004 were blinded. Another is DeepVariant, which withholds samples and chr20 from training.

### Cited Findings
**GA4GH best practices (Krusche et al. 2019, *Nature Biotechnology*)**
- Counts and labels:
  - TP, FP and FN require matching genotypes by default.
  - Allele-only matches are tallied separately as FP.GT or FP.AL, and each counts as both an FP and an FN.
  - Query calls outside the truth BED are labelled UNK.
  - True negatives and specificity are deliberately not used, because the number of true negatives is ill-defined for indels and complex variants.

  [Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)
- Definitions: precision = QUERY.TP/(QUERY.TP + QUERY.FP), and recall = TRUTH.TP/(TRUTH.TP + TRUTH.FN). Matching tiers, from strictest: genotype match (the default), allele match, local match. — [Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)
- hap.py implements the framework and produces a standardized report. vcfeval is the recommended comparison engine because it gives the most partial credit for matches. Stratification should be applied after comparison, so that variants represented differently are still matched. — [Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)
- Box 1 recommendations:
  - Use benchmark sets that include both calls and confident regions.
  - Stratify by variant type and by genomic region.
  - Compute confidence intervals for precision and recall.
  - Manually curate a subset of FPs and FNs, and report suspected errors in the benchmark.
  - Supplement with simulation and other approaches.

  [Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)
- Caveats:
  - Metrics mostly reflect easier regions.
  - The ranking of methods can flip with matching stringency, or between confident regions and the rest.
  - ALT loci in GRCh38 are poorly handled by most aligners and benchmarking tools.

  [Krusche et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/)

**GIAB benchmark sets**
- v4.2.1 (Wagner et al. 2022, *Cell Genomics*) used linked and long reads to extend benchmarks for seven samples into difficult-to-map regions and segmental duplications:
  - more than 300,000 extra SNVs and 50,000 extra indels
  - 16% more exonic variants
  - for HG002, 92% of the autosomal GRCh38 assembly, up from 85%

  Segmental duplications and low-mappability regions still make up 56.4% of the GRCh38 bases not covered. — [Wagner et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9706577/); [bioRxiv version](https://www.biorxiv.org/content/10.1101/2020.07.24.212712v5.full-text)
- A separate HG002 benchmark, CMRG, covers 273 challenging, medically relevant genes. It is built from a trio-based hifiasm assembly and covers genes that v4.2.1 left below 90% coverage. — [Bio-IT World report on Wagner et al., Nat Biotechnol 2022](https://www.bio-itworld.com/news/2022/02/14/giab-establishes-new-benchmark-for-several-medically-relevant-difficult-to-sequence-genes)
- Stratifications (Dwarshuis et al. 2024, *Nature Communications*, doi 10.1038/s41467-024-53260-y):
  - BED files that mark distinct genomic contexts, for GRCh37, GRCh38 and T2T-CHM13, generated by a Snakemake pipeline.
  - CHM13 adds many hard-to-map and GC-rich regions, and benchmarking against CHM13 carries a measurable performance penalty.

  [Dwarshuis et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11489684/)
- T2T/Q100:
  - GIAB and T2T are building near-perfect "Q100" diploid HG002 assemblies and associated benchmarks. — [Olson, GIAB update slides 2023](https://www.slideshare.net/slideshow/2023-giab-amp-update/263498422)
  - NIST lists HG002 v1.0 benchmarks for chrX/Y small variants and for tandem-repeat indels and SVs on GRCh38, plus preliminary v1.1 drafts still under evaluation. — [NIST GIAB program page](https://www.nist.gov/programs-projects/genome-bottle)
  - A news report (August 2026) says v1.1 adds 701.4 Mb (11.7%) of autosomal sequence and 216.8 Mb of sex-chromosome sequence missing from v4.2.1. This is secondary and needs verifying against the primary release. — [News-Medical 2026](https://www.news-medical.net/news/20260810/New-human-genome-benchmark-captures-both-parental-copies-with-near-perfect-accuracy.aspx)

**Held-out samples as community practice**
- precisionFDA Truth Challenge V2 (May–June 2020; Olson et al. 2022, *Cell Genomics*):
  - The HG002 v4.1 benchmark was public. The HG003 and HG004 benchmarks were withheld and served as the blind truth.
  - Scores were F1 for SNVs and indels combined, averaged over HG003 and HG004, across 64 submissions.

  [NIST/precisionFDA](https://www.nist.gov/node/1838346); [Olson et al. 2022 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9205427); [Bio-IT World 2020](https://www.bio-itworld.com/news/2020/09/02/truth-challenge-v2-latest-challenge-results-from-genome-in-a-bottle)
- DeepVariant v0.5 formalised a benchmarking-compatible training scheme: withhold an entire sample (HG002) and all chr20 data from training. — [Google Research blog](https://www.research.google/blog/deepvariant-accuracy-improvements-for-genetic-datatypes/)
- A DeepVariant metrics page in a GitHub fork reports HG003 results against v4.2.1 and describes HG003 as held out. This conflicts with the v0.5 statement above, so which sample is held out appears to depend on the version. — [DeepVariant metrics doc (fork)](https://github.com/listenw1/deepvariant/blob/r1.6/docs/metrics.md)

**Simulated reads with known origin**
- ART (Huang et al. 2012, *Bioinformatics* 28(4):593–594) simulates reads with empirical, technology-specific error and quality models, and was designed for benchmarking alignment, assembly and variant discovery. — [Huang et al. 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3278762)
- A comparison of six short-read simulators (Mason 2.0.9, ART, DWGSIM, InSilicoSeq, NEAT, wgsim) is available. Simulators differ in how realistic their errors are. — [Heredity, performance evaluation of six short-read simulators](https://pmc.ncbi.nlm.nih.gov/articles/PMC9905089/)
- Rabema notes that the genome itself defines how far from the simulated position a read may legitimately map. The simulated origin is a representative, not necessarily the unique optimum, which is why Rabema scores equivalence classes. — [Holtgrewe et al. 2011](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)
- minimap2's `paftools.js mapeval`:
  - Scores SAM/PAF output from simulated reads whose true position is encoded in the read name, via `mason2fq` or `pbsim2fq`.
  - Each "Q" line reports a MAPQ threshold, the reads in that MAPQ interval, the wrong mappings, the cumulative error rate and the cumulative number of mapped reads.
  - Papers often report the error rate on a phred scale (−10·log10 e).

  [minimap2 misc README (mirror)](https://lovelace.cluster.earlham.edu/mounts/lovelace/software/minimap2/misc/README.md); [Li 2021, arXiv 2108.03515](https://arxiv.org/pdf/2108.03515)

### Inferences
- **Calibrate on HG002, evaluate on HG001 and HG005** follows Truth Challenge V2 (public HG002, blind parents) and DeepVariant (held-out sample plus chr20). Two refinements would strengthen it:
  - HG001 and HG005 are not from the HG002 family, so the test also checks generalisation across ancestry: HG001 is CEU and HG005 is Han Chinese (well-known GIAB facts, not re-sourced here). Say so.
  - Benchmark coverage differs per sample under v4.2.1. Report the size of each sample's benchmark BED and its version, so readers can tell whether gains depend on the sample.
- **Watch for leakage through the caller.** If the downstream caller is a learned model (for example DeepVariant), its training may have included HG001, HG005 or HG002 outside chr20, depending on the version. That would compromise the hold-out at the variant-calling stage. Use a non-learned caller (bcftools or GATK HaplotypeCaller), or restrict the learned-caller evaluation to that model's held-out chromosome or sample, and state the caller's version and training set. Keep the caller and its parameters fixed across aligners, because variant F1 measures the aligner and caller together.
- **Accuracy reporting for an exact aligner** should have three layers:
  1. Simulated reads, with mapeval-style MAPQ curves, scored against simulated origin.
  2. Rabema-style all-best and any-best sensitivity against an exhaustive gold standard, which is the natural test of "certified complete/optimal".
  3. GIAB hap.py/vcfeval results, stratified (low mappability, segmental duplications, tandem repeats and homopolymers, GC extremes, "all difficult" and its complement), with confidence intervals.

**CERTA checklist for Q2**
- [ ] Name the GIAB benchmark version (v4.2.1, or T2T-Q100 v1.x if used), the reference build (GRCh38 with or without ALT; decoys), the stratification version, and the hap.py and vcfeval versions and flags (`--engine=vcfeval`).
- [ ] Report SNV and indel precision, recall and F1 separately for each sample, stratified, with confidence intervals (for example binomial or Wilson on TP/FN counts).
- [ ] Keep the downstream caller and its settings identical for all aligners, and document the caller's own training data if it is learned.
- [ ] Simulated reads: name the simulator and version, random seed, error and variant model, read lengths, and single- or paired-end. Define "correct" (overlap with origin, or within w bp, or same equivalence class). Report MAPQ-stratified error curves for every tool.
- [ ] Report Rabema all-best and any-best results with Yara or RazerS 3 as the gold-standard generator. This shows exactness empirically, not only in theory.
- [ ] Give the specifics of the sequencing data for HG001, HG002 and HG005 (accession, instrument, read length, coverage, any downsampling).

### Gaps
- I did not find a GIAB document that explicitly recommends calibrating on one sample and testing on another for aligner evaluation. The supporting precedent comes from the precisionFDA challenge and DeepVariant's training practice.
- I did not confirm which samples each DeepVariant release holds out. The official repo's documentation was not read, and the sources conflict (HG002 in v0.5, HG003 in one fork's doc).
- I did not verify HG001 and HG005 v4.2.1 coverage figures (the share of the genome inside each benchmark).
- The Mason technical report (Holtgrewe 2010, FU Berlin) could not be retrieved for exact citation details.

---

## Q3. Reporting speed and memory: wall vs CPU time, peak RSS, units, pinning, warm-up, replicates, statistics, hardware disclosure, artifact reproducibility

### Takeaway
The performance-reporting literature (Hoefler & Belli SC'15, Heiser's benchmarking crimes) asks for:
- absolute values alongside every ratio
- an explicit statement of whether measurements are deterministic, with nonparametric confidence intervals for the median when they are not
- no averaging of ratios with the arithmetic mean
- unambiguous units, including binary prefixes for GiB
- disclosure of warm or cold cache state
- full hardware and software documentation

Three runs, the median and the range is a common minimum. Nonparametrically, though, three runs cannot support a 95% confidence interval for the median, because the min–max range covers the median only 75% of the time.

### Cited Findings
**Hoefler & Belli 2015, "Scientific Benchmarking of Parallel Computing Systems: Twelve Ways to Tell the Masses when Reporting Performance Results", SC15 (doi 10.1145/2807591.2807644)**

The twelve rules, paraphrased from the authors' PDF — [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf):

| Rule | Paraphrase |
|---|---|
| 1 | For a parallel speedup, say whether the base case is one parallel process or the best serial run, and give the base case's absolute performance. Generalised: never report a ratio without the absolute values. |
| 2 | Justify reporting only a subset of benchmarks, or not using all system resources. Corollary: report all results, not just the best. |
| 3 | Use the arithmetic mean only to summarise costs; use the harmonic mean for rates. |
| 4 | Do not summarise ratios. Summarise the costs or rates behind them, and use the geometric mean only if those are unavailable. |
| 5 | Say whether the measurements are deterministic, and give confidence intervals for nondeterministic data. |
| 6 | Do not assume normality, for example because the sample is large, without checking it. |
| 7 | Compare nondeterministic data soundly, for example with non-overlapping confidence intervals or ANOVA. |
| 8 | Check whether mean or median is meaningful at all; some questions need other percentiles. |
| 9 | Document every varying factor and its levels, and the complete hardware and software setup. |
| 10 | For parallel timings, report the measurement, synchronisation and summarisation methods. |
| 11 | Show upper performance bounds where possible. |
| 12 | Plot enough to interpret the results, and connect points with lines only where interpolation is valid. |

Further points from the same paper:
- Survey of 120 papers from three top HPC conferences:
  - Most reported nondeterministic data, but the text says only 15 mentioned any measure of variance. The summary table gives "Measure of Variation" as 17/95, a small internal inconsistency.
  - Only 2 of 95 reported confidence intervals around the mean.
  - Processor model was reported by 79/95 papers, but RAM size, type and bus information by only 26/95.

  [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)
- Units: follow IEC 60027-2 and use binary prefixes (MiB, GiB), or state the base. Only 2 of 95 papers had fully unambiguous metrics. The paper also warns against cherry-picking, for example not using the whole node. — [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)
- Medians and robustness:
  - Computer-performance distributions are often multi-modal and heavily right-skewed, because of scheduling, congestion and cache misses.
  - Nonparametric statistics such as the median and other percentiles are the most robust, but they still need i.i.d. measurements.
  - Nonparametric CIs of the median follow Le Boudec's rank formula and may be asymmetric.
  - Data that are nearly deterministic can be summarised, for example by measuring until the 99% CI lies within 5% of the reported value.
  - The required number of measurements can be computed from the CI width.

  [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)
- Warm-up and caches: leave the first iteration out where programs build up working state on demand. Warm versus cold cache matters, and should match the intended use. — [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)

**Heiser, benchmarking crimes**
- Further crimes:
  - giving no indication of significance (no standard deviations or tests)
  - averaging normalised scores with the arithmetic mean instead of the geometric mean
  - confusing percent with percentage points, or choosing a convenient denominator
  - giving relative numbers only
  - leaving out the platform: CPU model, cores, clock, memory and cache sizes, OS, compiler versions

  [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)
- Heiser's measurement practice:
  - Run several times and check the standard deviation. In his systems context, typical deviations should be below 0.1%, and above 1% should raise concern.
  - Measure each point repeatedly, in varied order.
  - Drop outliers only after formal detection, and disclose it.

  [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)

**Bioinformatics-specific practice**
- Mangul et al.: report server specifications and processor counts; only 28% of benchmarks reported runtime and memory. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Aligner papers vary in how much they disclose:
  - Strobealign 2022 reports 16-thread alignment, single-thread indexing and peak memory, but the visible text does not say whether times are wall or CPU, how many runs were made, or which timing tool was used. — [Sahlin 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264/)
  - Alser et al. report CPU time and RAM over 10 samples. — [Alser et al. 2021](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8390189/fullTextXML)

**Artifact evaluation (ACM badging v1.1, effective 15 May 2020)**
- Results Replicated: another team validated the main results, at least partly using the authors' artifacts.
- Results Reproduced: the same, without any of the authors' artifacts.
- ACM swapped the two terms in 2020 to match NISO usage; badges issued up to 14 May 2020 follow v1.0.
- The badges, including Artifacts Available and Artifacts Evaluated (Functional and Reusable), are independent of each other.

[ACM badging terms](https://www.acm.org/publications/badging-terms); [ACM v1.0 (superseded)](https://www.acm.org/publications/policies/artifact-review-badging/)

### Inferences
- **Medians of three.** With n = 3 i.i.d. runs, P(min ≤ median ≤ max) = 1 − 2·(1/2)³ = 75%. This is order-statistics arithmetic I derived, not a quote. A nonparametric 95% CI for the median needs n ≥ 6 (1 − 2/64 ≈ 96.9%), and a 99% CI needs n ≥ 8. CERTA can either:
  - (a) keep n = 3, report median with min–max, say explicitly that this is a range and not a 95% CI, and add the coefficient of variation; or
  - (b) raise n to at least 6 for the headline comparisons, and give nonparametric 95% CIs for the median (Le Boudec, as cited by Hoefler).

  Either way, Rule 5 asks for a sentence on whether timings are deterministic. Peak RSS for a deterministic program is often nearly so.
- **Wall and CPU time.** Wall time is the user-facing quantity for multithreaded aligners, so it should be the headline. CPU time (user + sys) should also be reported, to show parallel efficiency and to allow fair comparison across thread counts.
  - Separate index or model loading from mapping, since CERTA's recent commits changed how the index directory is loaded.
  - For the A100 runs, count GPU initialisation and host-device transfers in end-to-end time, and report GPU memory separately from host RSS.
- **Peak memory.** Use one tool for every aligner (for example `/usr/bin/time -v`, "Maximum resident set size", in KiB), convert to GiB, and label it GiB as Hoefler's units section asks.
  - Say whether file-backed mmap pages count. RSS includes resident file-backed pages, so mmap-heavy aligners can look different under a cold page cache.
  - The repo's current convention ("the median run's peak") is reasonable, but should be stated in the paper.
- **Page cache.** State whether runs are warm (index already in the page cache after a discarded warm-up run) or cold (caches dropped before each run). Ideally report both for index-load-dominated tools, since warm-cache runs favour tools with large mmap'd indexes. Dropping caches via `/proc/sys/vm/drop_caches` is standard Linux practice; I did not source it here.
- **Pinning and environment.** Document:
  - thread count per tool (and why it is the same for all)
  - SMT on or off
  - NUMA policy (`numactl`/`taskset`)
  - CPU frequency governor and turbo
  - compiler, flags and SIMD level for each tool, for example BWA-MEM2's AVX-512 build on x86 versus NEON on ARM64
  - OS and kernel
  - CUDA driver and toolkit
  - storage type (input on local NVMe or network FS)
  - whether the machine was otherwise idle
- **Summaries across datasets.** Per Rule 4, do not average speedups arithmetically. Summarise the underlying times, or use the geometric mean of the ratios, and always show the absolute seconds.
- **Upper bounds (Rule 11).** Possible bounds for CERTA: the time to stream the FASTQ and write the SAM (I/O floor), or exact-match-only throughput. Either helps interpret how close CERTA is to the limit.

**CERTA checklist for Q3**
- [ ] Hardware table for both servers: CPU model, core and thread counts, base and boost clock, cache sizes, RAM size, type and channels, NUMA nodes, GPU model and memory (A100 40 or 80 GB), storage. Software table: OS, kernel, compiler and flags, CUDA, versions of every tool.
- [ ] Wall-clock time (primary) and CPU time, index loading separated, threads stated, the same thread count for all tools, plus a thread-scaling plot.
- [ ] Peak RSS in GiB, all tools measured with the same tool, page-cache state stated, and GPU memory reported separately.
- [ ] Number of runs, run order (interleaved or randomised across tools), warm-up policy, median plus range or nonparametric CI, and a statement on determinism (Rule 5).
- [ ] Absolute numbers next to every speedup; no arithmetic means of ratios.
- [ ] Results for every dataset and configuration, including those where CERTA loses (Rule 2 corollary).
- [ ] An artifact (code, exact commands, container or conda lockfile, raw timing logs, plotting scripts) with an archival DOI (for example Zenodo), aiming at ACM "Artifacts Available" and "Evaluated – Reusable" quality even for a journal.

### Gaps
- I found no aligner-specific standard on page-cache state, thread pinning or the number of replicates. These recommendations are carried over from HPC and systems guidance.
- Heiser's thresholds of 0.1% and 1% for standard deviation come from OS microbenchmarking. I found no equivalent threshold for whole-genome alignment runs.
- The ACM v1.1 policy page itself returned HTTP 403. The badge definitions come from ACM's terminology-change page.

---

## Q4. Determinism and reproducibility of aligner output: thread count, batching, floating point, GPU vs CPU

### Takeaway
Widely used aligners are "deterministic" only under conditions. BWA-MEM's output can depend on thread count and batching unless the batch size is fixed with `-K`, and on input order, because repeat placement is hashed on read index. Large consortia (CCDG) therefore mandate `-K 100000000`. Claims that GPU and CPU outputs are equivalent are mostly vendor statements, not published read-by-read audits. An exact, certified aligner can make, and should demonstrate, the stronger claim: byte-identical output regardless of threads, batching, hardware and adaptive settings.

### Cited Findings
- Heng Li's explanation, as relayed on Biostars, gives two sources of variation in BWA-MEM:
  - the insert-size estimate, which is computed per batch of reads
  - placement of repetitive hits, which is seeded by a hash of deterministic information, so the i-th input read always gets the same hash and reordering the input changes the results

  [Biostars: re-ordered FASTQ input](https://biostars.org/p/238637); [Biostars: different threads](https://www.biostars.org/p/249286)
- Users have reported different SAM output with `-t 6` versus `-t 16`, which they attributed to the per-batch insert-size estimate. Another user reported different output on two clusters at the same thread count, with no confirmed explanation. — [Biostars p/248777](https://www.biostars.org/p/248777); [Biostars p/9582822](https://biostars.org/p/9582822)
- The CCDG pipeline standard requires `bwa mem -K 100000000` "to achieve deterministic alignment results" and notes that this is a hidden option. — [Biostars p/350890 quoting CCDG PipelineStandard](https://www.biostars.org/p/350890)
- Regier et al. 2018 (*Nature Communications* 9:4038) defined "functional equivalence" standards for WGS pipelines across five genome centres. The resulting variant calls varied much less than sequencing replicates do. — [Regier et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6168605)
- BWA-MEM2 (Vasimuddin et al., IPDPS 2019, doi 10.1109/IPDPS.2019.00041) is described as producing alignments identical to BWA-MEM while being 1.3–3.1× faster. — [NIH Biowulf](https://hpc.nih.gov/apps/bwa-mem2.html)
  - An early DNAnexus evaluation found discrepancies, which went away only with a newer binary from the developers. The identity therefore depends on the version. — [DNAnexus 2020](https://blog.dnanexus.com/2020-03-10-bwa-mem2-review)
- NVIDIA says Parabricks' BWA-MEM gives the same results as the CPU version at up to 60× the speed. No concordance methodology is published there. — [NVIDIA Parabricks fq2bam docs](https://docs.nvidia.com/clara/parabricks/latest/Documentation/ToolDocs/man_fq2bam.html)
  - A third-party Galaxy report found only summary-level agreement (99.6% mapped, mean MAPQ 31.4, matching BWA-MEM2), and says fq2bam wraps BWA-MEM2 where NVIDIA's docs say BWA-MEM. — [Galaxy Hub 2025](https://galaxyproject.org/news/2025-12-15-fq2bam)
- Hoefler Rule 5: state whether the values are deterministic. — [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)
- Weber: report random seeds. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)

### Inferences
- For CERTA, determinism is both a selling point and a testable claim. Run a determinism matrix:
  - thread counts {1, 2, …, max}
  - batch sizes
  - input order (shuffled FASTQ, then sorted by read name for comparison)
  - ARM64 vs x86-64
  - GPU vs CPU-only path
  - each adaptive-parameter trajectory

  Report SHA-256 checksums of the SAM/BAM body (excluding `@PG` and timestamps), or a read-by-read diff of position, CIGAR, flags, MAPQ, NM and AS.
- **Ties.** State the deterministic tie-breaking rule among co-optimal hits, for example lowest reference coordinate first. If CERTA reports all co-optimal hits, state the order and whether any are capped. An exact aligner that breaks ties by thread completion order would lose the claim.
- **Adaptive parameters must not create nondeterminism.** Suppose run-time adaptation reads timing, queue depth or batch-level statistics, which is exactly the BWA insert-size mechanism. Then either:
  - (a) adaptation affects only *how* the certified optimum is found and never *what* is reported, and output identity across adaptation states is tested; or
  - (b) adaptation inputs are deterministic functions of the data.
- **Paired-end.** If CERTA estimates the insert size, estimate it once, globally or from a fixed prefix of the reads, not per batch. Otherwise it inherits BWA's dependence on threads and batching.
- **Floating point.** If scores and checks are integer, say so; that removes the floating-point question. If any GPU reduction uses floating point, the order of operations can change results, and should be pinned or avoided.
- **Competitors.** Document their determinism settings (`bwa-mem2 -K`, and minimap2 and strobealign thread behaviour), so speed is compared at equally reproducible settings.

**CERTA checklist for Q4**
- [ ] A determinism matrix (threads × batch × input order × host × GPU/CPU × adaptive state) with identity results, and the scripts in the artifact.
- [ ] A documented tie-break rule and co-optimal hit reporting policy.
- [ ] Insert-size estimation independent of batching, or tested for it.
- [ ] A statement that scoring and checking are integer-only, or a description of how FP order is fixed.
- [ ] Competitors' determinism flags listed (for example `-K` for BWA-MEM2).

### Gaps
- I found no primary documentation (man page or source comment) of how BWA's `-K` works. That it fixes batch size so batches no longer depend on thread count is widely understood, but I did not confirm it from a primary source here.
- I found no published read-level GPU-vs-CPU concordance study for Parabricks, and no determinism documentation for minimap2, strobealign, Bowtie 2 or Yara.
- I could not check the exact wording of the BWA-MEM2 IPDPS paper on output identity.

---

## Q5. Avoiding overfitting when parameters are tuned or adapted to data: calibration vs held-out evaluation, default vs tuned reporting, preregistration

### Takeaway
Calibrating and validating on the same data is listed as a benchmarking crime. Calibration and evaluation sets must be disjoint, the tuning procedure must be disclosed, and both untuned and tuned results should be reported. Leakage audits from ML-based science (Kapoor & Narayanan) give a ready checklist.

For CERTA the key point is structural. If certification guarantees identical output whatever the parameter values, adaptive parameters cannot overfit *accuracy*, only *speed*. The held-out samples then need to show that the *speed* gains generalise, and that output is byte-identical.

### Cited Findings
- Heiser lists "same dataset for calibration and validation" as a crime: fitting a model and then "validating" it on the same data proves only fit, not predictive power, so the two workloads must be disjoint. He also lists "evaluating a model against itself under the same assumptions". — [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)
- Weber: tuning a new method more extensively than competitors introduces bias, and the tuning strategy must be reported transparently. — [Weber et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584985/)
- Buchka: paying more attention to one's own method's parameters is a documented source of over-optimism. Introducing-paper win rates were 94.7% against 64.3% in neutral studies. — [Buchka et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8111726/)
- Nießl et al.: the many choices in datasets, measures, missing-value handling and aggregation leave room for over-optimistic conclusions. They show how much rankings vary across all combinations of these choices. — [Nießl et al. (arXiv 2106.02447)](https://arxiv.org/abs/2106.02447)
- Kapoor & Narayanan 2023 (*Patterns*, doi 10.1016/j.patter.2023.100804) found data leakage affecting 294 papers in 17 fields (329 in the preprint). They give an eight-type leakage taxonomy and propose "model info sheets" for checking each type. — [Kapoor & Narayanan 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499856/)
- The size of tuning effects in aligners:
  - Rabema: different settings of the same mapper changed sensitivity by more than 20 percentage points. — [Holtgrewe et al. 2011](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3128034/fullTextXML)
  - Mangul: tuning RNA-seq aligners added about 10% correctly mapped reads. — [Mangul et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6437167/)
- Precedents for held-out evaluation in genomics:
  - DeepVariant withholds a whole sample (HG002 in v0.5) and chr20 from training. — [Google Research blog](https://www.research.google/blog/deepvariant-accuracy-improvements-for-genetic-datatypes/)
  - Truth Challenge V2 kept the HG003 and HG004 truth blind while HG002 truth was public. — [NIST/precisionFDA](https://www.nist.gov/node/1838346)
- Hoefler Rule 2 corollary: report all results, not just the best. Rule 11: show upper bounds. — [Hoefler & Belli 2015](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf)
- Strobealign reported its default results and one non-default setting (`-M 40`) separately, rather than replacing the defaults. — [Sahlin 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264/)

### Inferences
- **Protocol to adopt (a lightweight preregistration):**
  1. Write down the adaptive-controller design, its meta-parameters, the search space and the objective (for example median wall time on HG002) before running anything held out.
  2. Tune on HG002 only, ideally on a defined subset such as a 2M-pair subsample or certain chromosomes, and keep the rest of HG002 as an internal validation split.
  3. Freeze the result as a git tag or commit hash, and record the hash and date in the paper, optionally also in an OSF registration or a signed tag.
  4. Only then run HG001 and HG005 on both hosts.
  5. Report every held-out run, including those where adaptation did not help.
- **What counts as "tuned".** A controller that adapts online to each input is not overfitting in itself. Its meta-parameters (thresholds, step sizes, warm-up length, candidate sets) were chosen on HG002, though, and must be treated as tuned hyperparameters subject to the hold-out rule. The same applies to the fixed constants chosen during development.
- **Four modes to report**, for each dataset and host:
  - (a) CERTA-static (shipped defaults, no adaptation)
  - (b) CERTA-adaptive (frozen controller)
  - (c) CERTA-oracle: the best static setting chosen *on the test sample*. This is an upper bound in the sense of Rule 11, shows how much room is left, and is clearly labelled as not a fair comparison.
  - (d) CERTA-calibrated-static: the best static setting from HG002, applied to the held-out samples

  The generalisation gap, (b) or (d) on HG001/HG005 versus on HG002, is the key number reviewers will look for.
- **Cross-host hold-out.** If the controller was calibrated on one server (for example x86-64), report the other (ARM64 + A100) as an additional held-out environment. Hardware is a second axis that speed parameters can overfit to.
- **Cross-data hold-out beyond GIAB.** If feasible, add one different read length, a different instrument or chemistry, or a non-human genome. The GIAB samples share library preparation and instruments across the consortium, so "held-out sample" is weaker than "held-out data-generating process" (an inference; not verified here for the specific HG001/HG005 runs used).
- **Sensitivity analysis.** Show speed as a function of each key parameter around the chosen value on both calibration and held-out data. A flat optimum argues against overfitting.
- **Exactness as a firewall.** State, and test (see Q4), that output is identical for every parameter value the controller can choose. Accuracy is then invariant by construction, and any overfitting can affect only reported speed. This removes most of Buchka's "selective tuning" concern on the accuracy axis.
- **Equal effort for competitors.** Where CERTA reports an adaptive or tuned mode, give each competitor its documented fast or sensitive presets too, or justify why not.

**CERTA checklist for Q5**
- [ ] An evaluation protocol written and timestamped before the held-out runs (commit hash or registration).
- [ ] Calibration data (HG002 subset) and evaluation data (HG001, HG005, the rest of HG002) explicitly disjoint, and both shown.
- [ ] Static, adaptive, calibrated-static and oracle results, with the generalisation gap quantified.
- [ ] Hold-out across hosts (ARM64 + A100 vs x86-64).
- [ ] Parameter-sensitivity plots.
- [ ] Output-identity test across all adaptive settings (accuracy cannot overfit).
- [ ] Equivalent presets for competitors, or a justification.
- [ ] A leakage self-audit, in the spirit of Kapoor & Narayanan's model info sheets: did any HG001 or HG005 run influence a design decision?

### Gaps
- I found no aligner-specific literature on validating *runtime-adaptive* parameters against held-out data. The recommendations are adapted from benchmarking, systems and ML-leakage guidance.
- I found no bioinformatics source explicitly recommending preregistration of benchmark protocols. Buchka et al. do not mention it, and the Nießl et al. summary did not confirm it.

---

## Q6. How papers with certified or verified algorithmic claims validate them: independent checkers, mechanized proofs, testing

### Takeaway
The established pattern comes from certifying algorithms (McConnell et al. 2011) and from proof checking for SAT solvers (DRAT/LRAT, verified checkers such as cake_lpr). It has four parts:
- The program emits a witness.
- A separate, much simpler checker validates (input, output, witness) against the *original problem specification*.
- Trust is reduced to the checker, ideally a formally verified one.
- The paper reports the checker's trusted base and overhead.

Mechanized verification of the alignment dynamic programs themselves exists (Isabelle/HOL), but it is early-stage and mostly at the level of algorithms, not production aligners.

### Cited Findings
- McConnell, Mehlhorn, Näher & Schweitzer 2011, "Certifying algorithms", *Computer Science Review* 5(2):119–161:
  - A certifying algorithm produces, with each output, a certificate or witness (an easy-to-verify proof) that the output was not compromised by a bug.
  - A checker C accepts (x, y, w) exactly when w is a valid witness that y = f(x).
  - The authors argue that for complex tasks only certifying algorithms are satisfactory, and prove that the concept is universal.
  - Certification is the design principle of the LEDA library.

  [ANU research portal](https://researchportalplus.anu.edu.au/en/publications/certifying-algorithms/); [Dagstuhl abstract](https://www.dagstuhl.de/fileadmin/redaktion/Programm/Event/lehrer/AbstractCertifyingAlgorithmsApril2014.pdf)
- Alkassar, Böhme, Mehlhorn & Rizkallah (*J. Automated Reasoning* 2014; preliminary version CAV 2011) give a framework for *formally verifying* the checkers of certifying computations, which shrinks the trusted base to the verified checker. — [ANU research portal (related work listing)](https://researchportalplus.anu.edu.au/en/publications/certifying-algorithms/)
- Beyer's group, working on verification witnesses, stresses a distinction: certificate checkers must check the *original problem*, whereas proof checkers only check that a proof is internally consistent. — [Beyer, Dagstuhl 2025 slides](https://www.sosy-lab.org/research/prs/2025-06-03_Dagstuhl-CertAlg_Witnesses-2.0_Dirk.pdf)
- SAT solving is the most mature large-scale precedent:
  - SAT competitions require the DRAT clausal proof format, and DRAT-trim checks that a proof certifies a formula unsatisfiable. — [Wetzler/Heule et al., "The DRAT format and DRAT-trim checker" (arXiv 1610.06229)](https://arxiv.org/pdf/1610.06229)
  - Efficient *certified* checking of RAT proofs followed. — [Cruz-Filipe et al., "Efficient Certified RAT Verification" (arXiv 1612.02353)](https://arxiv.org/pdf/1612.02353)
  - cake_lpr (Tan, Heule & Myreen 2021) is a checker verified in CakeML down to machine code. That removes unverified extraction and compilation tools from the trusted base, without significant loss of performance. — [Tan, Heule & Myreen (Chalmers full text)](https://research.chalmers.se/publication/531575/file/531575_Fulltext.pdf)
- Mechanized verification of alignment dynamic programs:
  - A 2025 Isabelle/HOL framework (seqAlign) verifies dynamic-programming alignment algorithms, including three-sequence alignment, and generates C++ code from them. — [J. Computer Research and Development](https://crad.ict.ac.cn/en/article/doi/10.7544/issn1000-1239.202330698)
  - A 2025 paper verifies an *imperative* edit-distance algorithm in Isabelle/HOL using generated verification conditions. It notes that earlier Isabelle work covered only functional edit distance. — [Wuhan Univ. J. Nat. Sci. 2025](https://wujns.edpsciences.org/10.1051/wujns/2025306576)
- Heiser names the ratio of proof lines to verified code lines as an example of a meaningless metric. Verification effort should not be offered as a quality measure. — [Heiser](https://gernot-heiser.org/benchmarking-crimes.html)

### Inferences
- **What CERTA must state precisely.** The certified property, for example:
  - every reported hit has optimal score under scoring scheme S, or
  - all hits with edit distance ≤ k are reported, or
  - no better alignment exists within the search space.

  Also state what lies *outside* the certificate: read trimming, paired-end pairing and insert-size logic, MAPQ, soft-clipping, and SAM formatting. Following Beyer, the checker should check against the original problem (read, reference, scoring), not just that the witness is internally consistent.
- **Independent checker.** Use a separate code path with no shared seeding or filtering code. Report its size (lines of code), its trusted dependencies, and whether it is mechanically verified, as in the DRAT and cake_lpr precedent.
  - If the A100 path produces results and a CPU checker validates them, that is also the cleanest GPU-vs-CPU equivalence argument (Q4).
  - Report checker overhead separately in the timings: certified-mode vs uncertified-mode wall time.
- **Testing the checker (fault injection / mutation-style).** Show that the checker rejects corrupted outputs:
  - off-by-one positions
  - altered CIGARs
  - a suboptimal alignment replaced for an optimal one
  - a dropped co-optimal hit
  - wrong strand
  - scores inconsistent with the CIGAR

  Also mutate the aligner itself (disable a filter, introduce a seed bug) and show that the checker catches the resulting wrong outputs. Report detection rates.
- **Cross-validation against independent exact tools.** On simulated and real subsets, compare CERTA's certified hit sets with RazerS 3 or Yara full-sensitivity output and with brute-force dynamic programming on a random sample, Rabema-style (Q2). Report any disagreements, and how each was resolved: a gold-standard flaw, as Rabema allows, or a CERTA bug.
- **Report rejection rates on real data**, which should be zero, and the fallback behaviour if a certificate fails.
- **Avoid meaningless verification metrics** (Heiser). Report the trusted base and test outcomes, not proof-to-code line ratios.

**CERTA checklist for Q6**
- [ ] A formal statement of the certified property and its exclusions.
- [ ] The checker's design, trusted base (LOC, dependencies), verification status, and runtime overhead.
- [ ] Fault-injection or mutation results showing that the checker rejects corrupted outputs and catches bugs in the aligner.
- [ ] Hit sets agree with independent exact tools (RazerS 3, Yara, brute-force DP on a sample), with every disagreement explained.
- [ ] Certificate rejection counts on HG001, HG002 and HG005 and on the simulated data.
- [ ] The checker and witness format released in the artifact, so third parties can re-check the published outputs.

### Gaps
- I found no published short-read aligner that emits per-read optimality certificates checked by an independent checker. The precedents are from general algorithmics and SAT solving.
- I did not retrieve sources on mutation testing (for example DeMillo et al. 1978; Jia & Harman 2011, *IEEE TSE*). The recommendation to fault-inject the checker is an inference and should be backed with citations before it is used.
- The Isabelle alignment-verification papers (2025) were seen only as abstracts and snippets. Their scope (which scoring schemes, whether code is extracted to production) is unverified.
