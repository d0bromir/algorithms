# Evaluating short-read placement accuracy and MAPQ: truth sets, correctness criteria, matched sensitivity, pitfalls

Scope: human short reads, 2012–2026 (Rabema, from 2011, is included because later tools build on it). Our setup, used to judge relevance: HG002 2x150 NovaSeq X real reads, labelled by certifying each read unique on one HG002 Q100 haplotype and projecting to GRCh38 through GIAB's dipcall alignments; plus reads simulated from the HG002 haplotypes. We score at the locus (150 bp) and at matched sensitivity (each aligner's N most confident mates, with ties at the boundary MAPQ split pro rata). We plan HG005 with HPRC release 2 assemblies as a held-out sample.

## Q1. How benchmark papers define a correct placement (window, overlap, CIGAR), and how they handle soft clipping, tandem-repeat shifts, PARs and ALTs

### Takeaway
There is no single standard. The criteria in use, from strict to lenient:
- An endpoint within a window: 20 bp in lh3's wgsim_eval.pl, 30 nt in Qtip, 50% of read length in Teaser.
- Interval overlap of at least 10%: minimap2, paftools mapeval and minibwa (2026). This is the most lenient and the most common in recent lh3 tools.
- An edit-distance-equivalence definition in Rabema.

Endpoint criteria need explicit corrections for clipping, and wgsim_eval tries two clip adjustments. Overlap criteria avoid the clipping problem. Our 150 bp "locus" criterion sits between Teaser and the overlap convention. PARs and ALTs are handled through the choice of reference (the no-alt analysis set, PARs hard-masked on chrY), not in the scoring itself.

### Cited Findings
**lh3 wgsim_eval.pl**
- The default window is 20 bp (`my %opts = (g=>20)`; help text: `-g INT correct if withint INT of the true coordinate`). The true chromosome, leftmost and rightmost positions are parsed from the read name (`/^(\S+)_(\d+)_(\d+)_/`). An alignment on a different chromosome is incorrect, and unmapped reads are skipped. — [wgsim_eval.pl (lh3/wgsim)](https://raw.githubusercontent.com/lh3/wgsim/master/wgsim_eval.pl)
- Strand decides which end is compared. Forward-strand reads are compared on the leftmost coordinate (POS, minus a leading S/H clip). Reverse-strand reads are compared on the rightmost coordinate (POS + M/D/N lengths − 1, plus a trailing S/H clip). — [wgsim_eval.pl](https://raw.githubusercontent.com/lh3/wgsim/master/wgsim_eval.pl)
- The script also computes alternative coordinates (`$left0 -= trailing clip`, `$rght0 += leading clip`) and calls a read correct if either version is within the window. The source comment says this handles "indels that are close to the end of a reads". — [wgsim_eval.pl](https://raw.githubusercontent.com/lh3/wgsim/master/wgsim_eval.pl)

**paftools mapeval and minimap2**
- paftools mapeval computes `r = o/(l1+o+l2)`, where o is the overlap and l1, l2 are the non-overlapping parts; this is a Jaccard index of the two intervals. "The reported mapping is considered correct if r>0.1 by default." — [minimap2 misc/README](https://raw.githubusercontent.com/lh3/minimap2/master/misc/README.md)
- The minimap2 paper (Fig. 1 legend) says a read is correct "if its longest alignment overlaps with the true interval, and the overlap length is ≥10% of the true interval length". The short-read test used 10 million pairs of 150 bp reads simulated with mason2 (`--illumina-prob-mismatch-scale 2.5`) and mapped to the GRCh38 primary assembly. — [Li 2018, minimap2 (arXiv:1708.01492)](https://arxiv.org/abs/1708.01492)
- The paper's wording (≥10% of the true interval) differs slightly from the mapeval README's Jaccard ratio > 0.1.

**minibwa (2026)**
- Li & Homer (arXiv, 13 June 2026) use "|Gs∩Ga|/|Gs∪Ga|≥10%" (simulated interval vs. aligned interval). They simulated 2x150 reads with mason v2.0.9 (`--illumina-read-length 150 --illumina-prob-mismatch-scale 2.5`) "from GRCh38, T2T-CHM13 or the HG002 diploid genomes". One panel excludes centromeres and acrocentric short arms. — [minibwa, arXiv:2606.15357 (HTML)](https://arxiv.org/html/2606.15357v1)

**Qtip**
- An alignment is correct if "the leftmost base involved in the alignment is within 30 nucleotides (nt)" of the true leftmost base, "with appropriate adjustments for soft clipping". The paper does not state a strand requirement. — [Langmead 2017, Qtip, Genome Biology (PMC5557537)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)

**Teaser**
- A mapping is correct if the reported position is within a set distance of the true start (default 50% of read length) and on the simulated strand. Reads below the MAPQ cutoff (default 0) count as unmapped. — [Smolka et al. 2015, Teaser, Genome Biology (PMC4618857)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4618857); [Teaser PDF](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-015-0803-1)

**Rabema**
- Rabema defines correctness by match sets rather than a coordinate window. A match is "feasible" if the read aligns with distance ≤ k; a "best" match has the smallest distance among feasible ones. Mappers are scored in three categories: all, all-best and any-best. The default distance is edit distance (Hamming is optional), and there is an oracle mode for simulated reads. — [Holtgrewe et al. 2011, Rabema (PMC3128034)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3128034); [rabema_evaluate man page](https://manpages.ubuntu.com/manpages/lunar/man1/rabema_evaluate.1.html)
- The authors argue that just counting mapped reads "only allows for a crude assessment of read mappers relative to each other but not to the best possible solution". — [Rabema (PMC3128034)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3128034)

**strobealign**
- The paper evaluates "by looking at the overlap of read alignments with the correct genomic location". The exact criterion is in Additional file 1, Note B, which I did not retrieve.
- Simulated sets: SIM1–SIM4 from hg38 at increasing divergence, and REPEATS (500 copies of a 100 kb sequence at about 90% identity). Read lengths of 50–500 nt were simulated, with analyses focused on 100–500 nt.
- Real reads: HG004 (BIO150 at about 32x, BIO250 at about 26x). — [Sahlin 2022, strobealign, Genome Biology (PMC9753264)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264)

**vg and Giraffe**
- vg (Garrison et al. 2018) plots ROC curves parameterised by MAPQ for 10 million read pairs simulated from the NA24385 (HG002) genome. Results are split into all reads, reads from segments that match the linear reference, and reads from segments that differ from it. — [Garrison et al. 2018 (PMC6126949)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6126949)
- Giraffe reports accuracy as recall vs. false discovery rate, stratified by MAPQ (per a search summary; the exact distance criterion was not retrieved). — [Sirén et al. 2021, Science, doi:10.1126/science.abg8871](https://zendy.io/title/10.1126/science.abg8871)

**PARs and ALTs**
- ALT contigs carry flanks nearly identical to the primary assembly, so reads there get MAPQ 0. lh3 recommends avoiding ALTs (`GCA_000001405.15_GRCh38_no_alt_analysis_set`). PARs duplicated on chrY make PAR variants uncallable: "The right solution is to hard mask PARs on chrY and those extra copies of alpha repeats." — [lh3 blog, "Which human reference genome to use?" (2017)](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
- minibwa does not support alternate contigs ("The only BWA-MEM feature missing from minibwa is the support of alternate contigs"). — [minibwa (arXiv HTML)](https://arxiv.org/html/2606.15357v1)
- dipcall "treats PARs the same way as autosomes". Male samples need `-x dipcall.kit/hs38.PAR.bed`. Outside PARs, chrX regions covered by father contigs and chrY regions covered by mother contigs are filtered out. "users should hard mask PARs on the reference chrY." — [lh3/dipcall README](https://github.com/lh3/dipcall)
- GIAB's Q100 draft README says non-PAR X/Y genotypes `.|1`/`1|.` are written as `1`, and hap.py should be run with `--gender male` for HG002. — [GIAB defrabb v0.019 README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113_README.md)

### Inferences
**How our 150 bp locus criterion compares with the field**
- Our criterion is in line with the field's lenient conventions. For two 150 bp intervals, Jaccard ≥ 0.1 means an overlap o with o/(300−o) ≥ 0.1, so o ≥ ~27 bp. That allows a start offset of up to ~123 bp, more lenient than Teaser's 75 bp and close to our 150 bp locus.
- Most published short-read accuracy numbers from lh3's group (minimap2, minibwa) are therefore "locus-level", not base-exact. Our earlier 20 bp leftmost-base window was the strict wgsim_eval convention.
- The soft-clip failure we observed (a read clipped at an insertion counted as misplaced) is exactly the case wgsim_eval handles. It tests both clip-adjusted endpoints, compares reverse-strand reads on their right end, and ignores I operations in the reference span.

**Recommendations**
- Report a strict and a lenient criterion side by side: an endpoint within 20 bp that is strand- and clip-aware in the wgsim_eval style, and an interval Jaccard ≥ 0.1. This shows that conclusions are not artefacts of the window.
- For tandem-repeat shifts, a Rabema-like equivalence class is the principled fix: any placement whose alignment to the true haplotype has the same (or best) edit distance counts as correct. A cheaper version is to accept any placement inside the same tandem-repeat interval from the GIAB stratifications.
- Count a read as correct when it is placed on chrY in a PAR whose truth is chrX, or vice versa, unless chrY PARs are masked in the reference. Otherwise PAR reads become spurious "misplaced" reads.
- Use the no-alt analysis set for every aligner, so that ALT-aware MAPQ behaviour does not confound the comparison.

### Gaps
- I did not retrieve the exact correctness criteria of GCAT (Highnam et al. 2015), dwgsim_eval (DWGSIM's convention), Giraffe/vg (I could not confirm whether `vg gamcompare` uses a 100 bp radius), or strobealign (its Additional file 1 Note B).
- I did not retrieve the accuracy part of Alser et al. 2021 ("Technology dictates algorithms", Genome Biology 22:249). The sources confirm it surveys 107 aligners and benchmarks 11 ([PSU record](https://pure.psu.edu/en/publications/technology-dictates-algorithms-recent-developments-in-read-alignm/)), but not its correctness definition.
- I did not verify how Rabema represents tandem-shift equivalence (its "interval" gold standard); only the feasible/best match-set definitions were confirmed.

## Q2. Matched-sensitivity / ROC-style comparisons: the curve, how ties are handled, best practice

### Takeaway
The field standard is a cumulative curve. Reads are sorted by MAPQ (descending), and at each MAPQ threshold the curve plots the number (or fraction of all input reads) at or above the threshold against the cumulative error rate. Points exist only at distinct MAPQ values (or MAPQ bins of 10 in wgsim_eval), so an aligner's operating points are discrete. Qtip handles ties explicitly by averaging the penalty within a tie group. That is equivalent to our pro-rata split, or linear interpolation between ROC vertices, so our matched-sensitivity scheme is a defensible formalisation of the standard curve.

### Cited Findings
**Curve definitions**
- wgsim_eval bins MAPQ in tens (`int($t[4]/10)`) and prints, from the highest bin down, the wrong/total counts per bin, the cumulative total and the cumulative error rate (cumulative wrong / cumulative total). With `-a` it iterates over exact MAPQ values instead. — [wgsim_eval.pl](https://raw.githubusercontent.com/lh3/wgsim/master/wgsim_eval.pl)
- paftools mapeval Q-lines report the MAPQ threshold, "the number of reads mapped with mapping quality equal to or greater than the threshold", the "number of wrong mappings", the "accumulative mapping error rate" and "the accumulative number of mapped reads". A U-line reports unmapped reads. — [minimap2 misc/README](https://raw.githubusercontent.com/lh3/minimap2/master/misc/README.md)
- minimap2 Fig. 1: "Read alignments are sorted by mapping quality in the descending order. For each mapping quality threshold, the fraction of alignments (out of the number of input reads) with mapping quality above the threshold and their error rate are plotted along the curve." The error-rate axis is logarithmic (10^-7 to 10^-1 for short reads). — [Li 2018, minimap2 (arXiv:1708.01492)](https://arxiv.org/abs/1708.01492)
- minibwa: "Each point in the figures gives the fraction of mapped reads and the fraction of correctly mapped" reads, and "A mapper of higher accuracy tends to yield a curve closer to the upper left corner." — [minibwa (arXiv HTML)](https://arxiv.org/html/2606.15357v1)
- vg uses "ROC curves parameterised by mapping quality". — [Garrison et al. 2018 (PMC6126949)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6126949)

**Qtip: rank metrics and tie handling**
- Qtip introduces rank-based MAPQ metrics. CID (cumulative incorrect difference) sorts alignments by MAPQ, computes the cumulative count of incorrect alignments, and plots the difference between two MAPQ predictors. CSED does the same with the cumulative squared error Σ(correct − p̂)². The x-axis is "the cumulative number of reads/ends passing the threshold". — [Qtip (PMC5557537)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
- Tie handling in Qtip: "For a group of alignments sharing the same Q, the penalty is averaged across the group's elements in C and C′". For CSED, "the corresponding elements of E and E′ equal the mean squared error of the group". — [Qtip, Endnotes 1–2 (PMC5557537)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557537)
- Qtip defines MAPQ as q = −10·log10(1 − p), where p is the probability that the reported alignment is correct. It also notes that a high alignment score does not imply a high MAPQ. — [Qtip search summary / Genome Biology 2017](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-017-1290-3)

**Parameter sweeps**
- Teaser compares mappers across parameter ensembles at a MAPQ cutoff (default 0) and reports mapping statistics. — [Teaser (PMC4618857)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4618857)

### Inferences
**The standard curve already uses a matched x-axis**
- The standard (lh3) curve already compares aligners at matched sensitivity implicitly, because the x-axis is "fraction of input reads mapped at ≥ threshold". Comparing each aligner at its own MAPQ ≥ 40 (our earlier error) is a single point on each curve, at different x positions. The literature compares whole curves, not fixed MAPQ cut-offs, for exactly this reason.

**Our pro-rata tie split**
- Splitting ties at the boundary MAPQ pro rata equals linear interpolation along the segment between two adjacent ROC vertices. It is also the "expected" error count if tied reads were ordered at random, the same principle as Qtip's within-group averaging. Report it as "interpolated at matched N".
- When the boundary tie group is large (many MAPQ 0 or MAPQ 60 reads), say so. Interpolation hides the fact that the aligner cannot actually be operated at that N.

**Recommendations**
- Plot the full curve on a log error axis, with the denominator = all input reads (unmapped included, as minimap2 does).
- Report the error count at several matched N values (e.g. N at the smallest aligner's MAPQ ≥ 1, ≥ 20 and ≥ 60 sensitivity). Optionally add a calibration-style metric (Qtip CSED / expected vs observed errors per MAPQ bin) to separate "ranking quality" from "calibration".
- Pair-level: mates are usually scored independently, as in minimap2's short-read evaluation. If we score "mates", say whether a pair is counted once or as two ends (Qtip's x-axis is "reads/ends").

### Gaps
- I found no published consensus document on tie handling beyond Qtip's endnotes. wgsim_eval's default 10-wide MAPQ binning coarsens ties further, and I found no discussion of interpolation in lh3's tools.
- I did not find a standard confidence-interval method for differences between ROC curves of aligners. Bootstrapping over reads or regions would be our own choice.

## Q3. Truth for real reads: assembly-derived truth sets, labelling real reads, held-out samples, and dipcall-style GRCh38 alignments

### Takeaway
Real-read placement truth is rarely attempted in the aligner literature. Aligner papers (minimap2, minibwa, strobealign) evaluate placement only on simulated reads, and validate real reads through variant calling against GIAB/syndip.

**HG002: the best-supported truth**
- The HG002 Q100 v1.1 assembly (Hansen et al. 2025) is the best basis for a real-read truth.
- GIAB's v5.0q release ships the exact GRCh38 dipcall alignments we need: hap1 and hap2 BAMs of about 1.0 GB each, run with minimap2 `-z200000,10000`.

**HG005 and HG001: what exists**
- For HG005, HPRC release 2 contains trio-phased hifiasm 0.19.9 assemblies.
- GIAB has no assembly-based benchmark for HG005 or HG001; their latest releases are NISTv4.2.1 (2021).
- Platinum Pedigree / CEPH 1463 (NA12878 family) is a second, pedigree-validated truth source.

### Cited Findings
**HG002 Q100 / T2T-HG002**
- The paper is "A complete diploid human genome benchmark for personalized genomics" (Hansen et al., bioRxiv 2025-09-21, doi:10.1101/2025.09.21.677443). A repository record and a news article report it as published in Cell in 2026 (doi:10.1016/j.cell.2026.06.016), but I did not verify the Cell listing directly. — [Hansen et al. 2025 (PMC12458380)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458380); [bioRxiv PDF](https://www.biorxiv.org/content/10.1101/2025.09.21.677443.full.pdf)
- Construction: about 170x HiFi plus 209x ONT ultra-long reads, Verkko, then 38,037 small corrections and 210 patches (14.18 Mb). QV went from Q63.1 (v0.7) to Q68.9 (v1.1), with v1.1 released July 2024.
- "99.35% of the diploid genome length is free of any detectable errors"; 86.3% of the 38.8 Mb of low-confidence sequence is rDNA.
- The benchmark adds 701.4 Mb of autosomal sequence and 216.8 Mb of sex-chromosome sequence absent from prior benchmarks. 99.9% of the genome (excluding rDNA) is in high-confidence genome-benchmark regions, vs. 85.2% for v4.2.1. — [Hansen et al. (PMC12458380)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458380)
- The v1.1 assembly has gaps (scaffolded Ns) in nine of ten rDNA arrays and is otherwise T2T; it was released on 22 July 2024. — [marbl/HG002 GitHub (search summary)](https://github.com/marbl/HG002)

**GIAB v5.0q (HG002, draft): files and release history**
- Release directory: `.../release/AshkenazimTrio/HG002_NA24385_son/v5.0q/` (last modified 2026-04-07). There is also a `latest/` alias (2026-05-27). — [GIAB HG002 release index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/)
- GRCh38 small-variant files: `HG002_GRCh38_v5.0q_smvar.vcf.gz` (40M), `HG002_GRCh38_v5.0q_smvar.benchmark.bed` (669K) and `HG002_GRCh38_v5.0q_smvar.exclusion_stats.txt`. Matching stvar files exist (VCF 46M, BED 225K), and the same set is provided for GRCh37 and CHM13v2.0. — [GIAB v5.0q directory](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/)
- `v5.0q/dipcall_output/` contains, for GRCh38:
  - `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap1.bam` (1.0G, .bai 1.5M)
  - `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap2.bam` (1.0G, .bai 1.5M)
  - `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed` (9.5K)
  - `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz` (40M)
- The same dipcall outputs exist for CHM13v2.0 (hap1 935M, hap2 1.0G) and GRCh37 (1.1G each). Identical files were in the earlier draft directory `.../data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/dipcall_output/`. — [v5.0q dipcall_output](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/dipcall_output/); [defrabb v0.019 dipcall_output](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/dipcall_output/)
- The v5.0q README calls itself the "draft 'v5.0q'" assembly-based benchmark built from T2T HG002 Q100 V1.1 with defrabb v0.020 (data dated 2025-01-17; README updated 2026-07-17). "We plan to publish this work soon"; until then, users are asked to cite doi:10.1101/2025.09.21.677443. — [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)

**GIAB v5.0q: method, exclusions and recommended evaluation**
- Method: dipcall v0.3 with minimap2 `-z` "modified from default to `-z200000,10000`". The README does not explain the "z2k" filename tag, and it does not state whether the GRCh38 used was the no-alt analysis set or whether chrY PARs were masked. — [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)
- Exclusions:
  - assembly gaps plus flanks
  - satellites, tandem repeats >10 kb and segmental duplications with a break in the assembly-to-reference alignment on either haplotype (segdups and satellites within 40 kb are merged first)
  - regions without a 1:1 alignment (except X/Y)
  - SV-affected regions (from smvar)
  - PAV/dipcall discrepancies and inversions
  - dipcall bugs (adjacent insertions/deletions)
  - TSPY2
  - VDJ
  - suspected assembly errors and mosaic/de novo variants

  — [defrabb v0.019 README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113_README.md); [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)
- Recommended evaluation: hap.py with `--engine vcfeval`, `--gender male` (hap.py "does not predict sex correctly"), and GIAB v3.6 stratifications. For SVs, use truvari bench `--refine` after filtering `ALT=*` records. Complex SVs are no longer excluded. — [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)
- v5.0q supersedes v4.2.1 for HG002, which is deprecated. — [search summary citing GIAB pages](https://www.nist.gov/node/1559076)

**dipcall confident regions (dip.bed)**
- A base is included if it is covered by one ≥50 kb alignment with mapQ ≥ 5 from each parent and is not covered by other ≥10 kb alignments. The BED filter is stricter than the VCF filters (HET1, GAP1). — [lh3/dipcall README](https://github.com/lh3/dipcall)

**GIAB stratifications**
- The latest is v3.6 (directory dated 2025-04-25; earlier versions v2.0–v3.5). — [GIAB genome-stratifications index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/)
- The resource paper is Dwarshuis et al., "The GIAB genomic stratifications resource for human reference genomes", Nat Commun 15:9029 (2024), doi:10.1038/s41467-024-53260-y. It covers GRCh37, GRCh38 and CHM13 and highlights the extra hard-to-map and GC-rich regions in CHM13. — [Dwarshuis et al. 2024 (PMC11489684)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11489684)

**CMRG**
- Wagner et al. 2022, Nat Biotechnol 40:672–680, doi:10.1038/s41587-021-01158-1, characterise 273 of 395 challenging autosomal genes using a haplotype-resolved assembly. The benchmark has over 17,000 SNVs, 3,600 indels and 200 SVs on GRCh37/38 for HG002. — [Wagner et al. 2022 (UKY record)](https://scholars.uky.edu/en/publications/curated-variation-benchmarks-for-challenging-medically-relevant-a/)
- The release directory is `CMRG_v1.00/` (2021-07-08). — [GIAB HG002 release index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/)
- A search summary of the CMRG data listing mentions files named like `HG002v11-align2-GRCh38.hap1.bam`. — [search result, data.amerigeoss.org listing](https://data.amerigeoss.org/dataset/challenging-medically-relevant-genes-benchmark-set)

**Platinum Pedigree and the CEPH 1463 assemblies (two separate papers)**
- "The Platinum Pedigree: a long-read benchmark for genetic variants" is in Nature Methods 22(8):1669–1676 (2025), doi:10.1038/s41592-025-02750-y.
  - It covers CEPH 1463, with Mendelian-inheritance filtering of Illumina, HiFi and ONT calls.
  - It reports >4.7 M SNVs, 767,795 indels, 537,486 tandem repeats and 24,315 SVs over 2.77 Gb of GRCh38.
  - Retraining DeepVariant on it cut genotyping errors by about 34%.
  - It is licensed CC BY 4.0 and hosted on the AWS Open Data registry.

  — [Platinum Pedigree on AWS Open Data](https://registry.opendata.aws/platinum-pedigree/); [PacBio press release](https://www.pacb.com/press_releases/nature-methods-paper-leverages-pacbio-sequencing-technology-to-develop-the-platinum-pedigree-benchmark-a-new-standard-for-accurate-characterization-of-variation-in-the-human-genome-that-improves-trai/)
- The separate Nature paper is Porubsky et al. 2025, "Human de novo mutation rates from a four-generation pedigree reference", Nature 643:427–436. It "phased and assembled more than 95% of each diploid human genome" in the 28-member CEPH 1463 family and calls the result a truth set. — [PubMed 40269156](https://pubmed.ncbi.nlm.nih.gov/40269156/)
- So "Platinum Pedigree (Porubsky et al. 2025, Nature)" in our plan conflates two papers: the variant benchmark is in Nature Methods, and the near-T2T family assemblies are in Nature.

**HPRC release 2 and HG005**
- Release 2 was announced on 12 May 2025 (pre-publication) and covers 232 individuals. Data are on AnVIL, AWS S3 and GCP, and in INSDC, under CC0. — [HPRC Data Release 2](https://humanpangenome.org/hprc-data-release-2/)
- The assembly index is `data_tables/assemblies_release2_v1.0.index.csv` (13 columns, including S3 path, md5, fai, gzi). "Reference-level assemblies for CHM13v2, GRCh38, and HG002 (Q100) are included in the index file." — [hprc_intermediate_assembly data_tables](https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables)
- The index has HG005 rows for `HG005_pat_hprc_r2_v1.0.1` (haplotype 1, GCA_018506945.2) and `HG005_mat_hprc_r2_v1.0.1` (haplotype 2, GCA_018506965.2). Both are trio-phased, assembled with hifiasm 0.19.9 and dated 2024-12, at `s3://human-pangenomics/working/HPRC_PLUS/HG005/assemblies/release2/HG005_{pat,mat}_hprc_r2_v1.0.1.fa.gz`. — [assemblies_release2_v1.0.index.csv](https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv)
- I found no HG001 rows in the first ~200k of the ~228k characters of that index.
- A JBrowse README points to a release-2 all-vs-GRCh38 PAF (`hprc465vsgrch38.aln.paf.gz`) and a Minigraph-Cactus multiple alignment projected onto GRCh38. These are not dipcall-style per-sample BAMs. I did not open the file. — [JBrowse HPRC multiway README (search summary)](https://jbrowse.org/demos/hprc_multiway/README.txt)

**GIAB HG005 and HG001 releases**
- GIAB's HG005 release directory has only NISTv3.3, v3.3.2, v4.2.1 (2021-09-29) and `latest/`. — [GIAB HG005 release index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/)
- HG001 likewise has nothing newer than NISTv4.2.1 (2021). — [GIAB HG001 release index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/)

**Labelling real reads with their origin**
- Aligner papers do not label real reads.
  - minibwa evaluates real reads only for speed and through DeepVariant + vcfeval against the Q100 v1.1 truth. Its placement-accuracy curves come from simulations, including from "the HG002 diploid genomes". — [minibwa (arXiv HTML)](https://arxiv.org/html/2606.15357v1)
  - minimap2 validated real reads via variant calling on syndip (ERR1341796, a CHM1+CHM13 mix with an assembly-derived truth). — [minimap2 (arXiv:1708.01492)](https://arxiv.org/abs/1708.01492)
- Mapping to the sample's own assembly is used for sequencing-error benchmarking. McGee et al. (bioRxiv 2026, "Short-Read Sequencing Benchmarking with Donor-Specific Assemblies") benchmark nine chemistries on seven sequencers using HG002 and COLO829BL donor-specific assemblies, attributing mismatches to sequencing error. The preprint is unreviewed, and the vendor blog summarising it has a conflict of interest. — [Sciety: McGee et al. 2026](https://sciety.org/articles/activity/10.64898/2026.06.24.734333); [Element Biosciences blog](https://www.elementbiosciences.com/blog/measuring-what-actually-matters-the-quest-for-unbiased-sequencing-benchmarks)
- Mapping to a diploid self-assembly depresses MAPQ in homozygous regions. In the HPRC draft-pangenome paper, aligning HiFi reads back to each sample's diploid assembly gave higher identity, but about 20% more reads with MAPQ < 10 (homozygous regions), so they did not filter on MAPQ. — [Liao et al. 2023, Nature (draft human pangenome)](https://www.nature.com/articles/s41586-023-05896-x)

### Inferences
**Recommendations for HG002**
- For our real-read truth, use the v5.0q `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap{1,2}.bam`. Restrict scored reads to `HG002_GRCh38_v5.0q_smvar.benchmark.bed` (or at least to the dipcall `dip.bed`). Outside those regions, GIAB itself found the 1:1 projection unreliable.
- Our "unique on one haplotype" certification is the right guard against the homozygous-region ambiguity seen in the HPRC paper. Reads unique on hap1 and hap2 separately, but identical on both, should project to the same GRCh38 locus through the two BAMs. A disagreement between the two projections is a useful flag for projection error.

**Recommendations for HG005 (held-out)**
- HPRC r2 provides HG005 pat/mat hifiasm 0.19.9 trio assemblies, but there is no GIAB-produced dipcall output and no assembly-based GIAB benchmark. We must run dipcall ourselves: v0.3, `-z200000,10000`, `-x hs38.PAR.bed` (HG005 is male), on the no-alt GRCh38 with chrY PARs hard-masked, mirroring GIAB.
- Downstream validation of HG005 can only use NISTv4.2.1, which covers less of the genome and is easier than v5.0q.
- Expect noisier labels than HG002. The hifiasm assemblies are not polished to Q100 and are not T2T, so restrict to dip.bed and report label quality (e.g. the fraction of certified reads in dip.bed).
- HG001/NA12878 via Porubsky et al.'s CEPH 1463 assemblies (and the Platinum Pedigree variant truth) is a possible third sample. HPRC r2 did not show HG001 in the part of the index I read.

### Gaps
- I did not obtain the v5.0q GRCh38 benchmark-region size, variant counts, or whether the GRCh38 used was the no-alt analysis set with masked chrY PARs; the README doesn't state these.
- The meaning of "z2k" in the filenames is not explained (the README says `-z200000,10000`).
- I could not confirm whether HPRC r2 publishes per-sample dipcall-style GRCh38 BAMs (as opposed to the all-vs-GRCh38 PAF), or HG001's status in r2 (the last ~28k characters of the index were not read).
- I did not verify the Platinum Pedigree author list, whether its assemblies (vs. variants) are distributed, or their QV.
- I did not retrieve the CHM13-reads-mapped-to-CHM13 analyses (Aganezov et al. 2022, Science) or k-mer-based read-origin labelling methods; no specific short-read placement-truth paper based on k-mers was found.

## Q4. Pitfalls: projection errors, easy-region label bias, circularity, coordinate conventions, sequence without a 1:1 counterpart

### Takeaway
The main documented risks:
- Assembly-to-GRCh38 projection fails in repeats, SVs and gene-conversion regions.
- Truth regions are biased toward easy sequence.
- Large SVs (extra gene copies) between the sample and the reference cause most high-MAPQ mismapping, so GRCh38-simulated reads overstate accuracy.
- Downstream metrics can be biased by callers tuned to a particular aligner's MAPQ.

### Cited Findings
**Projection between assemblies**
- In the Q100 paper, genome benchmarking found 66,544 errors in HPRCv1 assemblies. Only 50,843 could be lifted to GRCh38 and 30,505 fell in v4.2.1 high-confidence regions; "more than half of the total genome benchmarking errors" were invisible to variant benchmarking.
- Some regions "could not be reliably mapped between HG002 and GRCh38 due to repeats, structural variation, or gene conversion". An example is DUSP22 on maternal chr16, which is absent from GRCh38. — [Hansen et al. (PMC12458380)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458380)
- GIAB: "some regions likely remain where the assembly could be aligned to the reference in a different way"; "some errors remain in highly homozygous regions, homopolymers, and tandem repeats". Regions without a 1:1 alignment, and large repeats with alignment breaks, are excluded. — [defrabb v0.019 README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113_README.md); [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)
- Q100 assembly caveats: k-mer QV can overestimate PacBio assembly quality (by 7.9–11.7 QV) and misses loss of heterozygosity; homopolymers > 20 bp remain hard for all assemblies. — [Hansen et al. (PMC12458380)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458380)

**Sequence without a 1:1 counterpart**
- minibwa: "most false read alignment at high mapping quality is caused by large structural variants between individuals". If the reference has one gene copy and the sample an extra one, reads from the extra copy mismap at high MAPQ. Also, "the theoretical accuracy as is shown in Fig. 2 may not reflect accuracy in downstream data processing". BWA-MEM's simulated advantage comes mostly from centromeres and acrocentric arms and does not translate into a variant-calling difference. — [minibwa (arXiv HTML)](https://arxiv.org/html/2606.15357v1)
- CMRG: false duplications in GRCh37/GRCh38 cause reference-specific missed variants (e.g. CBS, CRYAA, KCNE1). Masking them raises recall in affected genes from 8% to 100%. — [Wagner et al. 2022 (UKY record)](https://scholars.uky.edu/en/publications/curated-variation-benchmarks-for-challenging-medically-relevant-a/)

**Bias toward easy regions**
- In GA4GH's analysis, SNV concordance between two methods was 99.7% inside high-confidence regions vs. 76.5% outside them. — [Krusche et al. 2019 (search summary; doi:10.1038/s41587-019-0054-x)](https://www.nist.gov/publications/best-practices-benchmarking-germline-small-variant-calls-human-genomes)
- The CHM13 stratifications show extra hard-to-map and GC-rich regions that lower benchmarking performance. — [Dwarshuis et al. 2024 (PMC11489684)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11489684)

**Aligner/caller coupling**
- "Variant callers use MAPQ scores, so callers tuned on a particular aligner's MAPQ may favor that aligner." — [strobealign (PMC9753264)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264)

**Simulation pitfalls**
- Reads from segments that differ from the linear reference behave very differently from reads that match it; vg's advantage over BWA-MEM is concentrated in the former. So results should be stratified by whether the read carries non-reference sequence. — [Garrison et al. 2018 (PMC6126949)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6126949)
- wgsim names encode the truth as `chr_left_right_...` with 1-based coordinates compared against SAM POS. Converting between 0-based BED/PAF and 1-based SAM is a classic off-by-one source; wgsim_eval subtracts 1 from the computed right end. — [wgsim_eval.pl](https://raw.githubusercontent.com/lh3/wgsim/master/wgsim_eval.pl)

### Inferences
**Pitfalls specific to our setup**
- Circularity: our truth projection uses minimap2-based dipcall alignments of whole contigs, not read alignments. This limits circularity at the read level, but projection choices in segdups and tandem repeats inherit minimap2's chaining decisions. Aligners whose seeding and chaining resemble minimap2 (minimap2 itself, minibwa, possibly strobealign) could be favoured in exactly the ambiguous regions. Mitigations:
  - restrict to the v5.0q benchmark BED or dip.bed;
  - report results stratified by GIAB v3.6 segdup, low-mappability and tandem-repeat strata;
  - optionally check robustness by re-projecting with a different assembly aligner (e.g. a PAV-style or Liftoff-style projection).
- Label bias toward easy regions: certifying reads as unique on a haplotype drops reads from sequence repeated within the haplotype. That is where aligners differ most, so the real-read truth systematically under-samples the hard cases. Report the certified fraction per stratum, and complement it with HG002-haplotype simulations, which include those regions.
- Reads with no 1:1 GRCh38 counterpart (extra gene copies, insertions absent from GRCh38, DUSP22-like loci) should get a separate truth class: "no correct GRCh38 placement". Count any confident placement of them as an error (the minibwa high-MAPQ failure mode), rather than dropping them. Unmapped or MAPQ-0 output is the correct behaviour for these reads.
- Coordinate conventions: dipcall BAMs give contig→GRCh38 alignments. When projecting a read's haplotype coordinate, use the CIGAR of the covering alignment record, handle reverse-strand contigs, and decide the convention for bases inside contig insertions (no reference counterpart) and for reads spanning a deletion. Document 0- vs 1-based handling explicitly.

### Gaps
- I found no published quantitative estimate of the projection error rate of dipcall-based liftover for short-read coordinates specifically. The Q100 paper reports projection failures only for assembly errors.
- I found no literature that labels real short reads by k-mer uniqueness on a personal assembly for aligner benchmarking. This appears to be relatively novel (our method), so there is no published validation to compare against.

## Q5. Downstream validation: how mapping changes are validated in variant calling, and what effect sizes are meaningful

### Takeaway
Mapping changes are validated by calling variants (DeepVariant, GATK HaplotypeCaller or bcftools) and comparing to GIAB truth with hap.py/vcfeval (GA4GH framework), stratified by genome context. For HG002 the current truth is v5.0q (or Q100 v1.1 directly) with v3.6 stratifications, plus CMRG for difficult medically relevant genes.

Effect sizes reported as "comparable" by aligner authors are about 1–3% of FN counts genome-wide (e.g. minibwa vs BWA-MEM). Differences presented as meaningful are larger or concentrated in difficult strata: a 2.5–2.7% SNP F-score difference for strobealign with bcftools, recall rising from 8% to 100% in falsely duplicated genes, and a ~34% error reduction from DeepVariant retraining.

### Cited Findings
**Framework and recommended settings**
- GA4GH best practices (Krusche et al., Nat Biotechnol 2019, doi:10.1038/s41587-019-0054-x) give guidance on matching calls with different representations, defining standard metrics, and stratifying by variant type and genome context. The tools are on GitHub and precisionFDA. — [NIST publication page](https://www.nist.gov/publications/best-practices-benchmarking-germline-small-variant-calls-human-genomes); [bioRxiv 270157](https://www.biorxiv.org/content/10.1101/270157v2)
- For v5.0q, GIAB recommends hap.py with `--engine vcfeval`, `--gender male`, and the v3.6 stratifications. — [v5.0q README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md)

**Reported effect sizes**

| Comparison | Truth and caller | Reported difference | Authors' reading |
|---|---|---|---|
| minibwa vs BWA-MEM vs strobealign (2026) | HG002 short reads to GRCh38, DeepVariant v1.10.0, RTG vcfeval v3.13 vs GIAB Q100 v1.1 | SNP FN 46,367 / 46,895 / 57,545; SNP FP 7,544 / 7,585 / 8,108; indel FN 36,321 / 37,425 / 37,828; indel FP 5,308 / 5,218 / 5,655 | minibwa has "528 fewer SNP FNs, 41 fewer SNP FPs, 1,104 fewer indel FNs, and 90 more indel FPs" than BWA-MEM and is "comparable" ([minibwa](https://arxiv.org/html/2606.15357v1)) |
| minimap2 vs BWA-MEM (2018) | syndip, GATK HC v3.5 | SNP FNR 2.6% vs 2.3%; SNP FPPM 7.0 vs 8.8; indel FNR 11.2% vs 11.3%; indel FPPM 6.4 vs 6.5 | "broadly comparable" ([minimap2](https://arxiv.org/abs/1708.01492)) |
| strobealign vs BWA-MEM (2022) | HG004 BIO150/BIO250, bcftools | SNP recall −1.0% / −1.9%; precision +5.4% / +5.3%; F-score +2.7% / +2.5%; indels similar (~0.1%) | presented as a gain ([strobealign](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264)) |
| strobealign vs BWA-MEM, simulated | SIM150 | SNP recall 0.45% lower at the same precision; indel recall +3.3% | — ([strobealign](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753264)) |
| DeepVariant retrained on Platinum Pedigree | Platinum Pedigree | about 34% fewer genotyping errors | — ([AWS registry](https://registry.opendata.aws/platinum-pedigree/)) |
| Masking GRCh38 false duplications | CMRG | recall in affected genes from 8% to 100% | — ([Wagner et al. 2022](https://scholars.uky.edu/en/publications/curated-variation-benchmarks-for-challenging-medically-relevant-a/)) |

- Relative to BWA-MEM's totals, the minibwa differences are about 1.1% of SNP FNs and about 2.9% of indel FNs (computed from the minibwa counts above).

**Graph mappers**
- Giraffe: DeepVariant called SNVs and indels more accurately with Giraffe pangenome alignments than with single-reference alignments. The exact numbers were not retrieved. — [Sirén et al. 2021 (Science, doi:10.1126/science.abg8871)](https://zendy.io/title/10.1126/science.abg8871)

### Inferences
**Recommendations for downstream validation**
- For CERTA, validate placement changes downstream with DeepVariant (or GATK, matching what reviewers expect) and hap.py `--engine vcfeval --gender male` against v5.0q GRCh38 smvar plus CMRG. Stratify with GIAB v3.6, at minimum: all benchmark regions, low-mappability, segdups, tandem repeats/homopolymers, and CMRG genes.
- Genome-wide differences of about 1% in FN counts are framed as "comparable" by the field's most prominent aligner author. A claim of meaningful improvement should therefore be either a multi-percent change in FN or FP counts, or a concentrated change in a difficult stratum.
- Each claim should be backed by consistency across HG002 and the held-out sample (HG005 against NISTv4.2.1) and, ideally, across two callers. strobealign's MAPQ-tuning caveat implies that one caller alone is insufficient.
- Report absolute FN/FP counts as well as F1. Small F1 deltas on ~4 M SNVs correspond to hundreds or thousands of variants, and readers of minibwa/minimap2-style papers expect counts.

### Gaps
- I found no published threshold or statistical test for "meaningful" differences in hap.py metrics (e.g. bootstrap CIs over regions). Effect-size norms above come from how authors describe their own results, not a formal standard.
- I did not retrieve the precisionFDA Truth Challenge V2 (Olson et al. 2022) effect sizes or the HPRC 2023 short-read variant-calling improvements as numbers.
