# T2T-CHM13, pangenome references and paralog-resolving methods for short-read placement (segdups, medically relevant genes)

Scope note: human short reads, 2018-2026. For every result the evaluation design (simulated vs real reads, truth source, coordinate system) is recorded where the source gave it. Context from our own HG002 NovaSeq X measurements is used only in the Inferences bullets.

## 1. What T2T-CHM13 measurably gains for short reads (Aganezov 2022, Nurk 2022, Vollger 2022, CMRG masking), and whether it survives projection to GRCh38

### Takeaway
T2T-CHM13 gives a small global gain in mapping (+0.97 % properly paired reads, 20-25 % lower mismatch rate). The large gains are local: in GRCh38 collapsed duplications, short-read precision rises from about 64-67 % to 98-99.7 %, and in GRCh38 false duplications recall rises from 1-2 % to 57-68 %. Mapping to CHM13 and lifting to GRCh38 (levioSAM2) keeps most of the gain: 23.9 % fewer GATK errors and 11.4 % fewer DeepVariant errors genome-wide, 35.2 % and 20.0 % fewer in CMRG genes. The headline "12-fold fewer false positives in 269 medically relevant genes" applies to HiFi, not to Illumina.

### Cited Findings
**Reference content (Nurk 2022; Vollger 2022)**
- T2T-CHM13 is 3.055 Gbp and gapless except chrY. It adds nearly 200 Mbp of sequence carrying 1,956 gene predictions, 99 of them predicted protein-coding. The new sequence includes centromeric satellites, recent segmental duplications (SDs) and all five acrocentric short arms — [Nurk et al. 2022, Science, abstract via UCSC](https://genomics.ucsc.edu/?p=222194)
- SDs make up 7 % of CHM13 (218 Mbp of 3.1 Gbp), and CHM13 adds 51 Mbp of SD sequence. Table 1 gives 207.6 Mbp of SD (6.67 %) in CHM13 v1.0 against 167.3 Mbp (5.37 %) in GRCh38. A further 81.3 Mbp of SD was previously unresolved or structurally variable: 16.5 Mbp is wholly absent from GRCh38 and 64.8 Mbp differs in copy number or structure — [Vollger et al. 2022, Science](https://pmc.ncbi.nlm.nih.gov/articles/8979283/)
- Short-read read-depth copy number for 268 SGDP genomes: CHM13 approximates the population median copy number for 94 % of SD bases (70.6 Mbp), GRCh38 for 57 % (42.8 Mbp). In non-syntenic SD regions, GRCh38 underestimates copy number by a mean of 9.2 copies (median 3.0). CHM13 better represents 119 protein-coding genes, GRCh38 better represents 65. Of 182 candidate previously unresolved genes, 166 (91 %) are in SD gene families — [Vollger 2022](https://pmc.ncbi.nlm.nih.gov/articles/8979283/)
- SNV density is 1.47 SNVs/kbp in SDs against 0.95/kbp in unique sequence. Only 124 Mbp of SDs (60 %) had a clear 1:1 orthology between GRCh38 and CHM13, which limits lift-over in exactly these regions. Only 25 % of SMN1/SMN2 haplotypes were fully resolved by HiFi assembly at that time — [Vollger 2022](https://pmc.ncbi.nlm.nih.gov/articles/8979283/)

**Mapping and variant calling, 3,202 1KGP samples (Aganezov 2022)**
- Design: all 3,202 1KGP samples (≥30x PE Illumina) were remapped with BWA-MEM to CHM13 v1.0, matching the GRCh38 1KGP pipeline (GRCh38 with ALT and decoys). There was no read-level truth; quality was measured by mismatch rate, pair orientation, coverage uniformity, Mendelian concordance and the GIAB HG002 CMRG benchmark — [Aganezov et al. 2022, Science (author manuscript)](https://escholarship.org/uc/item/18c391gh)
- BWA-MEM mapped an additional 7.4 x 10^6 (0.97 %) properly paired reads to CHM13, even against GRCh38 with ALT and decoys. The per-read mismatch rate was 20-25 % lower in every continental population, misoriented pairs decreased, and the standard deviation of coverage within genes fell 4-fold — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- GATK HaplotypeCaller PASS calls: 126,591,489 variants on CHM13 (per-sample mean 4,717,525) against 125,484,020 on GRCh38 (per-sample mean 5,101,897). Fewer calls per sample reflect fewer reference-private alleles and fewer false heterozygous calls from collapsed SDs. The AF ≈ 0.5 deficit is "largely explained by the corrections to collapsed SDs", since a false pileup of reads from duplicated regions creates heterozygous paralogous sequence variants (PSVs) — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- In 602 trios, the median Mendelian-discordant variant count fell from 8,879 (GRCh38) to 7,484 (CHM13), and the median number of low-quality variants from 890,701 to 682,609 (Wilcoxon P ≈ 10^-96) — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- Unique sequence: relative to GRCh38, CHM13 adds 14.9 Mbp of sequence unique at 50-mers, 23.5 Mbp at 100-mers and 39.5 Mbp at 300-mers. More than 106 Mbp absent from GRCh38 needs reads longer than 300 bp to map uniquely. GRCh38 has 28 large identical interchromosomal sequence pairs plus 42 pairs involving unplaced contigs; all but one are non-identical in CHM13 — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)

**Falsely collapsed and falsely duplicated regions in GRCh38 (Aganezov 2022)**
- Design, effectively an assembly truth: CHM13 is haploid, so any heterozygous call from mapping CHM13 reads to GRCh38 marks a collapse, copy-number variation or a mapping error. CHM13 HiFi reads and Illumina-like simulated 150 bp reads from CHM13 were mapped to GRCh38. This gave 368,574 heterozygous SNVs, of which 56,413 (15.3 %) were shared between the two technologies, clustered into 908 putative problematic regions (20.8 Mbp; 73.6 % SD-associated) — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- Using SGDP copy number, 203 loci (8.04 Mbp) have missing copies in GRCh38, affecting 308 gene features; 14 of 48 protein-coding genes are fully inside such a region (e.g., DUSP22, KMT2C). Separately, 12 regions (1.2 Mbp, 74 genes, 22 of them protein-coding) are falsely duplicated in GRCh38. Five of these are BAC clones on chr21q misplaced from the p-arm — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- Short-read accuracy in those regions, measured against HiFi calls for two GIAB/PGP trios:
  - GRCh38 collapsed duplications (copy number <10, ~910 kbp): Illumina recall was similar (20.1-28.3 % on CHM13, 21.5-25.4 % on GRCh38). Precision rose from 64.3-67.3 % on GRCh38 to 98.1-99.7 % on CHM13.
  - False duplications (1.2 Mbp): Illumina recall rose from 1.1-1.8 % (GRCh38) to 57.4-68.3 % (CHM13), and precision from 76.5-95.8 % to 98.5-99.3 % — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- Medically relevant genes: 622 genes with rare or erroneous structural alleles in GRCh38 appear corrected in CHM13. In KCNE1 most reads map to the false duplicate KCNE1B on chr21p, and KCNJ18 sits in a 17p11.2 collapse whose missing paralog is "KCNJ18-1". Seventeen medically relevant genes lie in erroneous duplicated or collapsed GRCh38 regions, with 58 variants/kbp in collapsed and 21/kbp in falsely duplicated genes against 40/kbp for the rest — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- Benchmark of 269 challenging medically relevant genes on HG002 (Illumina-BWA-MEM-GATK, HiFi-PEPPER-DV, ONT-PEPPER-DV): false positives and false negatives "substantially decreased for all three". The 12-fold reduction is for HiFi false positives (missing KMT2C-like sequence added; false duplications of CBS, CRYAA, H19 and KCNE1 removed). True positives dropped by about 14 %, because CHM13 carries fewer private alleles. The Illumina-specific counts are in Table S15, which I did not retrieve — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- On calling on CHM13 and lifting to GRCh38, the authors warn that this "requires special handling" where the two references carry different alleles, that some variants "become invisible", and that representation differences in low-complexity regions make lifted call sets differ. For ClinVar, 99.8 % of variants lift GRCh38 to CHM13, as do 99.6 % of pathogenic/likely pathogenic ones — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)

**Masked GRCh38 (GIAB CMRG; Wagner 2022)**
- The CMRG benchmark covers 273 of 395 challenging genes (HG002, from a haplotype-resolved assembly): 27,178 small variants and 217 SVs on GRCh38. Genes falsely duplicated in GRCh38 include CBS, CRYAA, KCNE1, H19 and U2AF1. Masking the false duplications (sequence set to N, coordinates kept) raised Illumina-BWA-MEM-GATK sensitivity in the affected genes "from 8% to 100%". A v2 masked GRCh38 adds further false duplications found by T2T, and the authors recommend it — [Wagner et al. 2022, Nat Biotechnol](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117392/)
- The hs37d5 decoy removed 1,272 false positives but caused 78 false negatives (52 in CYP4F12, 18 in LMF1). Decoys and masking both trade sensitivity for precision — [Wagner 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117392/)

**Map to CHM13, lift to GRCh38 (levioSAM2; Chen et al. 2024, Nat Methods)**
- Mechanism:
  - Reads are mapped to CHM13 and lifted to GRCh38 through a chain file (nf-LO, minimap2 mode).
  - "Suppress": reads in source-only regions with no target counterpart (CHM13 centromeres, sequence missing from a GRCh38 collapse) are dropped. Source segments over 5 kbp missing from the chain and unmapped or repetitive in the target are labelled unliftable.
  - "Defer": low-confidence lifts are re-mapped directly to GRCh38. For BWA-MEM a lift is committed only at MAPQ ≥30 or AS ≥100 (Bowtie 2: MAPQ ≥10 or AS ≥ −10), and lifts into "mappability-reduced" regions are also deferred.
  - "Reconcile": the version with lower edit distance to the target wins; ties go to higher MAPQ, then random. MAPQ is not recalculated after lifting — [Chen et al. 2024, levioSAM2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/)
- Real Illumina reads (30x HG001/HG002/HG005, GIAB v4.2.1, scored in GRCh38): CHM13-then-GRCh38 avoided 18,308 (23.9 %) GATK errors and 3,443 (11.4 %) DeepVariant errors. Mean GATK false positives fell by 11,551 (24.9 %). GATK errors in the LowMap, SegDups and OtherDifficult strata fell 27.5-42.4 %, and DeepVariant LowMap errors 15.7-22.3 %. In HG002 CMRG, errors fell 35.2 % with GATK (730 errors) and 20.0 % with DeepVariant (1,010) — [levioSAM2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/)
- Simulated reads (mason, 10M 100 bp pairs, chr20 and chr21, correct means leftmost position within 10 bp): only 0.08 % (chr20) and 1.37 % (chr21) more correct mappings than direct GRCh38. No high-MAPQ error rate was reported. For HiFi, the lift was neutral genome-wide (+0.8 % errors) but removed 22.5 % of CMRG errors — [levioSAM2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/)
- The abstract's "11.4-39.5 %" range mixes GRCh38 and GRCh37 targets (39.5 % is GATK to GRCh37). For GRCh38 the range is 11.4 % (DeepVariant) to 23.9 % (GATK) — [levioSAM2 search summary and PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11610747/)

### Inferences
- Our result that BWA-MEM2 on CHM13 places 38 % of the mates no GRCh38 aligner gets right, and that CHM13 as a second opinion halves confident errors, is consistent with and somewhat larger than the literature: levioSAM2 has 24-35 % fewer GATK errors (variant level), and Aganezov has a 3-fold drop in false-positive rate in GRCh38 collapses (precision 64-67 % to 98-99.7 %). The literature measures variants, not per-read placement, so our per-mate, matched-sensitivity metric is not directly comparable but has the same direction.
- levioSAM2's "suppress" category (no target counterpart) is the formal version of our "origin missing from GRCh38" class. Its policy is to drop such reads, the same abstention as our MAPQ 0 when CHM13 disagrees. Its "defer + reconcile by edit distance" is an alternative to our veto. A three-way rule worth testing: lift and commit, demote to MAPQ 0, or suppress. Applying levioSAM2's commit thresholds (MAPQ ≥30 or AS ≥100) to the CHM13 opinion would be a cheap ablation.
- The levioSAM2 simulation finds almost no gain (0.08 % on chr20) while real data shows large variant gains. Simulating reads from GRCh38 itself cannot show missing-sequence errors, which supports our choice of an HG002-assembly truth over GRCh38-based simulation.
- Masking GRCh38 false duplications (GIAB masked v2) is a near-free fix for KCNE1, CBS, CRYAA, H19 and U2AF1. If CERTA's GRCh38 index is not already the GIAB-masked analysis set, that is an obvious check.

### Gaps
- I did not retrieve Aganezov Table S15 (Illumina-specific FP/FN counts in the 269 genes) or the per-read MAPQ distributions. The paper reports mismatch rate and properly-paired fraction, not MAPQ.
- I did not fetch the full Nurk 2022 text beyond the abstract; Vollger's short-read mappability numbers (e.g., a 94.8 % short-read accessibility figure seen in a snippet) could not be confirmed.
- levioSAM2 reports no per-read high-MAPQ error rate on real data.

## 2. Pangenome references for short reads (HPRC 2023/2025, Giraffe 2021/2024/2025, haplotype sampling, DRAGEN multigenome, Minigraph-Cactus): measured gains, evaluation, survival of surjection, settings

### Takeaway
Published pangenome gains for short reads are variant-calling gains: 34 % fewer errors (HPRC 2023), 35 % fewer than BWA-MEM+DV with diploid haplotype sampling, 9.1 % fewer than the 10 %-AF-filtered graph. They were measured on GIAB v4.2.1 and CMRG in GRCh38 coordinates after surjection, with DeepVariant models retrained on graph alignments. The CHM13-based graph beat the GRCh38-based graph even when calls were projected to GRCh38, but only after excluding GRCh38 false-duplication and collapse regions, which "cannot be unambiguously projected". No published study scores per-read placement of real reads against an assembly truth. Giraffe is documented to leave more reads unmapped than linear mappers (reads the linear mappers misplace), and vendor benchmarks (DRAGEN 2024) dispute Giraffe+DV's lead.

### Cited Findings
**HPRC draft pangenome (Liao et al. 2023, Nature)**
- Abstract: the pangenome "reduced small variant discovery errors by 34%" and increased SVs detected per haplotype by 104 % against GRCh38-based workflows — [Liao et al. 2023](https://www.nature.com/articles/s41586-023-05896-x)
- Graph and processing:
  - Graph: the Minigraph-Cactus "clip.d9.m1000.D10M.m1000" allele-filtered graph, with nodes traversed by fewer than 9 haplotypes removed (minimum AF ~10 %). HG002, HG005 and NA19240 were held out of the graphs (90 haplotypes).
  - Mapping and realignment: Giraffe from vg 1.37.0, then `vg surject --prune-low-cplx` onto GRCh38 chr1-22, X and Y, then bamleftalign and ABRA realignment.
  - Calling: DeepVariant 1.3, retrained on Giraffe alignments with GIAB v4.2.1 truth (chr20 held out; HG003 held out for the evaluated model), with min MAPQ 1 in make_examples.
  - Comparators: BWA-MEM 0.7.17 on GRCh38 "with decoys but no ALTs", and DRAGEN 3.7.5 — [Liao 2023, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10172123/)
- Results:
  - HG005 at 30x, confident regions: Giraffe+DV made 21,700 errors (FP+FN), BWA-MEM+DV 36,144, DRAGEN 3.7.5 26,852, and DeepTrio with Giraffe 10,098.
  - CMRG v1.0: SNP F1 0.985 for Giraffe+DV against <0.976 for the others; indel F1 0.961 against <0.958. The largest gains were in regions with GRCh38 errors and large L1HS.
  - Simulated tandem-repeat reads: 6.6-8.5x more reads left unmapped by BWA-MEM/GRCh38 than by Giraffe — [Liao 2023, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10172123/)

**Minigraph-Cactus (Hickey et al. 2023, Nat Biotechnol): CHM13-based vs GRCh38-based graph, both projected to GRCh38**
- About 30x Illumina on HG001/HG002/HG005: 78.1 % (GRCh38-based graph) and 78.9 % (CHM13-based) of reads aligned perfectly, against 68.7 % for BWA-MEM/GRCh38. Reads with MAPQ >0 fell slightly, from about 94.9 % to 94.1 % — [Hickey et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906/)
- CMRG F1 after projecting to GRCh38: CHM13-based graph 0.9830, GRCh38-based graph 0.9823, DRAGEN 0.9777, BWA-MEM+DV 0.9756. The authors found "the CHM13-based graph to be superior" and that it "was able to improve upon the state-of-the art accuracy of the GRCh38-based graph, even when making calls on GRCh38" — [Hickey 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906/)
- Caveat: for the CHM13-based graph, regions falsely duplicated or collapsed in GRCh38 were removed from the benchmark because they "cannot be unambiguously projected from their corrected alleles in CHM13". The authors also say "such projections back to a linear reference can be lossy, especially in complex regions". The 10 % AF filter was chosen "to maximize variant calling sensitivity and mapping speed for the Giraffe/DeepVariant pipeline" — [Hickey 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906/)
- HPRC v1.0 resource notes: small-variant accuracy with Giraffe-DV improves if alleles with frequency <10 % are removed before indexing, and the released graphs contain spurious large deletion edges that slow Giraffe — [HPRC v1.0 MC notes (search snippet)](https://github.com/human-pangenomics/hpp_pangenome_resources/blob/main/hprc-v1.0-mc.md)

**Giraffe (Sirén et al. 2021, Science)**
- Giraffe maps short reads to haplotype paths in a GBWT-indexed graph about as fast as linear mappers, and DeepVariant was more accurate on Giraffe-to-pangenome alignments than on single-reference alignments. It was used to genotype about 167k SVs in 5,202 genomes. I did not retrieve the simulated per-read accuracy numbers — [Sirén et al. 2021](https://research.google/pubs/pangenomics-enables-genotyping-of-known-structural-variants-in-5202-diverse-genomes/)

**Personalized pangenomes / haplotype sampling (Sirén et al. 2024, Nat Methods)**
- Method: read k-mers are counted (KMC, k=29). The graph is split into ~10 kbp blocks, and in each block a greedy algorithm picks the best-scoring local haplotypes. "Diploid sampling" picks 2 haplotypes from 32 candidates — [Sirén et al. 2024, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760139/)
- Design: HG001-HG005 30x NovaSeq, GIAB v4.2.1, Giraffe+DeepVariant, all graphs GRCh38-based. "The constructed pangenomes include no GIAB sample or close relative" — [Sirén 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760139/)
- Small-variant results:
  - The diploid-sampled v1.1 graph had 9.1 % fewer total errors than the v1.1 AF-filtered graph and "consistently outperformed the other graph configurations".
  - Against BWA-MEM/GRCh38+DV it had 35 % fewer errors on average; against GATK best practice 74 % fewer (about 4-fold).
  - Sampling 4 or 8 haplotypes also beat the filtered graph, but accuracy fell as more haplotypes were sampled. Sampled graphs without diploid selection gave more exact and gapless alignments but fewer properly paired and fewer MAPQ-60 alignments, which the authors read as some extra alignments being wrong — [Sirén 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760139/)
- SVs (HG002, GIAB Tier1 v0.6, PanGenie): F1 rose from 0.7780 (v1.0 filtered) to 0.8926 (v1.1 diploid), against 0.29 for the best short-read SV discovery method (MetaSV). Overhead per sample: k-mer counting 432 s, sampling 341 s, index construction 997 s — [Sirén 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760139/)

**Giraffe 2025 (Chang et al.; HPRC release 2 graphs with more than 450 haplotypes; preprint)**
- Simulated short reads (1M Illumina and Element reads with real error profiles), linear comparators BWA-MEM and minimap2 on CHM13: Giraffe had "the fewest incorrectly mapped reads with a mapping quality (MAPQ) of 60" and "generally had better mapping quality calibration" — [Giraffe 2025 preprint, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633584/)
- Short-read Giraffe "could leave noticeably more reads unmapped than other mappers"; of the simulated reads it left unmapped, "all were incorrectly mapped by BWA-MEM and Minimap2". Haplotype sampling used 16 haplotypes. Short-read variant calling was not re-evaluated, because it "has not changed with the new release of Giraffe" — [Giraffe 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633584/); [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.29.678807.full.pdf)
- Discussion: graph construction can affect linear projections, and a "pathological 'spelling'" of a graph alignment can hurt downstream results — [Giraffe 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633584/)

**HPRC release 2 (2025) resources**
- Release 2 covers 232 individuals (464 assemblies; HG00272 excluded for a chrX misassembly), with improved polishing: about 3-fold fewer misassemblies and roughly halved base errors — [HPRC Release 2](https://humanpangenome.org/hprc-data-release-2/); [HPRC resources](https://github.com/human-pangenomics/hpp_pangenome_resources)
- The v2.1 Minigraph-Cactus graphs come in GRCh38-, CHM13- and GRCh37-based versions, each as default (clip), full, d46 AF-filtered (nodes/edges supported by <10 % of haplotypes dropped) and "eval" (HG002, HG005 and NA19240 excluded). The page says the AF-filtered indexes "give the best performance for Giraffe", and that `.hapl` files (default graphs only) let Giraffe infer a personalized pangenome as an alternative. It also says the files "have not been fully QC'd" — [HPRC pangenome resources](https://github.com/human-pangenomics/hpp_pangenome_resources)
- I found no published short-read accuracy benchmark specific to the release 2 graphs; the Giraffe 2025 preprint is the closest — [search; HPRC release 2 pages](https://humanpangenome.org/hprc-data-release-2/)

**Recommended Giraffe settings (vg wiki)**
- Haplotype sampling: `kmc -k29 -m128 -okff -t16 -hp sample.fq.gz ...`, then `vg giraffe -Z graph.gbz --haplotype-name graph.hapl --kff-name sample.kff -N sample -i -f reads.fq.gz`. Since vg 1.54.0 this does diploid sampling from 32 candidates with `--include-reference` by default. `--include-reference` is "necessary if the sampled graph is intended for a workflow based on a linear reference genome". Read coverage should be ≥20x, and the graph should be a Minigraph-Cactus clip graph in GBZ — [vg wiki: Haplotype Sampling](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
- With few haplotypes (e.g., 10) the default clip graph is usually best; with many more, the filter graph may be better. Giraffe can write BAM directly (`-o BAM`), projected to the nearest reference path (`--ref-paths` selects them). For performance supply an XG graph (`-x`). `--rescue-algorithm none` speeds mapping to complex graphs, with no accuracy statement — [vg wiki: Giraffe best practices](https://github.com/vgteam/vg/wiki/Giraffe-best-practices)
- After surjection, indel left-alignment (bamleftalign or `vg surject --left-align`) and ABRA realignment are recommended, as in the HPRC paper — [vg wiki: Mapping short reads with Giraffe](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
- vg autoindex applies proportional haplotype sampling by default when a graph has more than 192 haplotypes, and frequency filtering or sampling need many haplotypes to matter (J. Sirén on Biostars; forum, not peer-reviewed) — [Biostars](https://www.biostars.org/p/9614091)

**DRAGEN multigenome (Behera et al. 2024, Nat Biotechnol; vendor-authored)**
- Construction: GRCh38 plus 64 haplotypes (32 samples) of phased SNVs and indels, restricted to 366 Mbp of the hardest-to-map hg38 regions, excluding variants in <5 % of haplotypes.
  - Isolated SNVs are injected into hg38 as IUB ambiguity codes; indels and phased SNV clusters become alt contigs with known liftover alignments.
  - Alignments are grouped into "liftover groups", and the winning group's primary-contig member is reported, so the BAM is in hg38 coordinates.
  - MAPQ comes from the best vs second-best pair score — [DRAGEN preprint v2, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802302/)
- HG002, GIAB v4.2.1:
  - DRAGEN: F 99.85 %, 2,735 FP, 8,722 FN.
  - DeepVariant+Giraffe: 4,980 FP, 15,021 FN (1.74x DRAGEN's errors).
  - DeepVariant+BWA: 3,695 FP, 24,090 FN (2.42x).
  - GATK+BWA: 38,622 FP, 29,163 FN (5.91x).
  - CMRG F: DRAGEN 98.64 %, DV+Giraffe 98.10 %, DV+BWA 97.32 %, GATK 95.84 % — [DRAGEN](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802302/)
- DRAGEN applies corrections for falsely duplicated and collapsed genes in calling, e.g., 1,297 CBS variants and 341 H19 variants in 1kGP, where GATK found none for H19 — [DRAGEN](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802302/)
- Illumina says the second-generation multigenome mapper uses a 128-sample, 26-ancestry pangenome and has 40 % fewer errors than DRAGEN v4.2 (vendor claim, not peer-reviewed) — [Illumina](https://www.illumina.com/science/genomics-research/articles/second-gen-multigenome-mapping.html)
- Conflict: HPRC 2023 reports Giraffe+DV beating DRAGEN 3.7.5 (HG005: 21,700 vs 26,852 errors), while the DRAGEN paper reports DRAGEN 4.x with multigenome beating Giraffe+DV by 1.74x (HG002). The software versions, DeepVariant models and samples differ — [Liao 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10172123/); contradicted by [Behera 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802302/)

### Inferences
- Our Giraffe 1.78 result (94 vs 89 per million misplaced confident mates; 21,000 truth mates unmapped) does not contradict the literature. The published gains are variant-level, depend on DeepVariant retraining, surjection plus ABRA, and exclusion of GRCh38 false-dup/collapse regions, and nobody has published per-read real-data placement accuracy against an assembly truth. The 21,000 unmapped mates fit the documented behaviour that Giraffe leaves reads unmapped which linear mappers misplace. In a matched-sensitivity comparison those count against Giraffe even though they may be correct abstentions.
- Configuration choices likely left accuracy on the table:
  - We used the GRCh38-based eval graph. The literature (Hickey 2023) prefers the CHM13-based graph even for GRCh38-coordinate output.
  - We used d46 filtering rather than diploid haplotype sampling, which Sirén 2024 found has 9.1 % fewer variant errors and more properly paired / MAPQ-60 alignments than other sampling variants.
  - Candidate next run: `hprc-v2.1-mc-chm13-eval` (or the GRCh38 eval) default clip graph, `--haplotype-name .hapl --kff-name` with diploid sampling and `--include-reference`, then (a) surject to GRCh38 paths or (b) surject to CHM13 and lift with levioSAM2 or our own lift-over.
- Mechanism behind "surjection keeps some but not all gains": a mate drawn from sequence absent from GRCh38 can be placed correctly in the graph, but surjection must project it onto some GRCh38 path, typically the closest paralog, or leave it unmapped. In GRCh38 coordinates the only correct output for such mates is abstention (MAPQ 0 or unmapped), which is what our CHM13 second opinion does. The graph's value in GRCh38 space is therefore largely a better confidence signal, not a better coordinate. Hickey 2023 excluding those regions is a published admission of this.
- DRAGEN's "liftover groups" (score all alt haplotypes of a locus as one group, take MAPQ from the best vs second-best group) is a production-proven way to stop alt haplotypes deflating MAPQ while letting real paralogs deflate it. It maps directly onto a CERTA design in which T2T or HPRC haplotypes are alt contigs with known lift-over.

### Gaps
- No source gives per-read mapping accuracy of real short reads against an assembly truth for Giraffe, either surjected or in graph space. The Giraffe 2021/2025 per-read numbers are simulated, and the exact table values (Tables S1-S4 of the 2025 preprint) were not retrieved.
- No published HPRC release 2 short-read variant-calling benchmark was found. Which v2 graph (d46, clip+sampling, CHM13 vs GRCh38) gives the best short-read accuracy is undocumented beyond the resource page's "AF-filtered give the best performance" (likely speed and accuracy, unspecified).
- I did not retrieve the HPRC release 2 flagship paper text (if published) or its Giraffe section.

## 3. Paralog-aware methods (Parascopy, Locityper, ctyper, PSV-based gene callers, DRAGEN targeted callers, SDrecall): mechanisms for placing or reassigning reads, and accuracy

### Takeaway
Paralog-resolving tools do not trust the aligner's placement. They pool all reads, including MAPQ-0 reads, aligned to any copy of a duplication, then decide copy number and alleles from depth at a few paralog-differentiating sites (PSVs/SUNs) or from k-mers against pangenome haplotypes. Gene-specific callers (SMN1/2, CYP2D6, GBA) reach 96.5-99.8 % concordance. Genome-scale haplotype or k-mer methods (Locityper, ctyper) reach 95 %+ of haplotypes at QV ≥33 on up to 194 of 256 CMRG loci, against 22 loci for standard pipelines.

### Cited Findings
**PSV and depth-based single-gene callers (Illumina-authored)**
- SMN1/SMN2 (Chen et al. 2020, Genet Med): copy number from read depth at 8 informative SMN1/SMN2 differentiating sites chosen for consistency across African and non-African populations. SMN1 copy number was 99.8 % concordant with orthogonal methods and SMN2 99.7 %; SMA recall was 100 % and carrier recall 97.8 %, with 100 % precision. A 2019 preprint reported 48/48 SMN1 and 47/48 SMN2 agreement with digital PCR — [Illumina summary of Chen 2020](https://www.illumina.com/science/genomics-research/articles/sma-carrier-screening-wgs-data.html); [medRxiv](https://medrxiv.org/content/10.1101/19006635v2)
- Cyrius (CYP2D6; Chen et al. 2021): read counts include all reads mapped to CYP2D6 or CYP2D7, including MAPQ 0. Variants in homology regions are counted at both the CYP2D6 site and the matching CYP2D7 site, and unique regions use CYP2D6 reads only. Copy number is estimated at each of 117 CYP2D6/CYP2D7 differentiating bases, and hybrids are called where copy number changes along the gene. Concordance was 96.5 % against 84-86.8 % for existing methods (138 GeT-RM samples plus 8 HiFi samples, 40 star alleles), later improved to 99.3 % after truth corrections — [Cyrius, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805)
- Gauchian (GBA/GBAP1; Toffoli et al. 2022, Commun Biol): depth is measured over a 10 kbp unique region (GC-normalized; Gaussian mixture gives total GBA+GBAP1 copy number), and GBA copy number is estimated at each of 82 reliable differentiating sites, with breakpoints where it changes. In the exon 9-11 homology, 10 differentiating sites are phased by read and read-pair links to find hybrids and gene conversions. All 42 validation samples were consistent with ONT. In the AMP-PD cohorts Gauchian found 86 GBAP1-like variants against 44 for BWA-GATK, which produced 11 false positives for p.A495P. Runtime is about 1 min per sample. It "builds upon" the SMN1/2 and CYP2D6 callers — [Gauchian, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9259685); [GitHub](https://www.github.com/Illumina/Gauchian)
- DRAGEN targeted callers cover SMN, GBA, CYP2D6, CYP2B6, CYP21A2, HBA, LPA (KIV-2 copy number), RHD/RHCE and HLA. The paper gives mostly HG002 calls and 1kGP counts rather than accuracy (e.g., 1.53 % SMA carriers in 1kGP; GBA aggregate copy number 4 in 95.47 % of samples) — [DRAGEN](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802302/)

**Genome-scale duplication methods**
- Parascopy (Prodanov & Bansal 2022, Nat Commun) estimates aggregate and paralog-specific copy number for duplicated genes from WGS. It pools reads mapped to all copies of a low-copy repeat and uses PSVs as copy tags, with multiple samples used to learn which PSVs are reliable and to resist sequencing bias. It was evaluated on 2,504 1KGP genomes, where it is "robust to sequencing bias" with "higher accuracy compared to existing methods". The companion ParascopyVC calls pooled and locus-specific variants within SDs — [Parascopy, Nat Commun](https://www.nature.com/articles/s41467-022-30930-3); [IDEAS record](https://ideas.repec.org/a/nat/natcom/v13y2022i1d10.1038_s41467-022-30930-3.html); [eScholarship](https://escholarship.org/uc/item/7nx4c702)
- Locityper (Prodanov et al., Nat Genet 2025) — [bioRxiv full text](https://www.biorxiv.org/content/10.1101/2024.05.03.592358.full); [Nat Genet record](https://cris-portal.uni-muenster.de/portal/en/publication/191951483)
  - Inputs: candidate haplotypes come from a pangenome VCF or FASTA. Reads are recruited to all target loci at once by minimizers, with short reads needing matching minimizers at both ends.
  - Model: each candidate diploid genotype is scored on alignment error profile, insert size (negative binomial) and read depth, with down-weighted low-complexity windows.
  - Read placement: each pair's possible placements form a probability distribution over haplotype windows. Pairs with fewer than 5 target-specific k-mers are spread uniformly, and unplaceable reads go to a penalized null window. Assignment uses simulated annealing (default), greedy search or ILP.
- Locityper accuracy on 256 CMRG-type loci with 40 HPRC Illumina samples:
  - Full database: 94.9 % of haplotypes at QV ≥33.
  - Leave-one-out (the sample's own haplotypes removed): 59.9 % at QV >33 and 86.4 % at QV ≥23. Only 66.8 % of haplotypes even had a closest available haplotype at QV ≥33.
  - The 1KGP NYGC call set reached 27.4 % at QV ≥33; 194 loci reached 95 % at QV ≥33 with Locityper against 22 with standard calls. F1 was 0.939 against 0.842 for PanGenie.
  - Simulated reads: 98.6 % at QV ≥33. Runtime about 1 h 10 min per sample at 8 threads.
  - SMN1 and NCF1 were excluded (pangenome bubble over 300 kb) — [Locityper](https://www.biorxiv.org/content/10.1101/2024.05.03.592358.full); [search summary of Nat Genet abstract](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844405/)
- ctyper (Ma & Chaisson, Nat Genet 2025) is alignment-free. Read k-mer counts are fitted to a per-gene matrix of k-mer content across all similar pangenome sequences (paralogs and orthologs; 118 assemblies), then rounded to integer copy numbers along the phylogeny of allele sequences. Across 3,351 CNV genes (including HLA, SMN and CYP2D6) and 212 CMRG genes it recovers 96.5 % of phased variants with ≥99.1 % copy-number correctness, and 97.6 % F1 in 641 1KGP trios — [ctyper preprint, PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11326217/); [USC note on Nat Genet publication](https://www.qcb-dornsife.usc.edu/post/walfred-ma-published-a-paper-in-nature-genetics)
- SDrecall (Genome Biology, 2025/2026) builds an SD catalog, groups homologous copies, and realigns each group's reads against all homologous copies as a complement to GATK or DeepVariant. It reports about 95 % sensitivity against GIAB long-read-derived truth with low false positives. The authors say short-read SD ambiguity "cannot be completely resolved" — [SDrecall, HKU repository](https://hub.hku.hk/handle/10722/372497)

### Inferences
- Three distinct reassignment mechanisms emerge:
  - (a) Pool-then-decide at PSVs (SMN, Cyrius, Gauchian, Parascopy). The aligner's choice of copy is ignored, and allele evidence at differentiating sites decides copy number and haplotype.
  - (b) Haplotype likelihood with joint read redistribution (Locityper). Reads are probabilistically assigned across candidate haplotype windows, with depth and insert-size consistency as global constraints.
  - (c) Alignment-free k-mer projection (ctyper).
  All three treat MAPQ-0 or misassigned reads as recoverable evidence rather than noise. For CERTA --bio, a per-mate certificate could borrow (a): for mates in a known SD family, compute likelihood at the family's PSVs for each copy, with copies including CHM13/HPRC-only paralogs that are missing from GRCh38, and emit the copy only if the PSV likelihood ratio clears a threshold; otherwise abstain.
- Locityper's leave-one-out drop (94.9 % to 59.9 % at QV ≥33) shows that haplotype-database methods depend on the sample's haplotypes, or near-identical ones, being present. With HPRC eval graphs that exclude HG002, gains will be smaller than in published "full database" numbers. Matched-sensitivity measurements with HG002 held out, as ours were, are the honest setting.
- Depth constraints (copy number from a unique flanking region) are a global, cheap signal that per-read aligners ignore. Adding a depth-consistency check per SD family would catch "pile-ups" from collapsed copies, the mechanism Aganezov identifies for AF ≈ 0.5 false heterozygous calls.

### Gaps
- Parascopy's quantitative accuracy (concordance numbers, number of loci, comparison to QuicK-mer2 and others) was not retrieved; the Nature PDF redirected to a login.
- No paralog tool reports a per-read placement accuracy; all report genotype or copy-number accuracy.
- DRAGEN targeted-caller accuracy numbers were not in the parts of the paper retrieved.
- I found no 2018-2026 source quantifying PSV-likelihood read reassignment as a general mapper post-processor (as opposed to per-gene callers), apart from SDrecall's sensitivity figure.

## 4. Personalized or sample-specific references (personal diploid alignment, impute-first, VCF2diploid/AlleleSeq lineage) and reference bias

### Takeaway
Personalized diploid references give small, consistent variant-calling gains. Impute-first (Genome Research 2026) reaches HG002 F1 99.77 %, against 99.62 % for a linear aligner, 99.72 % for a pangenome graph aligner and 99.75 % for a personalized graph. Haplotype-sampled ("personalized") pangenomes beat the full or filtered graph. No study found quantifies per-read placement-error reduction in SDs from personalization.

### Cited Findings
- Impute-first (Vaddadi, Mun, Langmead; Genome Research 2026): reads are subsampled to genotype the donor, a phased diploid personal reference is imputed from a panel, and all reads are aligned to it with a linear (BWA-MEM, Bowtie 2) or graph (Giraffe) aligner. HG002 variant-calling F1 was 99.77 %, against 99.62 % (linear aligner), 99.72 % (graph pangenome aligner) and 99.75 % (personalized graph), with the best F1 coming from BWA-MEM on the personalized reference. The 2023 preprint had different numbers (F1 99.45 % vs 99.28 % for a graph pangenome) — [Genome Research](https://genome.cshlp.org/content/36/4/740); [preprint, PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10705441/)
- Biastools (Genome Biology 2024) built a diploid personalized HG002 reference (bcftools consensus) and simulated ~30x Illumina-like reads. More inclusive graph references produced fewer biased sites, and end-to-end alignment reduced bias at indels compared with local aligners — [biastools](https://link.springer.com/10.1186/s13059-024-03240-8)
- Mapping to a personalized diploid genome improved mapping accuracy and reduced reference-allele bias in allele-specific expression without modifying the aligner (Yuan & Qin; RNA-seq, pre-2018) — [PMC4304670](https://pmc.ncbi.nlm.nih.gov/articles/PMC4304670)
- Sirén 2024 is effectively a personalized-reference study in graph space: diploid sampling of 2 of 32 candidate haplotypes per 10 kbp block beat the full and filtered graphs, and accuracy fell as more haplotypes were included — [Sirén 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760139/)
- A near-complete diploid T2T HG002 benchmark (Q100) was reported in 2026, with the stated aim of less reference-biased benchmarking — [News-Medical report](https://www.news-medical.net/news/20260810/New-human-genome-benchmark-captures-both-parental-copies-with-near-perfect-accuracy.aspx)

### Inferences
- The consistent finding (Sirén 2024; impute-first) is that fewer, more relevant haplotypes beat more haplotypes: spurious alternative paths create false placements and MAPQ deflation. For a GRCh38-coordinate placer, an imputed personal diploid of the sample, plus CHM13/HPRC copies for SD families missing from GRCh38, is the natural "personal" second opinion. Its benefit for paralog placement depends on whether the imputation panel carries SD-copy haplotypes, which SNP-imputation panels mostly do not.
- Personalization addresses allele-level reference bias (diverged origin) more than missing copies (origin absent from GRCh38). Our error budget is dominated by both, so personalization alone would at best address the "diverged" part.

### Gaps
- AlleleSeq/vcf2diploid (2011) are outside the date window and were not fetched. No 2018-2026 study was found that maps real reads to the sample's own diploid assembly (e.g., HG002 Q100) to measure the residual placement error ceiling, though that is a natural upper-bound experiment.
- The impute-first paper's stratified (SD, CMRG) numbers were not retrieved; the Genome Research PDF required a login.

## 5. Studies quantifying placement errors in segdups with assembly-based truth, comparing linear vs graph vs T2T

### Takeaway
I found no published study that scores per-read placement of real short reads against an assembly-derived truth (e.g., HG002 Q100) across linear GRCh38, T2T-CHM13 and graph references. The nearest designs are: CHM13 self-mapping to GRCh38 as a haploid assembly truth (Aganezov), variant-level GIAB/CMRG comparisons of CHM13-lift vs GRCh38 (levioSAM2) and of graph vs linear (HPRC, Hickey, Sirén, DRAGEN), and simulated per-read accuracy (Giraffe, levioSAM2). Our measurement appears to fill this gap.

### Cited Findings
- Aganezov 2022 used CHM13's haploidy as truth: CHM13 HiFi and simulated 150 bp Illumina-like reads mapped to GRCh38 gave 368,574 heterozygous SNVs, all errors or cell-line mutations, delineating 908 problematic regions (20.8 Mbp) — [Aganezov 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- levioSAM2 scored per-read placement only on simulated reads (correct within 10 bp; 0.08-1.37 % more correct on chr20/21) and variants on real reads (GIAB v4.2.1, CMRG) in GRCh38 coordinates — [levioSAM2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/)
- Giraffe 2025 per-read accuracy uses simulated reads with CHM13 as the linear comparator, and notes that reads Giraffe leaves unmapped were all mis-mapped by linear mappers — [Giraffe 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12633584/)
- Hickey 2023 compared CHM13-based and GRCh38-based graphs against BWA-MEM on real reads at the variant level in GRCh38, excluding regions falsely duplicated or collapsed in GRCh38 — [Hickey 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906/)
- SDrecall states that short-read SD ambiguity "cannot be completely resolved due to the inherently limited information carried by short reads" — [SDrecall](https://hub.hku.hk/handle/10722/372497)
- A 2025 analysis of low-complexity regions attributes most errors to inconsistent read alignment there, and notes that assembly-to-reference alignment, as used for HG002-Q100, gives superior truth — [search result, arXiv 2507.03718](https://arxiv.org/html/2507.03718v2)

### Inferences
- Because no published per-read, real-data, assembly-truth comparison exists, our numbers (89 per million confident misplacements for minibwa; CHM13 rescuing 38 % of the never-correct mates; Giraffe 94 per million with 21k unmapped) are novel evidence. They should be reported with the scoring choices made explicit (within 150 bp, matched sensitivity, truth via dipcall projection), because they differ from all published designs.
- The literature's variant-level gains (HPRC 34 %, Sirén 35 %, levioSAM2 24 %) are consistent with our read-level finding that most residual errors are mates whose origin is missing from or diverged in GRCh38. Those are exactly the classes that T2T, graph and personal references target, and that GRCh38-based simulation cannot reveal (levioSAM2's simulated gain is only 0.08 % on chr20).

### Gaps
- I found no study using the HG002 Q100 assembly as read-origin truth for real short reads, and no published per-read error rate in SDs for Giraffe or DRAGEN multigenome on real data.
- I did not retrieve Giraffe 2021's simulated per-read precision-recall numbers or Giraffe 2025's Tables S1-S4, which would give simulated MAPQ-60 error rates per mapper.
