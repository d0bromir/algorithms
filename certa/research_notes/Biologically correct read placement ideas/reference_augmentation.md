# Reference augmentation for short-read alignment to GRCh38: decoys, ALT-aware mapping, masked false duplications, T2T-CHM13 / pangenome sequence in GRCh38 coordinates

Scope note for the report writer: the user's measurement (HG002 2x150 NovaSeq X, GRCh38 Ensembl primary with no decoys and no ALTs, Q100/dipcall truth, locus-level scoring at matched sensitivity) shows that about 89 confident misplacements per million mates come from copies that are missing from or diverged in GRCh38, mostly segmental duplications (SDs). A T2T-CHM13 second opinion lifted with hs1ToHg38 halves these errors. The notes below are organised around which published reference change reproduces that gain, ideally in one index and one pass.

Two citation corrections: (1) levioSAM2 is Chen, Paulin, Sedlazeck, Koren, Phillippy, Langmead, *Nature Methods* 2024 (DOI 10.1038/s41592-023-02069-6), not "Mun et al. 2023 Genome Biology". Mun et al. wrote the original levioSAM (Bioinformatics 2021). (2) Two different first authors share the surname Behera: FixItFelix (masked false duplications plus T2T decoys) is Behera et al. *Genome Biology* 2023, and the DRAGEN paper is Behera et al. *Nature Biotechnology* 2024.

---

## 1. Decoy sequences (hs37d5, hs38d1, analysis sets) and pangenome or CHM13 sequence used as decoys

### Takeaway
Classical decoys work by mechanism: reads from sequence missing from the primary assembly land on the decoy instead of being forced onto a paralog. They remove a large share of false heterozygous calls on GRCh37 (about 2x fewer with hs37d5) and a modest share on GRCh38 (about 15% fewer FPs with hs38DH, at a cost of about 0.6% of TPs). The only quantified pangenome-novel decoy study on GRCh38 (46.6 Mb of common Chinese pan-genome sequence) lowered GIAB FPs by only 3.3 to 5.2%. The decoys that target the user's failure mode (missing SD copies, resolved in T2T-CHM13) have been published only for a handful of genes (FixItFelix, Section 3).

### Cited Findings
**What decoys are and what they did (GRCh37 era)**
- hs37d5 (1000 Genomes phase 2 reference) consists of the GRCh37 primary assembly (chromosomes plus unlocalized and unplaced contigs), rCRS chrM (NC_012920), human herpesvirus 4 type 1 / EBV (NC_007605), and the concatenated decoy hs37d5cs.fa.gz, with 20 N bases between adjacent sequences. The decoys were built from GRCh37.p4 patch sequences, primate GenBank BAC/fosmid/HTG clones (excluding LOW.PASS), HuRef contigs, and the NA12878 ALLPATHS-LG assembly (AEKP). "Each sequence represents an unmapped segment," extracted by a pipeline that collected sequence unmapped to the GRCh37 primary — [1000G README_human_reference_20110707](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707)
- Li (2014) describes the hs37d5 decoy as 35.4 Mb of sequence "present in several de novo assemblies but likely to be missing from the primary assembly of GRCh37". These sequences "are supposed to attract many mismapped reads" — [Li 2014, arXiv 1404.0929 / Bioinformatics](https://arxiv.org/pdf/1404.0929)
- Quantified effect (Li 2014, haploid CHM1, so every heterozygous call is an error; BWA-MEM plus GATK HaplotypeCaller, autosomal heterozygous SNPs after the low-complexity filter, Fig. 6). Homozygous non-LC SNP counts were nearly identical across references: 2.408 M (GRCh37), 2.405 M (hs37d5), 2.412 M (GRCh38). "We called twice as many heterozygotes from GRCh37 in comparison to hs37d5," which shows that the 35.4 Mb decoy "indeed attracted many mismapped reads". The Venn regions are: GRCh37-only 132k; GRCh37 ∩ hs37d5, not GRCh38, 36.9k; all three 80.3k; GRCh37 ∩ GRCh38, not hs37d5, 36.8k; GRCh38 adds 24.6k new hets that map to v37 autosomes plus 73.0k hets in sequence unmapped to v37 autosomes. "GRCh38 further resolves 39.8k (=36.9k+2909) heterozygotes called from hs37d5," while it "retains 36.8k" GRCh37-only false hets that the decoy removed — [Li 2014](https://arxiv.org/pdf/1404.0929)
- Li (2014) also found that the maximum-depth filter (d+3√d to d+4√d) was the most effective filter against false heterozygotes. These "false positives are mostly caused by copy number variations (CNVs) or paralogous sequences not present in the human reference genome". The raw error rate was about 1 per 10–15 kb, falling to about 1 per 100–200 kb after filtering — [Li 2014](https://arxiv.org/pdf/1404.0929)

**GRCh38 decoys (hs38d1) and analysis sets**
- NCBI analysis sets (GCA_000001405.15) come in four files: `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz` (chromosomes, chrM, unlocalized, unplaced, EBV); `..._full_analysis_set.fna.gz` (adds the ALT_REF_LOCI alt scaffolds); `..._full_plus_hs38d1_analysis_set.fna.gz` (adds "human decoy sequences from hs38d1 (GCA_000786075.2)"); and `..._no_alt_plus_hs38d1_analysis_set.fna.gz`. In all four, the chrY PARs are hard-masked with N, and "duplicate copies of centromeric arrays and WGS on chromosomes 5, 14, 19, 21 & 22" are hard-masked with N. EBV is included "as a sink for alignment of reads". The BWA index tarballs of the full sets include the `.alt` file — [NCBI README_analysis_sets.txt](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt)
- hs38d1 is about 2,385 contigs (JTFH*, KN* accessions) according to a secondary source (Park lab decoy analysis repo). No total size was found — [parklab/reference_genomes_decoy_analysis](https://github.com/parklab/reference_genomes_decoy_analysis)
- bwakit's hs38DH builds its decoys "from GenBank clones and the de novo assembly of 254 public SGDP samples" (file `hs38DH-extra.fa`) and adds HLA alleles from IMGT/HLA v3.18.0. bwakit generates it with `bwa.kit/run-gen-ref hs38DH` — [bwa README-alt.md](https://raw.githubusercontent.com/lh3/bwa/master/README-alt.md)
- Quantified decoy effect on GRCh38 (bwa README-alt, CHM1/NA12878-style evaluation):

  | Reference | hs37 | hs38 | hs38DH | CHM1_1.1 | HuRef |
  |---|---|---|---|---|---|
  | FP | 255,706 | 168,068 | 142,516 | 307,172 | 575,634 |
  | TP | 2,142,260 | 2,163,113 | 2,150,844 | 2,167,235 | 2,137,053 |

  hs38DH (decoy plus HLA plus ALT) "reduces FP by ~25k but also reduces TP by ~12k". For scale, the difference between two NA12878 HiSeq X10 replicates is "~140k, an order of magnitude higher" — [bwa README-alt.md](https://raw.githubusercontent.com/lh3/bwa/master/README-alt.md)
- 1000 Genomes GRCh38 reference: `GRCh38_full_analysis_set_plus_decoy_hla.fa` (3.0 GB uncompressed), with prebuilt BWA index files (`.bwt` 3.0G, `.sa` 1.5G, `.pac` 767M, `.ann` 438K, `.amb` 20K), `.fa.alt` (476K), `.fai`, `.dict`, and `-extra.fa` (7.9M). It was unpacked from bwakit-0.7.12 and is "based on GRCh38 full analysis set" plus hs38d1 plus HLA. Contig suffixes: `_alt`, `_decoy` (e.g. chrUn_KN707606v1_decoy), `HLA-*` — [1000G GRCh38_reference_genome directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/); [README via same directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/README.20150309.GRCh38_full_analysis_set_plus_decoy_hla)
- That file has 3,366 contigs and is the baseline from which DRAGEN's hg38 hash tables are built — [Illumina, Demystifying GRCh38 versions in DRAGEN](https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html)
- The 1000 Genomes realignment to GRCh38 (alt-aware bwakit-0.7.12 with bwa-postalt.js) raised the mapping rate from 92.6% to 96.2% for low-coverage WGS and from 93.6% to 97.5% for WES, compared with GRCh37 — [Zheng-Bradley et al. 2017 GigaScience](https://pmc.ncbi.nlm.nih.gov/articles/PMC5522380)

**Pangenome-novel sequence added to GRCh38 as decoys**
- Chinese pan-genome (CPG; Li et al. 2021, Commun Biol): 486 deep-sequenced Han Chinese genomes yielded 276,116,245 bp (612,914 sequences) absent from GRCh38. Of this, 46,645,513 bp (45,836 sequences, shared by at least 2 of 486 individuals) were added to GRCh38. More than 90% of previously unaligned reads mapped to these common sequences. GIAB results with BWA-MEM plus GATK HC at 60x:

  | Sample | FP, GRCh38 | FP, GRCh38 + CPG | FP change | TP, GRCh38 | TP, GRCh38 + CPG |
  |---|---|---|---|---|---|
  | HG001 | 12,819 | 12,206 | −4.78% | 3,478,662 | 3,478,639 |
  | HG002 | 16,575 | 16,025 | −3.32% | 3,443,252 | 3,443,250 |
  | HG005 | 9,895 | 9,427 | −4.73% | — | — |
  | HG005, 300x | 6,900 | 6,542 | −5.19% | — | — |

  FN counts were essentially unchanged. Data: GSA-Human PRJCA003657; code at github.com/HKU-BAL/CPG and Zenodo 10.5281/zenodo.5155074 — [Li et al. 2021, Building a Chinese pan-genome of 486 individuals](https://pmc.ncbi.nlm.nih.gov/articles/PMC8405635/)
- African Pan-Genome (APG; Sherman et al. 2019 Nat Genet): 296.5 Mb in 125,715 contigs that are not in GRCh38, assembled from reads of 910 individuals that did not align. Only 1,548 contigs (4.4 Mb) could be anchored to GRCh38. Accession PDBU01000000. No decoy-mapping benchmark was found — [Sherman et al. 2019 (PubMed)](https://pubmed.ncbi.nlm.nih.gov/30455414/)
- HUPAN (Duan et al. 2019 Genome Biol) found 29.5 Mb of novel sequence (at least 188 novel protein-coding genes) in 185 plus 90 Han Chinese genomes. It states that adding these sequences "could improve" mapping and variant calling, without a benchmark in the snippet seen — [Duan et al. 2019](https://dx.doi.org/10.1186/s13059-019-1751-y)
- A 2018 bioRxiv on two Swedish long-read assemblies (>10 Mb non-GRCh38 sequence each) reports that adding these sequences to GRCh38 "radically improves" alignment and variant calling "at several genomic loci". This is qualitative only — [bioRxiv 267062](https://www.biorxiv.org/content/early/2018/02/18/267062.full.pdf)
- CHM13-derived decoys on GRCh38: the only published construction found is FixItFelix/GIAB (Section 3), with decoys for 4 to 8 collapsed genes. No study was found that adds all CHM13 non-syntenic or collapsed-SD sequence (e.g. from `chm13v2-unique_to_hg38.bed`) to GRCh38 as decoys and benchmarks it genome-wide.

### Inferences
- Most decoy-type sequence (HuRef, SGDP assemblies, CPG/APG novel contigs) is non-reference insertion or unique sequence. The user's residual errors are the opposite case: extra paralog copies that GRCh38 lacks (collapsed SDs, 8.04 Mb over 203 loci per Aganezov 2022, Section 5). This would explain why population-novel decoys move FPs by only 3–5% while the CHM13 second opinion halves confident errors. To capture the user's gain in one index, the decoy content should come from CHM13 (or HPRC haplotypes) copies of collapsed SDs, not from population-novel insertions.
- Decoys reproduce the user's "MAPQ 0 if T2T places it >150 bp elsewhere" rule only if they are ordinary (non-ALT) contigs. bwa computes the MAPQ of non-ALT hits "across non-ALT hits only" ([README-alt](https://raw.githubusercontent.com/lh3/bwa/master/README-alt.md)), so a decoy listed in the `.alt` file would not lower the MAPQ of the competing GRCh38 hit. Check with `grep -c decoy *.fa.alt` before reusing the 1000G/bwakit indexes. (This is an inference from the documented MAPQ rule; whether bwakit lists decoys in `.alt` was not verified.)
- The Ensembl primary assembly the user maps to lacks EBV, hs38d1 and the NCBI centromere and PAR masks. Moving to `no_alt_plus_hs38d1` is a cheap control, but the bwa README-alt numbers (about 15% FP reduction, partly bought with TP loss) suggest it will not approach a 2x reduction in confident errors.

### Gaps
- Total size of hs38d1 (GCA_000786075.2) was not confirmed. The count of about 2,385 contigs comes from a secondary source.
- No study was found that measures read-level misplacement (rather than variant FPs) with and without hs38d1 on real NovaSeq data. The Park lab repo runs this comparison (7 references, 8 1KGP samples, simulated reads from hs37d5/hs38d1) but publishes no numbers on its front page.
- No study benchmarks GRCh38 plus whole-genome CHM13-derived SD decoys (beyond FixItFelix's handful of genes) or GRCh38 plus HPRC-derived decoys.

---

## 2. ALT-aware mapping (bwa .alt plus bwa-postalt.js, GATK, DRAGEN alt-aware, alt-masked and multigenome graph)

### Takeaway
The GRCh38 ALT contigs (261 scaffolds, 109 Mb) mainly address reference bias at a few highly polymorphic loci. Naively including them collapses MAPQ, and bwa's ALT-aware MAPQ rule prevents that collapse. DRAGEN replaced liftover-based ALT handling with alt-masking, which hard-masks ALT segments similar to the primary so that divergent segments act as decoys. DRAGEN's multigenome graph adds population haplotypes from pangenome samples, and reports about 47–50% fewer SNP errors and about 24–27% fewer indel errors than its own linear reference. The v4.3 version (128 samples) reports a further 40% error reduction over v4.2. These are the largest published one-pass, GRCh38-coordinate gains, but they are vendor-reported and mostly variant-level.

### Cited Findings
- GRCh38 has 261 ALT scaffolds covering 178 regions, 109 Mb in total, and much of this sequence is nearly identical to the primary. "Without ALT awareness, multi-mapping reads get MAPQ 0". The 1000 Genomes project wrote ALT alignments as separate BAM lines so that ALT variants could be called independently. It recommends realigning and recalling rather than lifting variants — [Zheng-Bradley et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5522380)
- bwa ALT-aware mechanism: BWA-MEM reads ALT names from `idxbase.alt` and ignores the ALT-to-ref alignment during mapping. "The mapQ of a non-ALT hit is computed across non-ALT hits only". If there is no non-ALT hit, the best ALT hit becomes primary; otherwise ALT hits get flag 0x800. `bwa-postalt.js` lifts ALT hits in the XA tag to chromosomal positions, groups overlapping hits and re-estimates MAPQ. It also supports HLA typing by de novo assembly, with no accuracy figures given and the caveat that it "may not always be successful". The ALT contigs "are totaled 109Mb in length, spanning 60Mbp of the primary assembly" — [bwa README-alt.md](https://raw.githubusercontent.com/lh3/bwa/master/README-alt.md)
- DRAGEN alt-aware (up to v3.8) moved reads from ALT regions with a valid liftover to the primary. Illumina found that "the correct liftover between a long ALT haplotype and the primary assembly is often ambiguous", and that wrong liftovers create dense mismapping clusters and false calls. From v3.9, alt-masking hard-masks (as N, via `/opt/edico/fasta_mask/hg38_alt_mask.bed`) the ALT segments similar to the primary "so they don't compete for alignments or squash MAPQs", while divergent segments stay unmasked "and act essentially as decoys". Build example: `dragen --ht-reference hg38.fa --ht-alt-aware-validate=true ... --ht-mask-bed /opt/edico/fasta_mask/hg38_alt_mask.bed --output-directory hg38_alt_masked`. Alt-masking reduces SNP and indel FP+FN compared with alt-aware in both graph and non-graph builds; the numbers appear only in figures. DRAGMAP (open-source) supports alt-masked but not graph references — [Illumina, Demystifying GRCh38 versions](https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html); [DRAGEN v4.2 docs, ALT contig hash table options](https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/RefGenComALT.htm)
- DRAGEN graph vs non-graph (NIST v4.2.1, HG001–HG007): "DRAGEN graph reduces the SNP errors by ~50% and INDEL errors by ~27%" — [Illumina, Demystifying](https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html)
- DRAGEN multigenome generations: v3.7 used 32 haplotypes from 16 European samples added to the native ALTs, with a "47% reduction in SNP errors and a 24% reduction in indel errors" compared with the linear reference including native ALTs. v4.2 used 64 haplotypes from 32 globally diverse samples. v4.3 used 128 samples from 26 ancestries (256 haplotypes), with all variants phased. Population SNPs are IUB-encoded in the primary for seeding, and complex variants still use ALT-style liftover groups. The best liftover-group score wins and the primary-contig member is reported, so output is in GRCh38 primary coordinates. MAPQ is estimated per liftover group "using population information"; no formula is given. Reported v4.3 results: a further 40% error reduction over v4.2 on HG002 (v4.2.1 all regions); 61.61% fewer errors than Giraffe-DeepVariant averaged over 7 GIAB samples (SNP 63.8%, indel 53.53%); 65.51% fewer than Giraffe-DV and 37.77% fewer than v4.2 in difficult-to-map regions; 15% fewer FP+FN on CMRG (HG002) than v4.2. Hash table `hg38-alt_masked.graph.cnv.hla.rna_v4` is on the DRAGEN product files page — [Illumina, Second-generation multigenome mapping](https://www.illumina.com/science/genomics-research/articles/second-gen-multigenome-mapping.html)
- Peer-reviewed summary: Behera et al. 2024 Nat Biotechnol (DOI 10.1038/s41587-024-02382-1) describes multigenome mapping with pangenome references as outperforming the state of the art across variant types on 3,202 WGS datasets. A Nature Genetics News & Views piece says the multigenome mapper "outperforms existing pangenome implementations" — [Illumina feature on DRAGEN germline algorithms](https://www.illumina.com/company/news-center/feature-articles/DRAGEN-germline-algorithms.html); [X-MOL record of Nat Biotechnol paper](https://x-mol.com/paper/1849816231849041920)
- Illumina's roadmap (in the same article) lists "hg38-alt-masked-V2", which adds GIAB false-duplication masking, and "alt-masked-V3", which adds "new decoys (GIAB/Baylor College of Medicine)" and expanded graph contigs. Illumina says V2 recovers all CMRG FPs and FNs in the affected chr21 genes (CBS, CRYAA, KCNE1), with counts shown only in a table image — [Illumina, Demystifying](https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html)

### Inferences
- The user's reference (Ensembl primary) has no ALTs, so ALT handling is not a cause of their errors. The ALT and graph gains come mainly from population haplotype diversity, i.e. reference bias at polymorphic sites. Their 89/M residual is dominated by missing SD copies, which DRAGEN's graph captures only if the pangenome sample haplotypes include those copies as ALT-like paths. The Illumina articles do not say whether they do.
- Alt-masking is effectively "keep only divergent ALT segments as decoys". The same trick applied to CHM13 would be: align CHM13 to GRCh38, mask the CHM13 segments that are near-identical to GRCh38 1:1 orthologs, and keep the extra paralog copies and divergent haplotypes as decoys. This is the natural one-index generalisation of both FixItFelix and DRAGEN alt-masking. It is an untested design.

### Gaps
- GATK's own documentation of ALT handling (gatk.broadinstitute.org "Reference Genome Components" and the alt-contig mapping tutorial) returned HTTP 403, so the GATK specifics are unverified. Not verified: the hg38 bundle `Homo_sapiens_assembly38.fasta` contig composition, and whether the GATK Best Practices pipeline runs bwa-postalt.js (commonly reported as not run).
- No independent (non-Illumina) read-level evaluation of DRAGEN multigenome vs linear in SDs was found. The size and contents of `hg38_alt_mask.bed`, and whether the v4.3 multigenome contains SD copy-number variants, are not published in the sources read.
- Whether alt-masked-V2 or V3 shipped in a numbered DRAGEN release was not confirmed.

---

## 3. GRCh38 with masked false duplications plus T2T decoys (GRC exclusions, FixItFelix, GIAB GRCh38-GIABv3)

### Takeaway
GRCh38 contains about 1.2 Mb of false duplications (12 regions, 22 protein-coding genes, including 5 misplaced chr21 clones covering CBS, CRYAA, KCNE1, U2AF1) and about 8 Mb of collapsed duplications. Masking the false copies and adding T2T-CHM13 copies of a few collapsed genes as decoys (FixItFelix; GIAB's `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18`) has large local effects. CMRG SNV F1 in the duplicated regions rose from 0.012 to 0.980, and MAPQ-0 reads there fell 78%. Collapsed regions gain MAPQ-0 reads (+72%) as the decoys absorb reads from the missing copies. The construction covers only about 2 Mb and about 30 genes, so it cannot capture most of a genome-wide 2x effect.

### Cited Findings
- FixItFelix (Behera et al. 2023 Genome Biol): GRCh38 has "1.2 Mbp of falsely duplicated and 8.04 Mbp of collapsed regions", which affect variant calling in 33 protein-coding genes, 12 of them medically relevant. The authors masked the extra copies of falsely duplicated regions and added T2T-CHM13 sequence for missing genes as decoys, keeping GRCh38 coordinates. A ~30x BAM is fixed in about 4–5 CPU minutes by remapping only the affected reads, versus about 24 CPU hours for full remapping — [FixItFelix, Genome Biology 2023 (PubMed 36810122)](https://www.pubmed.gov/36810122/); [DOAJ record](https://doaj.org/article/ddf9214ce77c4538831718d8d9803918)
- Details from the FixItFelix preprint (bioRxiv 10.1101/2022.07.18.500506v2):
  - **Regions masked:** 12 falsely duplicated regions containing 22 protein-coding genes and 18 pseudogenes, 2,032,012 bp in total (1,010,809 bp false copy). Genes include CBS, KCNE1, CRYAA, TRAPPC10, DNMT3L, U2AF1, ICOSLG, SIK1, SMIM11A, TRPM3.
  - **Collapsed regions:** 9 regions with 9 genes, 843,139 bp. The candidates were TAS2R46, MAP2K3, KCNJ18, KATNAL2, FANCD2, LPA, MUC3A, KMT2C, GPRIN2, DUSP22; TAS2R46, KATNAL2, LPA and MUC3A were excluded for lacking simple decoys.
  - **Decoy versions:** the main reference has decoys for KCNJ12, KCNJ18, KMT2C, MAP2K3; a second version adds FANCD2, FANCD2P1, GPRIN2, DUSP22.
  - **HG002 short-read MAPQ:** MAPQ=0 reads in duplicated regions fell from 358,644 to 103,392 (−78%). MAPQ=0 reads in collapsed regions rose from 111,685 to 192,280.
  - **CMRG, GATK, duplicated-region SNVs:** recall 0.007→1.0, precision 0.063→0.961, F1 0.012→0.980. Indels went from 0 to 1.0 (16 variants).
  - **KMT2C:** SNV recall 0.949→0.974, precision 0.064→0.098.
  - **dipcall-based truth:** duplicated-region SNV recall 0.13→0.85; collapsed-region precision and F 0.22/0.32→0.452/0.453; the FANCD2 pseudogene dropped from 0.941 to 0.318.
  - **WES:** duplicated-region SNV F1 0.029→0.945.
  - **4,174 samples (1KGP, HGDP, TOPMed):** mean depth in collapsed regions fell from 51.5 to 38.4 (genome background 33.7); depth in duplicated regions rose from 17.7 to 34.5; Hardy–Weinberg departures were reduced.
  - **eQTL:** 10,437 vs 10,417 eGenes.
  - **Files:** modified references on bcm.box.com (two links in the paper); code at github.com/srbehera/FixItFelix.

  Source: [FixItFelix bioRxiv v2 full text](https://www.biorxiv.org/content/10.1101/2022.07.18.500506v2.full)
- GIAB/GRC chr21 masking: GIAB "worked with the GRC to mask 5 GRCh38 false duplications on chr21", which "improves short read accuracy from 8% to 100%" in CBS, CRYAA and KCNE1 (CMRG benchmark announcement) — [GIAB Google group, CMRG benchmark](https://groups.google.com/g/genome-in-a-bottle/c/ppjHR5eDPJ8)
- GIAB FTP files (ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/):

  | File | Size | Date |
  |---|---|---|
  | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz` | 845M | |
  | `GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta.gz` | 845M | 2021-04 |
  | `..._maskedGRC_exclusions_v2.fasta.gz` | 845M | 2021-07 |
  | `GRCh38_masked_v2_decoy_gene.fasta.gz` (FixItFelix) | 845M | 2023-02 |
  | `GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.gz` | 787M | 2023-02 |
  | `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz` | 755M | 2023-09-27 |
  | `GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed` | 296 B | |
  | `GRCh38_reference_update_to_GIABv3.ipynb` (build notebook) | | |

  Each FASTA comes with `.fai` and `.gzi` files — [GIAB references/GRCh38 directory](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/)
- The GIABv3 build notebook starts from `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz`. It hard-masks false duplications and contaminating non-human sequence using `GCA_000001405.15_GRCh38_GRC_exclusions.bed`, adds decoys for MAP2K3, KMT2C and KCNJ18 taken from the FixItFelix file, converts non-ACGTN bases to N, and shortens headers to `>chr1` etc. — [GRCh38_reference_update_to_GIABv3.ipynb](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_reference_update_to_GIABv3.ipynb)
- T2T/Aganezov 2022 supporting numbers: "12 regions affecting 1.2 Mbp and 74 genes (including 22 protein-coding genes)" are falsely duplicated in GRCh38. "The five largest duplications in GRCh38, affecting 15 protein-coding genes on the q-arm of chromosome 21," were BAC clones misplaced from the acrocentric p-arm. Short-read Illumina recall in the corrected false duplications was 57.4–68.3% on T2T-CHM13 vs 1.1–1.8% on GRCh38 (precision 98.5–99.3% vs 76.5–95.8%). In collapsed duplications (<10 copies, ~910 kbp), precision was 98.1–99.7% on CHM13 vs 64.3–67.3% on GRCh38 at similar recall (20.1–28.3% vs 21.5–25.4%) — [Aganezov et al. 2022 Science (author manuscript)](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)

### Inferences
- The FixItFelix MAPQ pattern is the same mechanism as the user's T2T second opinion, restricted to about 2.9 Mb. Masking false copies turns MAPQ-0 reads into confident correct ones. Decoys for collapsed copies turn wrongly confident reads into MAPQ-0 ones: +80,595 MAPQ-0 reads in 843 kb, i.e. reads whose true origin is the missing paralog. Scaling the decoy set from 4–8 genes to the 203 collapsed loci (8.04 Mb) of Aganezov 2022 is the obvious one-index, one-pass, GRCh38-coordinate candidate for capturing most of the CHM13 second-opinion gain.
- The FANCD2-pseudogene recall drop (0.941→0.318) and the "excludes_GPRIN2_DUSP22_FANCD2" variant show the risk of over-broad decoys: a decoy that is closer to a true GRCh38 locus than to the real missing copy steals real reads. A genome-wide CHM13 decoy set needs per-locus validation, e.g. against HG002 Q100 truth, which the user already has.
- `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz` (755 MB) is the cheapest drop-in test. Map the HG002 NovaSeq X set with each aligner and count the drop in confident misplacements per million against the CHM13 second opinion. Expect a small fraction of the 2x gain, given its coverage of about 2 Mb.

### Gaps
- The exact coordinates and bp count in `GRC_exclusions.bed` and `GRC_exclusions_T2Tv2.bed`, the decoy contig names and lengths in GIABv3, and why GIABv3 dropped the KCNJ12/FANCD2/GPRIN2/DUSP22 decoys are not shown in the readable part of the notebook.
- No published genome-wide read-placement accuracy (as opposed to CMRG-local variant metrics) for the GIABv3 reference was found.

---

## 4. Multi-reference and lift-over consensus mapping (levioSAM/levioSAM2, reference flow, vg Giraffe surjection to GRCh38)

### Takeaway
levioSAM2 is the published form of the user's second opinion: map to T2T-CHM13, lift to GRCh38, and selectively re-map only low-confidence reads to GRCh38. With BWA-MEM on real 30x NovaSeq data it cut small-variant errors by 23.9% (GATK) and 11.4% (DeepVariant) compared with direct GRCh38. In CMRG genes the cut was 35.2% (GATK) and 20.0% (DeepVariant). It used about 55% of the CPU of direct mapping. Graph mapping with vg Giraffe on the HPRC pangenome, surjected to GRCh38, gives about 34% fewer small-variant errors on average than BWA-MEM/GRCh38 (HG005: 21,700 vs 36,144 errors).

### Cited Findings
- levioSAM2 workflow (Chen et al. 2024 Nat Methods): reads are aligned to the source (CHM13), lifted to the target (GRCh37 or GRCh38), and classified as committed, suppressed or deferred.
  - **Suppressed reads** are those in source-only regions (T2T centromeres, missing collapsed copies identified by aligning >5 kb source-only sequence to the target with Winnowmap2). They are left unaligned in the target.
  - **Deferred reads** are low-confidence lifts, including GenMap2 mappability-reduced regions (100-mers, 0.01 mismatch). They are re-mapped to the target.
  - **Commit thresholds:** BWA-MEM MAPQ ≥30 or AS ≥100; Bowtie 2 MAPQ ≥10 or AS ≥ −10. The mate of a deferred read is also deferred.
  - **Reconcile:** the lower edit distance wins, then the higher MAPQ, then a random choice.
  - **Chains:** built with nf-LO (minimap2 mode). The text names no file.

  Source: [Chen et al. 2024, Improved sequence mapping using a complete reference genome and lift-over (PMC11610747)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/)
- levioSAM2 accuracy with real Illumina 30x NovaSeq and BWA-MEM, against GIAB v4.2.1 (HG001, HG002, HG005):
  - **GATK HC, CHM13→GRCh38 vs direct GRCh38:** 18,308 fewer errors (−23.9%) and 11,551 fewer FPs (−24.9%). For CHM13→GRCh37 vs direct GRCh37: 49,770 fewer errors (−39.5%) and 38,981 fewer FPs (−52.2%).
  - **DeepVariant:** −3,443 errors (−11.4%) vs GRCh38; −8,778 (−16.7%) vs GRCh37. SNP FNs −22.3% vs GRCh38.
  - **HG002 CMRG:** GATK −730 errors (−35.2%) vs GRCh38; DeepVariant −1,010 (−20.0%) vs GRCh38.
  - **GATK in difficult strata (LowMap, SegDups, OtherDifficult):** −27.5% to −42.4% vs GRCh38.
  - **Simulation (chr20/chr21):** +0.08% and +1.37% correct mappings, where chr21 contains 771 kbp of false duplications.
  - **HiFi:** CHM13→GRCh38 had 673 more whole-genome errors (+0.8%) but 337 fewer in CMRG (−22.5%).
  - **Runtime, HG002 30x:** 47.6 CPU h vs 86.3 CPU h for direct GRCh38 (55.4%).

  Source: [Chen et al. 2024 (PMC11610747)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610747/); [bioRxiv 2022.04.27.489683](https://www.biorxiv.org/content/10.1101/2022.04.27.489683.full.pdf)
- Reference flow (Chen et al. 2021 Genome Biol) aligns to a set of population-major references (a global major plus 5 super-population majors from 1KGP) with Bowtie 2, picks the best alignment, and lifts to GRCh38. It has "a similar level of accuracy and bias avoidance" as vg "with 14% of the memory footprint and 5.5 times the speed" — [Reference flow, Genome Biology 2021 (PMC7780692)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7780692)
- HPRC draft pangenome (Liao et al. 2023 Nature): 47 diploid assemblies; 119 Mb of euchromatic polymorphic sequence relative to GRCh38 (about 90 Mb from SVs) and 1,115 gene duplications. Giraffe on the Minigraph-Cactus graph plus DeepVariant gave 21,700 errors (FP+FN) on HG005 30x, vs 36,144 for BWA-MEM/GRCh38 + DeepVariant and 26,852 for the Dragen pipeline. CMRG F1 was 0.985 for SNPs (others <0.976) and 0.961 for indels (others <0.958). The average reduction was "34% ... in false-positive and false-negative errors", and the largest gains were in "regions with errors in GRCh38 and large L1HS sequences". Read mappings are "projected from the pangenome space back onto an existing linear reference" — [Liao et al. 2023 (PMC10172123)](https://pmc.ncbi.nlm.nih.gov/articles/10172123)
- Illumina's own comparison used the Giraffe graph `hprc-v1.1-mc-grch38.d9` from the AWS human-pangenomics bucket. DRAGEN v4.3 reported 61.61% fewer errors than Giraffe-DeepVariant (vendor benchmark) — [Illumina, Second-gen multigenome](https://www.illumina.com/science/genomics-research/articles/second-gen-multigenome-mapping.html)
- Li (2014) already proposed mapping "to the ensemble of multiple genomes" as a direction, while noting that "a practical solution is yet to be concluded" — [Li 2014](https://arxiv.org/pdf/1404.0929)

### Inferences
- The user's rule (T2T confident placement >150 bp from the GRCh38 placement → MAPQ 0) is a stricter, flag-only variant of levioSAM2's suppress/defer/reconcile. levioSAM2 commits the CHM13 placement when it lifts, so it rescues reads rather than only flagging them. The user's 2x reduction in confident errors (scored at the locus) is larger than levioSAM2's 24–35% reduction in variant errors (scored per variant after calling). Both point to CHM13 as the source of the gain, consistent with "a second GRCh38 aligner does nothing".
- levioSAM2 is "one index" (CHM13) but not strictly one pass: deferred reads are re-mapped to GRCh38, so two indexes are needed. A one-pass CHM13-only variant would send deferred and suppressed reads to MAPQ 0 or unmapped in GRCh38 coordinates. This trades sensitivity in unliftable regions, which matters under the user's matched-sensitivity scoring.
- A Giraffe-on-HPRC graph (the GRCh38-referenced MC graph, `hprc-v1.1-mc-grch38`) followed by `vg surject` to GRCh38 is the other one-index, one-pass route that can carry the missing copies, if the graph's haplotypes include them. Whether CHM13 is one of its paths is believed but not verified here. Its surjected MAPQ in collapsed SDs is the key thing to measure.

### Gaps
- levioSAM2's per-sample F1 is in its supplementary tables only. The exact GRCh38 flavour used as the "direct" baseline (decoys/ALT or not) was not stated in the text read. levioSAM2 published no read-level "confident misplacement per million" metric.
- No source was found on whether `hprc-v1.1-mc-grch38` contains CHM13 as a path, on its index sizes, or on surjected-MAPQ behaviour in collapsed SDs. Not read: HPRC release notes, the vg wiki, and Sirén et al. 2021 Science (Giraffe).
- No published comparison of levioSAM2 (CHM13→GRCh38) against DRAGEN multigenome, or against GRCh38 plus CHM13 decoys, was found.

---

## 5. T2T-CHM13 as reference (chm13v2.0_maskedY_rCRS) and GRCh38 vs CHM13 mapping and false positives

### Takeaway
On 3,202 1KGP samples, CHM13 maps 0.97% more properly paired reads than GRCh38 with ALTs and decoys, with a 20–25% lower mismatch rate and about 4x more uniform coverage in genes. It removes the spurious AF≈0.5 PSV calls caused by collapsed SDs, and cuts FPs in 269 CMRG genes by up to 12x. The gains come from 8.04 Mb of collapsed loci (203 regions) and 1.2 Mb of false duplications that GRCh38 cannot fix in place. That makes CHM13 (or CHM13-derived sequence) the right donor for any GRCh38-coordinate augmentation aimed at the user's SD errors.

### Cited Findings
- **Assembly comparison:** CHM13v1.0 "adds or revises 229 Mbp" (non-syntenic over 1 Mb windows); 189 Mbp are "not covered by any primary alignments from GRCh38". Non-syntenic sequence is mostly centromeric satellites (190 Mbp) and SD copies (218 Mbp). Unique sequence increases by 14.9 Mbp (50-mers), 23.5 Mbp (100-mers) and 39.5 Mbp (300-mers). More than 106 Mbp of CHM13 sequence absent from GRCh38 "require reads of more than 300 bp to uniquely map" — [Aganezov et al. 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- **GRCh38 errors detected with haploid CHM13 data:** 368,574 heterozygous SNVs (15.3% shared between HiFi and 150 bp simulated Illumina) fell in 908 problematic regions (20.8 Mbp); 73.6% intersect SDs. gnomAD excess-heterozygosity sites are enriched there, accounting for 23.6% of the CHM13 heterozygous variants. 203 loci (8.04 Mbp) show missing copies in GRCh38, affecting 308 gene features, with 14 protein-coding genes fully inside (e.g. DUSP22, KMT2C). 12 false-duplication regions cover 1.2 Mbp — [Aganezov et al. 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- **1KGP short reads (BWA-MEM, 3,202 samples):** CHM13 maps "an additional 7.4 × 10^6 (0.97%) of properly paired reads ... compared to GRCh38, even when considering the alternative (ALT) and decoy sequences used in the original analysis". The per-read mismatch rate is 20–25% lower, there are fewer misoriented pairs, and the SD of coverage in genes is about 4x smaller. PASS variants per sample: mean 4,717,525 (CHM13) vs 5,101,897 (GRCh38). Mendelian-discordant variants in 602 trios: median 7,484 vs 8,879; low-quality variants: median 682,609 vs 890,701. The AF≈0.5 reduction "is largely explained by the corrections to collapsed SDs", i.e. "erroneous heterozygous variants ... caused by the false pileup of reads from the duplicated regions to a single location" — [Aganezov et al. 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- **Medically relevant genes (HG002 CMRG, 269 genes):** FP and FN counts "substantially decreased for all three call sets" (Illumina-BWA-MEM-GATK, HiFi, ONT). HiFi FPs fell 12x, "primarily because of the addition of missing sequences similar to KMT2C ... and removal of false duplications of CBS, CRYAA, H19, and KCNE1"; TPs decreased about 14%. In KCNE1, most variants are missed on GRCh38 "because many reads incorrectly map to a likely false duplication (KCNE1B on the p-arm of chromosome 21)". KCNJ18 sits in a GRCh38 collapsed duplication (17p11.2) with the missing paralog KCNJ17 (KCNJ18-1 in CHM13) — [Aganezov et al. 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- **Lifting CHM13 calls back to GRCh38** "requires special handling" when the two references carry different alleles, and representation differences in low-complexity regions make lifted sets differ from direct calls. ClinVar liftover GRCh38→CHM13 succeeded for 99.8% of variants (800,942 of 802,674) — [Aganezov et al. 2022](https://escholarship.org/content/qt18c391gh/qt18c391gh.pdf)
- **SweGen, 1,000 Swedish genomes, short reads:** mismatch rate 0.54% (CHM13) vs 0.65% (GRCh38) vs 0.78% (GRCh37). Proper pairing 99.4% vs 98.2% vs 97.9%. About 10 M more variants in the cohort (5.5 M singletons). The fraction of ambiguous alignments was higher on CHM13, reflecting SDs collapsed in GRCh37/38 — [Schmitz, Ameur, Johansson 2025 Genome Res (PMC12581843)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12581843/)
- **CHM13 v2.0 analysis-set files** (s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/, https mirror https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/):
  - `chm13v2.0.fa.gz`: soft-masked.
  - `chm13v2.0_noY.fa.gz`: identical to v1.1.
  - `chm13v2.0_maskedY.fa.gz`: chrY PARs hard-masked to N.
  - `chm13v2.0_maskedY_rCRS.fa.gz`: PARs masked and chrM replaced with rCRS NC_012920.1. The README spells it `chm13v2.0_maskedY.rCRS.fa.gz`. ChrY in v2.0 is from HG002.
  - `chm13v2.0_PAR.bed`.
  - Index files are on AWS, with md5s in `analysis_set/README.txt`.

  Source: [marbl/CHM13 GitHub](https://github.com/marbl/CHM13); [UCSC hs1 bigZips](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips)
- **Chains and synteny** (…/T2T/CHM13/assemblies/chain/v1_nflo/): `grch38-chm13v2.chain`, `chm13v2-grch38.chain`, `hg19-chm13v2.chain`, `chm13v2-hg19.chain`, `grch38-chm13v2.paf`, and `chm13v2-unique_to_hg38.bed` (regions of CHM13 non-syntenic with GRCh38). The synteny directory has `chm13_v1.1_plus38Y.no_snyteny_1Mbp.bed`. The short-read accessibility mask is `annotation/accessibility/combined_mask.bed.gz`. 1KGP recalls on CHM13v2.0 (3,202 or 2,504 samples) are under `variants/1000_Genomes_Project/chm13v2.0/` — [marbl/CHM13 GitHub](https://github.com/marbl/CHM13)

### Inferences
- `chm13v2-unique_to_hg38.bed` plus the Aganezov collapsed-duplication and false-duplication stratifications (GIAB genome-stratifications FTP, cited in the Aganezov data statement) give the coordinates needed to build "GRCh38 + CHM13 SD-copy decoys + GRC-exclusion masking" as a single bwa/minibwa index. The 190 Mbp of satellites should probably be excluded or kept as a single sink, because satellite decoys would absorb many reads without changing GRCh38 placements in the user's confident set.
- Given Aganezov's "0.97% more properly paired reads, even when considering ALT and decoy", hs38d1 and ALTs do not reproduce CHM13's mapping gain. This matches the user's observation that the gain is about CHM13-specific copies.
- The user's lift uses UCSC `hs1ToHg38.over.chain.gz`, while levioSAM2 and T2T use the nf-LO `chm13v2-grch38.chain`. Checking that the two chains give the same lifted loci in SDs is a cheap robustness test. (The UCSC liftOver path was not re-verified here.)

### Gaps
- File sizes for the CHM13 analysis-set FASTA and indexes were not stated on the pages read (no size given).
- No study was found that maps GRCh38-truth-projected real reads and reports read-level misplacement per million for CHM13 vs GRCh38, which would be the direct analogue of the user's metric. Aganezov, levioSAM2 and SweGen report mapping summary statistics or variant metrics.
- Aganezov's Table S15 (per-call-set CMRG FP/FN on GRCh38 vs CHM13 for Illumina) was not read. Only the "up to 12x" (HiFi) figure is in the main text.
