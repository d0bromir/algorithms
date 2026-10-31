# Large-scale truth for read placement (recent data only)

Where a read *came from* is known in two ways: it was simulated from a known
genome, or it is a real read whose origin can be read off the sample's own,
near-perfect diploid assembly. Both now exist at whole-genome scale for
recent Illumina instruments. This note lists the resources (checked
2026-10-29; nothing downloaded yet) and how to use them to test `--bio`
(maximum-likelihood placement, [THEORY.md §12](../THEORY.md)) against the
default mode and the other aligners.

## Resources

| Resource | What it gives | Where |
|---|---|---|
| **T2T-HG002 Q100 assembly** v1.1 (July 2024, Merqury QV 77.8) and v1.2 (August 2026), both haplotypes (`chr1_MATERNAL`, `chr1_PATERNAL`, …) | the sample's own genome: true origin of every HG002 read, including segmental duplications, centromeres, satellites | `s3://human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz`, `…/hg002v1.2.fasta.gz` ([github.com/marbl/HG002](https://github.com/marbl/HG002), [Q100 project](https://genomeinformatics.github.io/HG002v1)); haplotype-to-haplotype chain `…/assemblies/changes/hg002v1.2_to_other_haplotype.chain.gz` |
| **GIAB HG002 v5.0q benchmark** (small variants and SVs, from the v1.1 assembly) | downstream truth: 701 Mb (11.7 %) more autosomal sequence than v4.2.1, error-free on 99.4 % of the diploid genome per the release | [NIST GIAB](https://www.nist.gov/programs-projects/genome-bottle); [Cell paper, PMC12458380](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458380) |
| **HG002 NovaSeq X 25B** reads (already used here: calibration) | real recent reads, 2/9/24/40 quality bins, 95 % of bases in Q40 | SRA SRR37356339 ([bench/RESULTS.md](RESULTS.md)) |
| **HG005 NovaSeq 6000** 2×151 PCR-free (already used: held out) | a second sample and instrument | SRA SRR14724528 |
| **HPRC Data Release 2** (May 2025): phased assemblies of 232 individuals, with Illumina WGS for over 200 | the same personal-genome truth across ancestries | [humanpangenome.org/hprc-data-release-2](https://humanpangenome.org/hprc-data-release-2/), [data.humanpangenome.org](https://data.humanpangenome.org/), BioProject PRJNA730822, `s3://human-pangenomics` |
| **Platinum Pedigree** (CEPH-1463, 28 members, four generations; Nature Methods 2025): phased assemblies and a Mendelian-consistent truth set (4.7 M SNVs, 768 k indels, 537 k tandem repeats, 24 k SVs over 2.77 Gb of GRCh38) | inheritance-checked truth; five platforms | `s3://platinum-pedigree-data` (us-west-1, CC BY 4.0, [registry](https://registry.opendata.aws/platinum-pedigree/)); the Illumina instrument and year are not stated there and must be checked before use |
| **NovaSeq X quality bins** (software v1.2): 2–17, 18–29, ≥ 30, the top bin reported as Q40 | the simulation profile | [Illumina, data quality](https://www.illumina.com/science/genomics-research/articles/data-quality-q-scores.html) |

## Three tests

**A. Simulated reads from the HG002 assembly (exact truth, any region).**
Simulate NovaSeq X-profile reads and pairs from both haplotypes
(`tools/simulate.py reads … --quality novaseqx --mates …`: the bins and
3′ decline measured on SRR37356339; sample variants are already in the
haplotypes, so `--snv`/`--indel` stay small, for somatic-like noise). Each
read's origin on its haplotype is exact. Project it to GRCh38 (and
T2T-CHM13) through each haplotype's whole-genome alignment to the reference
(`minimap2 -x asm5`, or the alignments GIAB's dipcall run made for v5.0q).
Measure, per GIAB stratification (segmental duplications, low mappability,
tandem repeats, CMRG genes):
- the fraction of reads in SAM and how many are at the projected truth;
- MAPQ calibration: misplaced reads against Σ 10^(−MAPQ/10)
  (`simulate.py eval` prints it);
- for `--bio` against the default mode: reads moved to the truth, reads
  moved away from it.

**B. Real HG002 NovaSeq X reads, truth from the sample's own genome.**
Map the real reads to each haplotype of the HG002 assembly. A read with a
unique certified placement (CERTA's edit certificate with the haplotype as
the reference: unique best, no second locus within the radius) that is
better than any placement on the other haplotype, or the same homologous
placement on both, has a known origin. Project the origin to GRCh38 as in A. This gives a truth placement for
most real reads, in difficult regions too, and every aligner's GRCh38
placement can be scored against it: `certa --bio`, `certa` default,
BWA-MEM2, minibwa, DRAGMAP, strobealign, minimap2 (the set in
`bench/final_bench.sh`). Reads whose origin falls where the haplotype does
not align to GRCh38 (SVs, sequence missing from GRCh38) are reported
separately: no linear-reference aligner can place them correctly.

**C. Downstream.** Variant calls from each aligner's BAM (DeepVariant or
GATK) against v5.0q, stratified, with hap.py/vcfeval: the placement
differences of A and B as variant-level errors.

Repeat B (and C where a benchmark exists) on HG005 (NovaSeq 6000, held out)
and on HPRC Release 2 samples with Illumina WGS of recent instruments, chosen
across ancestries, to show the result is not one sample's.

## Scale and cost

A 30–40× HG002 run is 600–900 M reads. The diploid assembly is about twice
GRCh38, above CERTA's 32-bit positions (a reference must be shorter than
2^32 bases), so each haplotype gets its own index (about GRCh38's size each)
and a read is mapped to both: a unique certified placement on one, or the
same homologous placement on both, fixes its origin. Test A needs no
download beyond the assembly and GRCh38; test B also the NovaSeq X run
already on galaxy. Runs on galaxy and a2 follow the host coordination with
the paper session.
