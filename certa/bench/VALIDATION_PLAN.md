# Validating `--bio` on real data: plan

Status, 2026-10-31: the pilot and phase 1 (tests A and B on HG002) are done;
the results are in [VALIDATION_RESULTS.md](VALIDATION_RESULTS.md). Phases 2
and 3 have not started. The resources are described in
[TRUTH_DATA.md](TRUTH_DATA.md); this note fixes the files, the definitions of
truth, the code and the order of the runs.

## Claims and how each is tested

1. **Soundness at scale.** No aligner finds an alignment that the `--bio`
   model scores below the certified `YP` of that read. The expected count of
   violations is 0, over tens of millions of real reads.
2. **Accuracy.** On reads whose origin is known, `--bio` places more reads
   at the origin than CERTA's default mode and the other aligners,
   especially in difficult regions. Comparisons are made at matched numbers
   of confident reads.
3. **Calibration.** In each MAPQ band, the fraction of misplaced reads
   matches 10^(−MAPQ/10).

**Freeze rule.** The build and the `--bio` defaults are committed and tagged
before the first real run. HG002 is the calibration sample: its results may
motivate changes. HG005 is held out and is run once with the frozen build.
Its result is the verdict.

## Already on the hosts (checked 2026-10-31)

| | galaxy | a2 |
|---|---|---|
| State | idle apart from resident ollama servers; 128 cores, 246 GiB RAM, 2 × A100 | idle; 64 cores, 376 GiB RAM |
| Free disk | 232 GB | 66 GB (97 % full) |
| HG002 NovaSeq X (SRR37356339; 444.5 M pairs, 44×, public 2026-03) | 20 M and 2 M pairs, uncompressed | the same 20 M / 2 M pairs; full mate 2 (28.8 GB gz) |
| HG005 NovaSeq 6000 (SRR14724528, 362.9 M pairs) | full run (51.9 GB gz); 2 M pairs | — |
| GRCh38 | Ensembl primary assembly (194 contigs named `1`…`22`, `X`, `Y`, `MT`, `KI270728.1`, …); indexes for CERTA, BWA-MEM2, minibwa, DRAGMAP, strobealign, minimap2, bowtie2 | the same, except DRAGMAP |
| Tools | minimap2, bgzip/tabix, pigz, python3, Java 21, apptainer, micromamba | minimap2, k8, python3, micromamba, sra-tools |

Neither host has samtools, bcftools, rtg-tools or hap.py. galaxy also holds
GIAB v4.2.1 for HG001/HG002/HG005; that release is legacy and is not used
for HG002. The full mate 1 of SRR37356339 was not kept. a2 is not needed
except to copy mate 2 in phase 3.

## Downloads (each phase needs approval)

**Phase 1: HG002 truth, ≈ 5.0 GB to galaxy `~/certa-data/truth/`**

| File | Size | Source | Used for |
|---|---|---|---|
| `hg002v1.1.mat_Y_EBV_MT.fasta.gz` | 0.91 GB | `s3://human-pangenomics/T2T/HG002/assemblies/` | maternal haplotype plus chrY: one genome, one index |
| `hg002v1.1.pat_X_EBV_MT.fasta.gz` | 0.91 GB | same | paternal haplotype plus chrX |
| `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap1.bam` (+ `.bai`) | 1.04 GB | [GIAB `latest/dipcall_output/`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/dipcall_output/) | haplotype → GRCh38 projection |
| `GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap2.bam` (+ `.bai`) | 1.10 GB | same | haplotype → GRCh38 projection |
| `HG002_GRCh38_v5.0q_smvar.{vcf.gz,vcf.gz.tbi,benchmark.bed}`, `…_stvar.{vcf.gz,benchmark.bed}`, README | 92 MB | [GIAB `latest/`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/) | benchmark regions; variants at the origin |
| v3.6 GRCh38 stratifications: `segdups`, `lowmappabilityall`, `AllTandemRepeatsandHomopolymers_slop5`, `alldifficultregions`, `alllowmapandsegdupregions`, `notinalldifficultregions`, `MHC`, `CMRGv1.00_falselyduplicatedgenes` | ≈ 75 MB | [GIAB v3.6 `GRCh38@all/`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/) | strata |
| samtools, bcftools (bioconda, via micromamba) | ≈ 50 MB | bioconda | reading the dipcall BAMs; CRAM later |

The plan uses assembly v1.1, not v1.2. The v5.0q benchmark and the
dipcall alignments were built from v1.1, so the assembly, the projection
and the benchmark stay consistent.

**Phase 2: held-out sample, 1.8 GB.**
`HG005_pat_hprc_r2_v1.0.1.fa.gz` (0.87 GB) and `HG005_mat_hprc_r2_v1.0.1.fa.gz`
(0.90 GB) come from HPRC Release 2 (trio-phased hifiasm, December 2024) under
`s3://human-pangenomics/working/HPRC_PLUS/HG005/assemblies/release2/`. Their
GRCh38 alignment is computed on galaxy with minimap2, using dipcall's settings
(`-x asm5 -z200000,10000`), so it needs no download. The reads are
already there.

**Phase 3: full scale and variant calls, ≈ 30 GB plus tools.**
`SRR37356339_1.fastq.gz` (29.4 GB) comes from ENA; mate 2 is copied from
a2. The variant caller and the comparison tool are chosen after phase 2
(see Test C).

## How truth is defined

**Test A: simulated reads.**
- Reads are simulated from both haplotypes: 5 M pairs each, with
  `tools/simulate.py reads … --quality novaseqx --mates`. The sample's
  variants are already in the haplotypes, so `--snv`/`--indel` stay at
  somatic-like levels.
- The origin on the haplotype is exact. It is projected to GRCh38 through
  that haplotype's dipcall alignment, the alignment v5.0q was built from.
- A read has a GRCh38 truth when both ends fall in one chain of aligned
  blocks, in a consistent orientation, and the projected span differs from
  the read's by at most the gap cap (50).
- Otherwise the read is labelled **no linear truth**: an SV, sequence missing
  from GRCh38, or a break in the alignment. It is reported separately.
- **Correct** means: right contig and strand, leftmost aligned base within
  20 bp of the projected one (wgsim_eval's convention). Exact agreement is
  reported too. Mates are scored separately.

**Test B: real reads, truth from the sample's own genome.**
1. Adapters and poly-G tails are trimmed with fastp-gpu, and every aligner
   gets the same trimmed FASTQ. `--bio` aligns end to end, so untrimmed
   reads would fall back.
2. The pairs are mapped to each haplotype with CERTA's default edit
   certificate. Each haplotype is 3.15–3.17 Gbp, below 2^32, so it fits in
   its own index.
3. The origin is known when the pair is certified unique on a haplotype,
   with no second placement within *m* edits, and either:
   - it beats every placement on the other haplotype by *m*; or
   - both haplotypes give certified unique placements that project to the
     same GRCh38 locus (homologous sequence).

   The primary run uses *m* = 2, with *m* = 1 and 3 as a sensitivity
   check.
4. The origin is projected as in Test A.

*Independence.* The labels are checked with a second aligner on the
haplotypes: minibwa, BWA-MEM's algorithm, whose index builds in 5 minutes
rather than BWA-MEM2's 25. Reads where its placement disagrees are excluded
and counted, so CERTA is not scored against labels only CERTA produced.

*Bias.* The labelled set leaves out reads that are ambiguous even within
the sample's own genome; no aligner can be scored on them. Their fraction
and each aligner's MAPQ distribution on them are reported.

**Test C: downstream variant calls.** It needs full coverage (20 M pairs is
about 2×), so it runs in phase 3. galaxy is ARM64:
- DeepVariant has no ARM64 image.
- GATK HaplotypeCaller and rtg-tools vcfeval are Java and run there.
- hap.py is x86-64 only and would run on a2.

The tools are chosen after phase 2.

## Measures

Each aligner is measured in each stratum: the v5.0q benchmark regions,
`notinalldifficult`, segmental duplications, low mappability, tandem repeats
and homopolymers, MHC, the CMRG falsely duplicated genes, and no linear
truth.

- **Placement:** the fraction in the SAM file, the fraction at truth, and
  the misplaced reads among those with MAPQ ≥ 1, 20, 40 and 60. A curve
  plots misplaced against placed reads as the MAPQ threshold falls.
- **Calibration:** in each MAPQ band, observed misplaced reads against
  Σ 10^(−MAPQ/10), as `simulate.py eval` prints it.
- **`--bio` against the default mode:**
  - reads moved to the truth, moved away from it, or wrong in both;
  - offsets inside tandem repeats (`YA > 0`) and ties (`YB`), counted
    separately.
- **Foreign reads:** reads whose origin is not in GRCh38 that an aligner
  still places with MAPQ ≥ 20.
- **Soundness (claim 1):**
  - Every aligner's end-to-end alignment of a read is rescored under the
    `--bio` model. A cost below `YP`, for a read `--bio` certified, is a
    violation; the expected count is 0.
  - Soft-clipped alignments are skipped, since the model is end to end.

The aligners are `certa --bio`, `certa` default, BWA-MEM2, minibwa, DRAGMAP,
strobealign and minimap2: the set in [final_bench.sh](final_bench.sh), run
on the same input with the same threads.

## Code to write

| File | Does |
|---|---|
| `bench/truth_lift.py` | builds the block map from the dipcall BAMs (via `samtools view`) and projects positions; translates contig names from Ensembl to UCSC (`1` ↔ `chr1`, `KI270728.1` ↔ `chrUn_KI270728v1`); flags the PARs |
| `bench/truth_label.py` | derives the Test B labels from the per-haplotype SAM files |
| `bench/truth_eval.py` | computes the measures and tables per aligner and stratum (generalizes `simulate.py eval`) |
| `tools/bio_score.cpp` | computes the `--bio` model cost of any SAM record (`bio_rescore` from `bio.h`), for the soundness check |
| `bench/truth_run.sh` | runs the pipeline on galaxy |

Before any result is used:
- **Projection:** error-free simulated reads from variant-free blocks must
  equal GRCh38 at their projected position, for 100 % of reads.
- **Rescoring:** `bio_score` must reproduce `YP` on `--bio`'s own records,
  for 100 % of records.
- **Disk:** SAM output is reduced to per-read TSV while it streams, so
  phases 1–2 stay below about 90 GB on galaxy.

## Order of runs (galaxy, `taskset -c 0-63`)

0. **Pilot (no download).**
   - Build `main` on galaxy and run `--bio` on the existing 2 M HG002 pairs
     against GRCh38.
   - Measure throughput and the fallback share. Estimate the insert size
     from the default mode's pairs; it sets `--bio-insert`.
   - Throughput is the unknown. The default mode maps 20 M reads in 44 s
     on 64 cores. `--bio` cost 2.2–6.7× the default on a 20 Mbp
     simulated genome, and human repeats will cost more.
1. **Phase 1:** index the two haplotypes (about 10 GB each), then Test A
   (10 M pairs) and Test B on the 20 M HG002 pairs.
2. **Phase 2:** Test B on 20 M HG005 pairs from SRR14724528, with the
   frozen build.
3. **Phase 3:** Test B on the full HG002 run (444.5 M pairs), then Test C.
   Storing alignments needs CRAM to fit the disk.

The hosts are shared with the paper session. Before each phase, ask that
session for galaxy (SendMessage), and send "free" afterwards.

## Later, optional

- **More HPRC Release 2 samples.** 231 have Illumina WGS, all NovaSeq 6000
  2×150. Check the run dates against the recent-data rule before use.
- **Platinum Pedigree.** Illumina data are in
  `s3://platinum-pedigree-data/data/illumina/`. The instrument isn't stated
  there; check it before use.
