# Components of the WES workflows

Everything the two workflows use: programs, Galaxy tools, reference and knowledge data, and the analyses built on
them. Update this file with every release. `install.sh` installs all of it (column *Stage*: the `install.sh` stage
that provides it); README section 1 explains the stages, section 2 the workflow inputs, section 7 the validation.

Release 0.15 (2026-11-21). GRCh38 throughout, Ensembl chromosome names (`1`, `MT`) unless a file says otherwise.

## 1. Analyses, by release

| Release | Date | Added | Validation (README section 7) |
|---|---|---|---|
| FASTQ v1 / VCF v2 | 2026-11-07 | QC (FastQC, fastp, MultiQC), BWA-MEM alignment, duplicate marking, FreeBayes calling, normalisation, SnpEff annotation, ClinVar and gnomAD frequencies, filtering, coverage (mosdepth, HS metrics) | HG001 end-to-end run |
| report 0.1 | 2026-11-16 | Clinical interpretation report (HTML + findings TSV), runs of homozygosity, complex-allele atomisation | Clinical report on HG001 |
| report 0.2 | 2026-11-16 | gnomAD v4.1 joint frequencies (FAF95, grpmax), ACMG secondary findings v3.3 | |
| report 0.3 | 2026-11-17 | Recurrent-artefact list of the laboratory | |
| report 0.4 | 2026-11-17 | AlphaMissense, SpliceAI (local tool `spliceai_score`), PharmCAT (local tool `pharmcat_pgx`), GenCC new export | Recent GIAB exomes |
| report 0.5 | 2026-11-18 | ClinGen dosage, ClinVar residue evidence (PS1/PM5), HPO disease match, read-backed phase | |
| release 0.6 | 2026-11-19 | Trio analysis: de novo, parental origin, compound heterozygotes, parentage | GIAB trios |
| release 0.7 | 2026-11-19 | Exome CNV calling (local tool `wes_cnv`) and CNV interpretation (ACMG/ClinGen CNV points); `install.sh` | GIAB IDT exomes |
| report 0.8 | 2026-11-20 | Automated ACMG/AMP classification | ClinVar ≥ 2 stars (30,000 variants) |
| report 0.9 | 2026-11-20 | Uniparental disomy and parent of origin from trios | GIAB trios |
| report 0.10 | 2026-11-20 | Reanalysis: changes since a previous report | GIAB exomes |
| release 0.11 | 2026-11-20 | MANE naming, curated gene-disease table, ClinGen VCEP rules and Evidence Repository, REVEL, MaveDB PS3/BS3, Collins dosage, gnomAD exome CNVs, UniProt features, provenance; mitochondrial DNA (local tool `wes_mito`) | Knowledge sources; mtDNA |
| release 0.12 | 2026-11-21 | gnomAD exome coverage gates PM2; pext next to loss-of-function variants | Coverage and pext |
| release 0.13 | 2026-11-21 | Couple carrier screening (local tool `wes_couple`), carrier gene table; sex check in the report | Couple screening |
| release 0.14 | 2026-11-21 | HLA typing (local tool `wes_hla`, T1K): class I/II alleles, CPIC HLA pharmacogenomics through PharmCAT outside calls, HLA disease associations | HLA typing |
| release 0.15 | 2026-11-21 | mtDNA haplogroup (`wes_mito` 0.2.0, HaploGrep3) and Y-chromosome haplogroup (local tool `wes_yhaplo`, Yleaf) as maternal/paternal lineages in the report, with a Y-versus-sex check | Haplogroups |

## 2. Programs on the host

Used by the data scripts (`data/`), not by Galaxy jobs (those get their own Conda environments).

| Program | Version (tested) | Licence | Stage |
|---|---|---|---|
| bcftools, samtools, tabix/bgzip (htslib) | 1.21 (apt) | MIT/Expat | `system` (apt, or Conda with `--conda`) |
| bedtools | 2.31.1 | MIT | `system` |
| curl, gawk, python3 (+ venv) | distribution | various | `system` |
| Ephemeris, BioBlend (Python) | 0.10.11, 1.9.0 | AFL-3.0, MIT | `venv` (`~/wes-venv`) |
| gxformat2 (Python, only to rebuild the `.ga` files with `build.py`) | 0.29.2 | AFL-3.0 | a throwaway venv |

## 3. Galaxy tools

Tool Shed repositories (`tools.yml`, stage `tools`; latest installable revisions):

| Tool | Licence | Used for |
|---|---|---|
| FastQC (devteam), fastp (iuc), MultiQC (iuc) | GPL-3.0, MIT, GPL-3.0 | read QC and the QC report |
| BWA (devteam), Picard (devteam), samtools stats / faidx | GPL-3.0, MIT, MIT | alignment, duplicates, HS metrics, statistics |
| mosdepth (iuc) | MIT | per-target and per-exon coverage, CNV depth |
| FreeBayes (devteam) | MIT | variant calling |
| bcftools norm / filter / annotate / stats / roh / view, bedtools (iuc) | MIT | normalisation, annotation, filtering, ROH, regions |
| SnpEff, SnpSift (iuc) | MIT | consequences and HGVS; prioritisation filter |

Local tools (`tools/`, stage `local-tools`; Galaxy installs their Conda requirements):

| Tool | Version | Conda requirements | Licence of the main dependency | Added |
|---|---|---|---|---|
| `wes_clinical_report` | 0.15.0 | python 3.12, samtools 1.22 | this repository | 0.1 |
| `spliceai_score` | 1.3.1 | spliceai 1.3.1, numpy 1.26.4, setuptools 80.9, bcftools 1.24 | SpliceAI: GPL-3.0 | 0.4 |
| `pharmcat_pgx` | 3.4.0+galaxy1 | pharmcat3 3.4.0 | MPL-2.0 (guidelines: CPIC CC0, DPWG, FDA via ClinPGx) | 0.4 |
| `wes_cnv` | 0.1.0 | python 3.12, numpy 2.2.6, samtools 1.22 | this repository | 0.7 |
| `wes_mito` | 0.2.0 | python 3.12, samtools 1.22, haplogrep3 3.2.2 | HaploGrep3: MIT | 0.11 (haplogroup 0.15) |
| `wes_couple` | 0.1.0 | python 3.12 | this repository | 0.13 |
| `wes_hla` | 0.1.0 | t1k 1.0.10, samtools 1.22, python 3.12 | T1K: MIT | 0.14 |
| `wes_yhaplo` | 0.1.0 | yleaf 4.1.7, samtools 1.22, python 3.12 | Yleaf: GPL-3.0 | 0.15 |

Tested on aarch64 (Ubuntu, Galaxy 26.2); every requirement also has a bioconda/conda-forge build for x86-64.
SpliceAI's TensorFlow 2.19 is the only aarch64 build and needs numpy 1.26.

Data inside the Conda packages (no separate download): HaploGrep3's mtDNA tree PhyloTree 17 Forensic Update 1.2
(6,380 haplogroups; installed by the package's post-link step from GitHub, so the Galaxy dependency install needs
internet) and Yleaf's YFull v14 and ISOGG trees with their GRCh38 marker positions. HaploGrep3's FASTA input needs an
x86-64 native library; the tools pass HSD text, which works on aarch64.

## 4. Reference and knowledge data

Stage `data` downloads and builds them into `--data-dir` and lists each file with its source and date in
`MANIFEST.tsv`; stage `upload` puts them into one Galaxy history.

| File | Source | Licence | Used for | Added |
|---|---|---|---|---|
| `GRCh38.ensembl115.primary_assembly.fa` | Ensembl 115 | no restrictions | reference genome | v1 |
| `canonical_cds.GRCh38.115.bed` | Ensembl 115 GTF | no restrictions | per-gene coverage | v1 |
| SnpEff `GRCh38.115` | SnpEff download tool (stage `upload`) | MIT / Ensembl | consequences | v1 |
| `clinvar.GRCh38.vcf.gz` | NCBI ClinVar (weekly) | public domain | ClinVar classifications | v1 |
| `gnomad.v4.1.joint.targets.GRCh38.bcf` | gnomAD v4.1 joint sites, `data/build_gnomad_extract.sh` (needs `--regions`) | CC0 | population frequencies | 0.2 |
| `gnomad.v4.1.constraint_metrics.tsv` | gnomAD v4.1 | CC0 | LOEUF / pLI | v1 |
| `gencc-submissions.tsv` | GenCC | CC0 | gene-disease validity | v1 |
| `phenotype_to_genes.txt`, `phenotype.hpoa` | HPO release | HPO licence (free, attribution) | phenotype ranking, disease match | v1, 0.5 |
| `ClinGen_gene_curation_list_GRCh38.tsv`, `ClinGen_region_curation_list_GRCh38.tsv` | ClinGen FTP | ClinGen terms (free) | dosage sensitivity, recurrent CNV regions | 0.5, 0.7 |
| `clinvar_protein.tsv.gz` | ClinVar `variant_summary`, `data/build_clinvar_protein.py` | public domain | PS1 / PM5 | 0.5 |
| `alphamissense.hg38.bcf` | Zenodo 10813168, `data/build_alphamissense.sh` | CC BY 4.0 | missense predictions | 0.4 |
| `missense_scores.hg38.bcf` | AlphaMissense + REVEL v1.3 (Zenodo 7072866), `data/build_revel.sh` | CC BY 4.0 + REVEL non-commercial | PP3 / BP4 (REVEL first) | 0.11 |
| `GRCh38_segdups.bed` | GIAB stratifications v3.6 | public domain | CNV confidence | 0.7 |
| `MANE.GRCh38.v1.5.summary.txt.gz` | NCBI MANE 1.5 | public domain | MANE naming | 0.11 |
| `gene_disease.tsv` | ClinGen validity, G2P, Orphanet; `data/build_gene_disease.py` | CC0, CC0, CC BY 4.0 | disease context, inheritance, mechanism | 0.11 |
| `vcep_rules.tsv`, `erepo_variants.tsv` | ClinGen CSpec registry, Evidence Repository; `data/build_vcep.py` | ClinGen terms | gene-specific ACMG rules, expert classifications | 0.11 |
| `Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz` | Zenodo 6347673 | CC BY 4.0 | pHaplo / pTriplo | 0.11 |
| `mave_evidence.tsv.gz` | MaveDB bulk download (Zenodo 20840937), `data/build_mavedb.py` | CC0 | PS3 / BS3 | 0.11 |
| `gnomad.v4.1.cnv.all.vcf.gz` | gnomAD v4.1 exome CNVs | CC0 | CNV frequencies (4O) | 0.11 |
| `uniprot_features.tsv.gz` | UniProtKB/Swiss-Prot, `data/build_uniprot_features.py` | CC BY 4.0 | protein features | 0.11 |
| `gnomad.v3.1.chrM.vcf.gz` | gnomAD v3.1 chrM | CC0 | mtDNA frequencies | 0.11 |
| `gnomad_exome_coverage.tsv`, `gnomad_pext.tsv` | gnomAD v4.0 coverage, GTEx v10 pext; `data/build_gnomad_sites.py` (needs `--regions`) | CC0 | PM2 gate, pext | 0.12 |
| `carrier_genes.tsv` | ClinVar × gnomAD v4.1 × gene-disease table × ClinGen actionability; `data/build_carrier_genes.py` (needs `--regions`) | as its sources | couple screening | 0.13 |
| `hla.dat` | IPD-IMGT/HLA (current release; the file states it, e.g. 3.66.0) | CC BY-ND 4.0 | HLA typing | 0.14 |

Not used, and why: OMIM (licence), HGMD (commercial), dbNSFP (registration pending), MITOMAP (no bulk download),
COSMIC/OncoKB (somatic, commercial licences).
