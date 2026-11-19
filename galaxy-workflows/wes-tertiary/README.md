# Whole exome sequencing (WES): tertiary analysis workflows for Galaxy

Two Galaxy workflows for germline whole exome sequencing. Import the `.ga` files.

| File | Version | Starts from | What it does |
|---|---|---|---|
| [`wes-tertiary-from-fastq.ga`](wes-tertiary-from-fastq.ga) | **v1** | paired-end FASTQ (`list:paired` collection, one element per sample) | Read QC and trimming, alignment, duplicate marking, capture/coverage QC, variant calling, then the full v2 tertiary analysis, plus one MultiQC report for all samples |
| [`wes-tertiary-from-vcf.ga`](wes-tertiary-from-vcf.ga) | **v2** | a germline VCF/BCF from any caller (FreeBayes, GATK, DeepVariant, DRAGEN, …) | Normalisation, restriction to the exome targets, quality filtering, ClinVar, gnomAD and SnpEff annotation, prioritisation, tables, runs of homozygosity and a **clinical interpretation report** |

Both end in a human-readable HTML report of the medically significant findings (section 4a): diagnostic
findings, candidate variants, ACMG secondary findings, carrier status, pharmacogenomics, risk alleles,
runs of homozygosity, coverage of clinically relevant genes, and what the analysis cannot detect.

v1 embeds v2 as a subworkflow, so both workflows apply exactly the same tertiary logic.

```mermaid
flowchart LR
  subgraph v1["v1: from FASTQ (per sample)"]
    fq[FASTQ pairs] --> fastqc[FastQC] & fastp[fastp trimming]
    fastp --> bwa[BWA-MEM] --> md[Picard MarkDuplicates]
    md --> stats[samtools stats] & hs[Picard CollectHsMetrics] & mos[mosdepth per-target coverage] & cds[mosdepth coding-exon coverage]
    md --> fb[FreeBayes]
  end
  subgraph v2["v2: from VCF (tertiary analysis)"]
    vcf[VCF] --> norm[bcftools norm<br/>split multiallelics, atomize complex alleles,<br/>left-align, REF check]
    norm --> filt[bcftools filter<br/>padded exome targets + quality]
    filt --> clin[bcftools annotate<br/>ClinVar] --> gnomad[bcftools annotate<br/>gnomAD v4.1 AF, grpmax, FAF95,<br/>homozygotes as gnomAD_*]
    gnomad --> am[bcftools annotate<br/>AlphaMissense AM_score] --> sai[SpliceAI<br/>rare variants]
    sai --> eff[SnpEff<br/>consequence, HGVS, LOF/NMD]
    eff --> prio[SnpSift filter<br/>prioritisation] --> ctab[candidate table]
    eff --> atab[annotated table] & bstats[bcftools stats] & roh[bcftools roh]
    vcf --> pcat[PharmCAT<br/>diplotypes, CPIC/DPWG report]
    eff & roh & pcat --> rep[clinical report<br/>HTML + findings TSV]
  end
  fb --> vcf
  md & mos & cds -.-> rep
  fastqc & fastp & md & stats & hs & bstats & eff --> mqc[MultiQC]
```

## 1. Install the tools and import the workflows

**Quick install.** [`install.sh`](install.sh) runs every step of this section and section 2 on the Galaxy server.
It installs the command-line programs, the Tool Shed and local tools with their dependencies, and the workflows.
It also downloads and builds the reference and knowledge files and uploads them into a history:
```bash
./install.sh check
```
```bash
./install.sh all --regions kit.bed
```
| Stage | What it does |
|---|---|
| `check` | Reports what is installed and what is missing; changes nothing |
| `system` | Installs bcftools, bedtools, samtools, tabix, curl, gawk and python3-venv with apt (`--conda`: into a conda environment instead) |
| `venv` | Creates a Python environment with Ephemeris and BioBlend (`$VENV`, default `~/wes-venv`) |
| `tools` | Installs the [`tools.yml`](tools.yml) repositories with Ephemeris; already-installed repositories are skipped |
| `local-tools` | Copies the three local tools to the local tool directory and checks `local_tool_conf.xml`, `galaxy.yml` and the HTML allowlist. Then Galaxy installs the tools' Conda dependencies |
| `workflows` | Imports both `.ga` files; with `--update`, adds a new version to workflows already present |
| `data` | Downloads and builds the section 2 files into `--data-dir`, with a `MANIFEST.tsv` (source, date, size) |
| `upload` | Uploads the `MANIFEST.tsv` files into a new history and runs *SnpEff download* there |

- **API key:** read from `$GALAXY_API_KEY` or `~/.galaxy_admin_key`, and passed to the Python helper ([`install/galaxy_setup.py`](install/galaxy_setup.py)) in its environment, not on a command line.
- **Galaxy configuration:** changed only with `--configure`, which writes a missing `local_tool_conf.xml` and adds the report tools to the allowlist. Otherwise the script prints what to add, since Galaxy servers managed by Ansible should change it in their templates.
- **gnomAD extract:** built only with `--regions` (the capture kit BEDs), because it reads about 2 hours from the gnomAD bucket. The capture kit BED itself is the laboratory's own and is not downloaded.
- **Re-running:** every stage can be run again; what is already in place is left alone. On the test server (November 2026 releases), the `data` stage without gnomAD took 16 min, most of it the Ensembl FASTA at 1.5 MB/s. It wrote 3.8 GB, with a peak memory of 2.1 GB (the AlphaMissense build).

The manual steps follow.

1. **Tools**: as a Galaxy admin, install the repositories listed in [`tools.yml`](tools.yml). You can use
   *Admin → Install and Uninstall*, or [Ephemeris](https://ephemeris.readthedocs.io). Install Ephemeris into a
   virtual environment; Debian/Ubuntu refuse a system-wide `pip install` with "externally-managed-environment":
   ```bash
   python3 -m venv ~/ephemeris-venv && ~/ephemeris-venv/bin/pip install ephemeris
   ```
   ```bash
   ~/ephemeris-venv/bin/shed-tools install -g http://192.168.191.12 -a <ADMIN_API_KEY> -t tools.yml
   ```
   If `python3 -m venv` reports that `ensurepip` is not available, install `python3-venv` (or `python3-full`)
   with `apt` first. The API key belongs to a Galaxy admin user (*User → Preferences → Manage API Key*).
2. **Local tools** (not in the Tool Shed): copy [`tools/wes_clinical_report`](tools/wes_clinical_report),
   [`tools/spliceai_score`](tools/spliceai_score) and [`tools/pharmcat_pgx`](tools/pharmcat_pgx) to the server,
   list them in a local tool config, and let Galaxy show the report HTML unsanitised:
   ```xml
   <!-- /etc/galaxy/local_tool_conf.xml -->
   <toolbox tool_path="/srv/galaxy/local_tools">
       <section id="wes_analysis_local" name="WES analysis">
           <tool file="wes_clinical_report/wes_clinical_report.xml"/>
           <tool file="spliceai_score/spliceai_score.xml"/>
           <tool file="pharmcat_pgx/pharmcat_pgx.xml"/>
       </section>
   </toolbox>
   ```
   Add that file to `tool_config_file` in `galaxy.yml`, set `sanitize_allowlist_file` to a file containing the
   line `wes_clinical_report` (or use *Admin → Manage Allowlist*), restart Galaxy, and install the tools'
   dependencies from *Admin → Manage Dependencies*: Python 3.12 and samtools 1.22 (report), SpliceAI 1.3.1 with
   setuptools 80.9.0 and bcftools 1.24 in one environment (TensorFlow 2.19 CPU, about 2 GB; setuptools < 81 is
   needed because SpliceAI imports `pkg_resources`), PharmCAT 3.4.0 (`pharmcat3`; on its first run it downloads its
   GRCh38 reference, 0.9 GB, into its environment). Without the allowlist entry the report still opens, but
   without its styling; the same applies to the PharmCAT report (`pharmcat_pgx`).
3. **Workflows**: *Workflows → Import*, then upload `wes-tertiary-from-vcf.ga` and `wes-tertiary-from-fastq.ga`
   (or paste the raw GitHub URL of each file).

**Galaxy version.** The workflows pin the current wrapper versions (bcftools 1.24, SnpEff 5.4, SnpSift 5.4.0c,
FreeBayes 1.3.10, BWA 0.7.19, MultiQC 1.35). The SnpSift and MultiQC wrappers declare Galaxy profile 25.0, so the
pinned versions need **Galaxy ≥ 25.0**. On older servers, install older revisions of those repositories. Galaxy
then uses the installed version when it opens the workflow and shows a "tool version changed" notice. The
parameters used here are present in older versions too, but check those steps in the editor.

The `.ga` files are prebuilt; nothing needs to be compiled to use them. The `.gxwf.yml` files in [`src/`](src)
are the human-readable sources (Galaxy Format 2). Rebuild only after editing them, which needs Python ≥ 3.10 and
`gxformat2` (a virtual environment avoids the "externally managed environment" error of system pip):
```bash
python3 -m venv .venv && .venv/bin/pip install "gxformat2>=0.29"
```
```bash
.venv/bin/python build.py
```
Import the `.ga` files, not the sources: the v1 source references v2 with `@import`, which only `build.py`
resolves.

**Complete tool states.** The sources set only the parameters that matter. Galaxy would fill every other
parameter with its default when loading a step, and then mark the workflow with "may contain tools which have
changed". To avoid that notice, `build.py` takes each tool step's complete state from
[`src/galaxy-tool-states.json`](src/galaxy-tool-states.json), which holds the states Galaxy computed for these
exact steps. A recorded state is used only if the tool version is the same and the state agrees with every value
the source sets. Otherwise `build.py` keeps the source state and prints a warning. After changing a source or a
tool version, refresh the file:
1. Run `build.py`.
2. Import both `.ga` files into a Galaxy that has the tools installed.
3. Run `python tests/refresh_tool_states.py --galaxy http://192.168.191.12`. It reads the API key from
   `$GALAXY_API_KEY` or `~/.galaxy_admin_key`.
4. Run `build.py` again.

## 2. Reference data (prepare once, keep in a shared history or data library)

Every reference must use **the same genome build and the same chromosome naming** (`chr1` vs `1`) as the reads
and VCFs. A naming mismatch does not raise an error. ClinVar/gnomAD annotations are simply missing.

| Workflow input | What to use | Notes |
|---|---|---|
| **Reference genome (FASTA)** | GRCh38 without ALT contigs, e.g. the *GRCh38 no-alt analysis set* (UCSC naming `chr1…`) | Upload as datatype `fasta` (uncompressed). For VCF input (v2) it must be the genome the variants were called on. |
| **Exome target regions (BED)** | The capture kit's target (or "covered") BED for the same build, from Twist, IDT, Agilent, etc. | No `track`/`browser` header lines. v1 uses it unpadded for coverage QC. The workflows pad it themselves for variant restriction. |
| **Target padding (bp)** | Default `50` | Keeps canonical and extended splice-site variants just outside the targets. |
| **SnpEff database** | Run the **SnpEff download** tool once, e.g. an Ensembl database such as `GRCh38.115` or a RefSeq/MANE one such as `GRCh38.mane.1.2.refseq` | The **SnpEff databases** tool lists the exact names available for SnpEff 5.4. MANE Select databases give clinically conventional NM_ HGVS. SnpEff itself tolerates `chr1` vs `1`. The default Java heap is 8 GB (`GALAXY_MEMORY_MB`). |
| **ClinVar (VCF)** | `clinvar.vcf.gz` from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ | ClinVar uses `1, 2, …, MT`. For a `chr`-named reference, run **bcftools annotate** on it once with *Rename chromosomes* = [`data/ncbi_to_ucsc_chromosomes.tsv`](data/ncbi_to_ucsc_chromosomes.tsv) and output type *compressed BCF*. |
| **Coding regions (BED)** (v1) | Coding exons of one transcript per gene, gene symbol in column 4 (recipe below) | Per-gene coverage in the report: which clinically relevant genes have gaps. |
| **Gene-disease validity (GenCC TSV)** (optional) | [GenCC submissions export](https://search.thegencc.org/download), TSV in the current format (`https://search.thegencc.org/download/action/submissions-export-tsv?format=new`, keyed on `sgc_id`; CC0). The legacy `uuid` format, retired after 30 Sep 2026, is read too | Diseases and mode of inheritance per gene; decides diagnostic vs carrier vs review. |
| **Gene constraint (gnomAD TSV)** (optional) | `gnomad.v4.1.constraint_metrics.tsv` from the gnomAD downloads page | LOEUF/pLI for loss-of-function intolerance (PVS1 applicability). |
| **HPO gene annotations** (optional) | `phenotype_to_genes.txt` from the latest [HPO release](https://github.com/obophenotype/human-phenotype-ontology/releases) | With *Patient phenotype (HPO terms)*, ranks candidates by phenotype match, and gives each term's information content for the disease-level match. |
| **HPO disease annotations (phenotype.hpoa)** (optional) | `phenotype.hpoa` from the same HPO release (35 MB) | Disease names, age of onset, clinical course and inheritance next to each gene; with HPO terms, the disease that best explains the phenotype. |
| **Dosage sensitivity (ClinGen TSV)** (optional) | [`ClinGen_gene_curation_list_GRCh38.tsv`](https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv) (250 kB) | Haploinsufficiency/triplosensitivity scores: decides whether loss of function is a disease mechanism (PVS1, LoF candidates in dominant genes). |
| **ClinVar missense by residue (TSV)** (optional) | Built from ClinVar's `variant_summary.txt.gz` with [`data/build_clinvar_protein.py`](data/build_clinvar_protein.py) (below; 4 MB) | ACMG PS1 (same amino-acid change pathogenic) and PM5 (another pathogenic missense at the residue); rare missense variants with either become candidates. |
| **Mother / Father: variant calls (VCF)** (optional) | Each parent's calls on the same reference, e.g. the *Raw variant calls* (or *Normalised variants*) of the FASTQ workflow run on the parent's reads | Trio analysis (section 4a): de novo variants, parental origin of every finding, compound heterozygotes phased by origin, parentage check. Raw (unfiltered) calls are best: a parental call missing because of a filter would look like a de novo variant. |
| **Mother / Father: aligned reads (BAM)** (optional) | The parents' *Aligned reads (duplicates marked)* | Reads at the proband's variants that a parent's calls lack: a de novo variant is *confirmed* only when both parents have ≥ 15 reads and at most one alt read. |
| **Affected parents** | `unknown` (default), `none`, `mother`, `father` or `both` | A dominant candidate inherited from a parent marked unaffected does not count as fitting the inheritance; one inherited from an affected parent is marked as segregating. |
| **Recurrent artefacts (VCF)** (optional) | The laboratory's own list, built with [`data/build_artefact_list.sh`](data/build_artefact_list.sh) from earlier samples (below) | Alleles called in many of the laboratory's samples but rare or absent in gnomAD: systematic errors of the assay. The report marks matching calls low confidence. |
| **Population allele frequencies (VCF)** | gnomAD v4.1 joint (exome + genome) frequencies reduced to the padded targets, as an indexed **BCF** built with [`data/build_gnomad_extract.sh`](data/build_gnomad_extract.sh) (below) | Transferred fields: `gnomAD_AF`, `_AN`, `_nhomalt` (homozygotes), `_grpmax` / `_grpmax_AF` / `_grpmax_AN` (genetic ancestry group with the highest AF), `_FAF95` / `_FAF95_grp` (highest filtering allele frequency, the statistic of the ACMG BA1/BS1 rules), `_filter` (gnomAD PASS / EXOMES_FILTERED / GENOMES_FILTERED / BOTH_FILTERED) and `_discordance_p` (exome vs genome frequency test). The caller's own `INFO/AF` is untouched. A file without these fields (e.g. an AF-only one) stops the annotation step with "tag not defined". |
| **Missense predictions (AlphaMissense BCF)** | AlphaMissense hg38 predictions as an indexed **BCF** with `INFO/AM_score`, built with [`data/build_alphamissense.sh`](data/build_alphamissense.sh) (below) | Missense PP3/BP4 at the ClinGen-calibrated strengths. CC BY 4.0. |

The gnomAD file is built once with [`data/build_gnomad_extract.sh`](data/build_gnomad_extract.sh) on any Linux
machine with `bcftools` (≥ 1.17) and `bedtools`. It reads only the indexed blocks of the public gnomAD v4.1 joint
VCFs (877 GB in all) that overlap the regions, over plain HTTP from Google Cloud Storage. `regions.bed` uses Ensembl
names and should cover every capture kit you use plus the coding exons (the coding-region BED below), and the
reference is the workflow's FASTA with its `.fai`:
```bash
cat kit1.bed kit2.bed canonical_cds.GRCh38.115.bed | cut -f1-3 > regions.bed
```
```bash
build_gnomad_extract.sh regions.bed GRCh38.fa gnomad.v4.1.joint.targets.GRCh38.bcf 64
```
For two exome kits plus all coding exons (110 Mb with 100 bp padding) this took 1 h 51 min with 64
parallel readers (about 116 MB/s from the bucket) and produced a BCF of 1.6 GB with 65.0 M allele records (the
padding and the gnomAD genomes add many intronic and rare alleles). gnomAD data are CC0. The joint VCF carries no CADD or
SpliceAI scores, so the extract has no non-commercial terms.
Upload the resulting `.bcf` with datatype `bcf`.

The AlphaMissense file is built once with [`data/build_alphamissense.sh`](data/build_alphamissense.sh), which
streams the public 643 MB table (Zenodo record 10813168, CC BY 4.0) into a BCF with Ensembl chromosome names:
```bash
build_alphamissense.sh GRCh38.fa alphamissense.hg38.bcf
```
It took 3.5 min and 2.1 GB of memory and produced 71.7 M alleles in 466 MB. Upload it with datatype `bcf`.
SpliceAI needs no data file: its models and GENCODE annotation come with the tool. Its code (PolyForm Strict) and
models (CC BY-NC 4.0) are for non-commercial use; Illumina's precomputed scores need a BaseSpace account, so the
workflow scores the rare variants itself.

Recipe for the coding-region BED (Ensembl canonical transcripts; use the Ensembl release of the SnpEff database).
For a `chr`-named reference, prefix column 1 with `chr`:
```bash
curl -sL https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz | zcat |
  awk -F'\t' '$3=="CDS" && $9 ~ /tag "Ensembl_canonical"/ && $9 ~ /gene_biotype "protein_coding"/ {
    match($9, /gene_name "[^"]+"/); g = substr($9, RSTART+11, RLENGTH-12); print $1"\t"$4-1"\t"$5"\t"g }' |
  sort -k1,1 -k2,2n -u > canonical_cds.GRCh38.115.bed
```

The `bcftools norm -m -any -f` step is required for any population or ClinVar file: split multiallelic records and
left-align/trim indels against the same reference the reads are aligned to. The workflow normalises the calls, and
annotations are transferred only on an exact `CHROM/POS/REF/ALT` match. If a source keeps multiallelic or padded
records, common indels look absent and pass the rarity filter. In the HG001 test run, before normalising the
population file, homozygous deletions with gnomAD AF 0.94 were reported as rare loss-of-function variants.

The ClinVar residue table is built from the same ClinVar release as the ClinVar VCF (ClinVar is public domain):
```bash
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
python3 data/build_clinvar_protein.py variant_summary.txt.gz clinvar_protein.tsv.gz
```
It keeps the GRCh38 single-nucleotide missense variants classified pathogenic/likely pathogenic or benign/likely
benign (219,540 in the November 2026 release, 49 s, 4.3 MB) with their residue on ClinVar's RefSeq transcript. Upload
it with datatype `tabular`. The HPO and ClinGen files are used as downloaded (HPO: free with attribution; ClinGen:
free to use with attribution).

The recurrent-artefact list is built from the called VCFs of earlier samples sequenced with the same kit, sequencer
and pipeline (the *Raw variant calls* of the FASTQ workflow), with the gnomAD extract above:
```bash
build_artefact_list.sh GRCh38.fa gnomad.v4.1.joint.targets.GRCh38.bcf artefacts.vcf.gz run*/calls.vcf.gz
```
Calls are normalised as in the workflow and must pass its default quality filter. An allele is listed when it is
called in at least 2 samples and in at least 10 % of them (`MIN_SAMPLES`, `MIN_FRACTION`) while its gnomAD AF and
FAF95 are below 0.1 % (`MAX_GNOMAD`) or absent: a real allele carried by one in ten unrelated people cannot be that
rare in gnomAD. Use unrelated individuals only; relatives and a cohort from one founder population share real rare
alleles. Rebuild the list as samples accumulate (three samples took 33 s and 51 MB) and upload it with datatype
`vcf_bgzip`.

## 3. Running

**v1 (FASTQ).** Upload the FASTQ files (`fastqsanger.gz`) and build a *List of Dataset Pairs*. Name each pair
with the sample ID: the element identifier becomes the read-group sample (`SM`) and the VCF sample name. Use
identifiers that end in a letter or digit other than `a f g q s t z` (e.g. `S01`, `NA12878`). The BWA wrapper
strips trailing characters from that set while removing `.gz`/`.fastq` suffixes, so `Patient_tag` would become
`Patient_`. One sample is simply a list with one pair. Samples are processed independently, and the result
collections are keyed by sample.

**Trios.** Run the FASTQ workflow on the parents' reads too (or call them with any caller), then run v2 on the
proband's VCF with the parents' *Raw variant calls* (or *Normalised variants*) and *Aligned reads (duplicates marked)* as the *Mother/Father* inputs and set *Affected
parents*. The FASTQ workflow passes the same inputs to its embedded v2 run.

**v2 (VCF).** Select a VCF/BCF, or select a list collection of VCFs to process a batch. Multi-sample VCFs work.
The tables show the genotype columns of the first sample (`GEN[0]`); change them to `GEN[*]` for all samples.

**Compute notes for v1.**
- BWA-MEM builds the GRCh38 index inside each mapping job: about 1 h and 5–6 GB RAM per sample. If your server
  has a BWA index installed through a data manager, open the *mapping* step and switch *Reference genome* to
  the built-in index to skip this.
- FreeBayes runs in parallel per chromosome using the job's `GALAXY_SLOTS`.

## 4. Outputs

| Output | v1 | v2 | Content |
|---|:-:|:-:|---|
| `candidate_table` | ✔ | ✔ | **Prioritised candidates**, one row per variant: genotype, depth, allelic depths, gene, transcript, consequence, impact, HGVS c./p., LOF/NMD, gnomAD AF, ClinVar significance, review status and disease |
| `candidate_vcf` | ✔ | ✔ | The same candidates as a fully annotated VCF |
| `annotated_table` / `annotated_vcf` | ✔ | ✔ | All on-target, QC-passed variants with the same annotation, for browsing or for custom filters |
| `snpeff_summary`, `variant_stats` | ✔ | ✔ | SnpEff HTML summary and `bcftools stats` |
| `multiqc_report`, `multiqc_stats` | ✔ | | One report covering FastQC, fastp, Picard duplicates and HsMetrics, samtools, bcftools and SnpEff |
| `hs_metrics` | ✔ | | Capture QC: on-target rate, fold-80, `PCT_TARGET_BASES_20X` and others |
| `target_coverage_summary`, `target_coverage_thresholds` | ✔ | | mosdepth summary, and bases covered at 1/10/20/30/50× **per target**, for finding poorly covered exons |
| `deduplicated_bam`, `raw_vcf` | ✔ | | Alignments (duplicates marked) and unfiltered FreeBayes calls |
| `clinical_report` | ✔ | ✔ | **Human-readable HTML report** of the medically significant findings (section 4a) |
| `clinical_findings` | ✔ | ✔ | Every finding of the report as a table: category, gene, variant, genotype, depth, allele fraction, HGVS, gnomAD AF, ClinVar, QC flags, evidence |
| `runs_of_homozygosity` | ✔ | ✔ | bcftools roh regions |
| `coding_coverage_thresholds` | ✔ | | Bases covered at 1/10/20/30/50× per coding exon |

## 4a. The clinical report

The report is produced by explicit, versioned rules from public knowledge bases; there is no statistical or
language model in it, so every statement can be traced to its source (ClinVar record and review status,
GenCC classification, gnomAD frequency, CPIC allele function). It is **decision support for a qualified
reviewer**: every finding that may be acted on must be confirmed by an orthogonal method in an accredited
laboratory.

| Section | Rule |
|---|---|
| Sample and data quality | Ts/Tv, het/hom, het allele-fraction distribution (contamination hint), genetic sex from chrX heterozygosity and chrY coverage, target coverage at 10/20/30× |
| Diagnostic findings | ClinVar pathogenic/likely pathogenic (not conflicting) **and** a genotype that fits the gene's mode of inheritance in GenCC: one allele for dominant genes, hemizygous for X-linked genes in males, two alleles for recessive genes. Genes with both dominant and recessive disorders, X-linked heterozygous females and genes without a GenCC entry go to a *review* table |
| Candidate variants | Not classified in ClinVar (or VUS/conflicting), gnomAD frequency below the cut-off (default 1 %), in a GenCC gene: predicted loss-of-function; missense with AlphaMissense ≥ 0.906 (PP3_Moderate or stronger); missense at a residue with a pathogenic ClinVar missense variant (PS1/PM5); predicted splice-altering outside the canonical sites (SpliceAI ≥ 0.5, or ≥ 0.2 in genes matching the HPO terms); missense/in-frame variants in genes matching the HPO terms. Automated ACMG/AMP evidence codes (PVS1, PS1, PM2_supporting, PM5, BA1, BS1, PP2, PP3, BP4, PP4) for triage, inheritance fit (incl. compound heterozygotes, phased from the reads where possible) and ranking by phenotype match, disease-level phenotype match and evidence |
| Phenotype match by disease | With HPO terms: for every gene with a finding, the HPO-annotated disease that explains the largest share of the patient's terms, each term weighted by its information content; which terms it explains and which it does not, with the disease's onset, course and inheritance |
| Gene knowledge | Next to every gene: GenCC diseases, GenCC disputed/refuted/no-known-relationship assertions, ClinGen haploinsufficiency and triplosensitivity, and (HPO) onset, course and inheritance of its diseases |
| ACMG SF v3.3 | 84 genes (v3.2 plus ABCD1, CYP27A1, PLN; Genet Med 2025;27:101454); biallelic only for MUTYH, ATP7B, BTD, GAA, RPE65, CYP27A1; ABCD1 hemizygous, homozygous or two heterozygous variants (not a single heterozygous variant in a female); HFE p.Cys282Tyr homozygotes only; TTN truncating only. Rare loss-of-function variants not in ClinVar are listed as *review* |
| Carrier status | One P/LP allele in a recessive gene; two P/LP variants shown by the parents or the reads to be in cis count as one allele |
| Trio (with the parents' calls) | Parental origin of every finding (inherited from the mother or father, both parents, de novo, possibly de novo, only one parent for a homozygote); compound heterozygotes phased by origin first, then from the reads; a dominant candidate inherited from a parent marked unaffected does not fit the inheritance; a **de novo table**: rare (gnomAD < 0.01 %) protein-altering or splice-altering variants of good quality absent from both parents' calls, confirmed in both parents' reads (≥ 15 reads, ≤ 1 alt read) or marked *possibly de novo*; three or more parent-absent changes within 10 bp (reads of a paralogous copy) and heterozygous calls without reference reads (repeat or complex loci the parents' calls write differently) stay *possibly de novo* without de novo points; a confirmed de novo variant with an alt fraction well below 0.5 (binomial p < 0.001) is noted as possibly post-zygotic (mosaic); PS2_Moderate (confirmed) / PM6_Supporting (assumed) in phenotype-matched genes after the ClinGen SVI de novo points; a **parentage check** on common SNVs (each parent must carry the proband's homozygous alleles), which withholds PS2/PM6 and puts a warning over the de novo table when it fails |
| Pharmacogenomics | CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, DPYD, TPMT, NUDT15, CYP3A5, CYP2B6, UGT1A1 at fixed GRCh38 positions (checked against ClinVar and the reference). A missing variant counts as the reference allele only with ≥ 15 reads at the position (BAM in v1); otherwise the gene is *indeterminate*. CPIC phenotypes, ClinVar drug-response variants (≥ 1 star). Next to it, the PharmCAT 3.4.0 diplotypes and phenotypes (all defining positions of each star allele, with the positions missing from the VCF listed); PharmCAT's own HTML report carries the CPIC, DPWG and FDA recommendations |
| Risk alleles | Factor V Leiden, prothrombin G20210A, HFE p.Cys282Tyr/p.His63Asp, APOE (flagged: disclose only with counselling and consent), ClinVar risk/association variants (≥ 1 star) |
| Runs of homozygosity | Total length of runs ≥ 1 Mb (consanguinity), runs concentrated on one chromosome (uniparental disomy), homozygous findings inside runs |
| Coverage | Per-gene coding coverage at 10/20×; clinically relevant genes below 95 % at 20× are listed, because a negative result there is not reliable |
| Limits | Lists what an exome pipeline cannot detect (CNVs, SVs, repeat expansions, mtDNA, CYP2D6, HLA, pseudogene regions, mosaicism…) |

**Population frequency.** Every frequency rule (BA1 > 5 %, BS1 above the cut-off, PM2_supporting below 0.01 %
for dominant and 0.1 % for recessive genes, "too common for a penetrant disorder", the candidate cut-off) uses the
higher of the overall gnomAD AF and the highest group filtering AF (FAF95: the lower 95 % confidence bound of the AF
in the most frequent non-bottlenecked group; Whiffin et al. 2017). The FAF stops a variant that is rare overall but
common in one ancestry group from passing as rare, and the overall AF still catches founder variants of the
Ashkenazi or Finnish groups, which FAF excludes. A homozygous candidate in a recessive gene seen homozygous in two or
more gnomAD individuals does not count as fitting the inheritance. Frequencies from sites that failed gnomAD's own
filters, or whose exome and genome frequencies disagree (p < 10⁻⁴), are flagged.

**In-silico predictors.** Following ClinGen SVI, one calibrated tool per mechanism, chosen in advance:

| Evidence | Tool | Rule (strength) |
|---|---|---|
| Missense PP3 | AlphaMissense | ≥ 0.792 Supporting, ≥ 0.906 Moderate, ≥ 0.990 Strong (Bergquist et al., Genet Med 2025; the 0.972–0.989 "+3" interval is reported as Moderate) |
| Missense BP4 | AlphaMissense | ≤ 0.169 Supporting, ≤ 0.099 Moderate, only if SpliceAI ≤ 0.1 or unscored |
| Splicing PP3 / BP4 | SpliceAI 1.3.1, raw, max delta, ±50 bp | ≥ 0.2 PP3, ≤ 0.1 BP4, Supporting only, outside the canonical ±1,2 positions (those go to PVS1; Walker et al., AJHG 2023) |

**Other automated evidence.**

| Evidence | Source | Rule |
|---|---|---|
| PVS1 (applicability) | GenCC, ClinGen dosage, gnomAD LOEUF | Loss of function is a disease mechanism if the gene is recessive or X-linked, has ClinGen haploinsufficiency score 3 or 2, or has LOEUF < 0.6 unless ClinGen scores haploinsufficiency unlikely (40) (Abou Tayoun et al. 2018; Riggs et al. 2020) |
| PS1 | ClinVar residue table | Same amino-acid change, by a different nucleotide change, pathogenic/likely pathogenic with ≥ 1 star, and SpliceAI < 0.2 for this variant |
| PM5 | ClinVar residue table | A different missense change at the same residue pathogenic/likely pathogenic with ≥ 1 star |
| (shown, not scored) | ClinVar residue table | Pathogenic and benign missense variants within ±5 residues; the same change classified benign |

Residues are compared only where the reference amino acid of SnpEff's transcript and ClinVar's RefSeq transcript
agree (98.7 % of HG001's missense variants at a ClinVar residue). A gene-specific PM1 hotspot needs a VCEP
definition, so the ClinVar density is only displayed.

**Read-backed phase.** With the BAM, two heterozygous variants in a recessive gene are phased from the reads or read
pairs that cover both (within 1 kb, MAPQ ≥ 20, duplicates excluded, ≥ 4 informative reads with ≥ 90 % agreeing).
Variants in *cis* are one allele: the patient is a carrier, not affected. Otherwise the phase is reported as unknown
and parental testing decides.

The missense and splicing PP3 are never summed (the stronger one is used), and PP3 is not combined with PVS1.
AlphaMissense's own "likely pathogenic" cut-off (0.564) is inside the calibrated indeterminate range, so it is
not used. Variant-specific VCEP rules (ClinGen CSpec) override these defaults and are not applied automatically.

Each reported variant carries QC flags: low depth, allele imbalance, low QUAL, indel, homopolymer or tandem
repeat context (from the reference), difficult/polymorphic region (HLA, KIR, mucins), segmental duplication / pseudogene genes, recurrent artefact of the laboratory. A call is marked **low confidence** when it has a het alt fraction < 0.25, QUAL < 30, an indel in a homopolymer or tandem repeat of ≥ 8, lies in a paralogous/hypervariable gene, or is on the laboratory's recurrent-artefact list; low-confidence findings stay listed, after the confident ones. In the GIAB validation below every false-positive finding was low confidence and no confirmed finding was.

Intermediate datasets are hidden in the history. Use *show hidden* to inspect them.

## 5. Filters and how to change them

Both filters are workflow parameters, so you can change them on the run form without editing the workflow.

**Quality filter** (bcftools `--include`, applied to calls inside the padded targets):
```
(FILTER="PASS" || FILTER=".") && QUAL>=20 && FMT/DP>=10
```

**Prioritisation filter** (SnpSift):
```
(CLNSIG =~ 'athogenic') | (((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')) & ((! exists gnomAD_AF) | (gnomAD_AF < 0.01)) & ((! exists gnomAD_FAF95) | (gnomAD_FAF95 < 0.01)))
```
This keeps:
- ClinVar *Pathogenic*, *Likely pathogenic* or *Conflicting classifications of pathogenicity*, at any frequency;
- rare (gnomAD AF and highest group FAF95 < 1 %, or absent) variants with HIGH (stop-gain, frameshift, canonical splice, start-loss, …) or
  MODERATE (missense, in-frame indel) impact in any transcript.

Variations:
- Stricter frequency for dominant disorders: replace both `0.01` with `0.001`.
- Only homozygous / hemizygous candidates (recessive model): wrap the expression in parentheses and append the
  genotype test: `(<default expression>) & (GEN[0].GT = '1/1')`.
- Also report LOW-impact splice-region variants: add `| (ANN[*].EFFECT has 'splice_region_variant')` inside the
  impact group.
- **Virtual gene panel**: run **SnpSift Filter** on `candidate_vcf` with *Filter type* = complex,
  expression `ANN[*].GENE in SET[0]`, and a text file of gene symbols (one per line) as the set. Then run
  **SnpSift Extract Fields** with the same field list as the `candidate_table` step.

To add columns, edit the *Fields to extract* of the `candidate_table` / `annotated_table` steps. For example, add
`CLNSIGCONF` and `ALLELEID` (already in the VCF), or `ANN[*].GENE` with one-effect-per-line.

## 6. Design notes

- **Normalise before annotating.** Multiallelic records are split and indels left-aligned against the reference,
  so ClinVar and gnomAD are matched on exact `CHROM/POS/REF/ALT`. Records whose REF does not match the reference
  (a wrong build) are dropped.
- **`gnomAD_*` names.** The population frequencies carry their own prefix, so the caller's `INFO/AF` (genotype
  allele frequency) is not overwritten.
- **Target restriction after calling.** The FreeBayes wrapper starts one FreeBayes process per BED line, which is
  impractical for ~200 000 exome targets. The workflow calls per chromosome and then keeps calls inside the
  padded targets (`bcftools filter -R`, overlapping records included).
- **BWA-MEM instead of BWA-MEM2.** Indexing GRCh38 with BWA-MEM2 needs ~85 GB RAM, so the classic aligner is used
  because it runs on modest servers. Swap in BWA-MEM2 (same read-group options) if you have a prebuilt index.
- **Duplicates are marked, not removed.** FreeBayes ignores reads flagged as duplicates.

## 7. Validation

The workflows were first checked offline:

- **Structure and tool parameters**: [`tests/validate_workflows.py`](tests/validate_workflows.py) checks every
  step of both `.ga` files, including the embedded subworkflow, against the actual wrapper XML from
  [galaxyproject/tools-iuc](https://github.com/galaxyproject/tools-iuc). It uses Galaxy's own parameter models
  (`galaxy-tool-util`) and datatype registry (`galaxy-data`): tool ids and versions, every parameter value and
  conditional case, every connection target and source output, and datatype compatibility (including implicit
  conversions). Result: 0 errors. `--self-test` injects 6 kinds of faults and confirms each is reported.
- **Annotation and filter logic**: [`tests/test_annotation_logic.py`](tests/test_annotation_logic.py) rebuilds the
  bcftools and SnpSift command lines from the parameters stored in `wes-tertiary-from-vcf.ga`. It runs them on
  synthetic data covering:
  - the quality, depth, FILTER, on/off-target, padding and REF-mismatch rules;
  - multiallelic splitting;
  - allele-specific ClinVar/gnomAD matching;
  - 12 prioritisation cases, including the exact 1 % boundary.

  All pass.

```bash
pip install -r tests/requirements.txt
```
```bash
git clone --depth 1 https://github.com/galaxyproject/tools-iuc.git /tmp/tools-iuc
```
```bash
python tests/validate_workflows.py --tools-iuc /tmp/tools-iuc --self-test
```
```bash
python tests/test_annotation_logic.py --snpsift-jar /path/to/SnpSift.jar
```

### End-to-end run in Galaxy (HG001)

v1 was run once in Galaxy 26.2 (aarch64 server) on a GIAB HG001/NA12878 exome. All 23 jobs finished. Inputs:
- **Reads:** GIAB Garvan HiSeq 2500 exome, library NIST7035, lane 1 only (38.2 M reads). The full sample has two
  libraries × two lanes; one lane is a deliberately thin, quick test.
- **Targets:** Illumina Nextera Rapid Capture Expanded Exome BED, lifted from hg19 to GRCh38 with CrossMap
  (201,442 targets, 62.1 Mb).
- **Reference and annotation:** Ensembl GRCh38 primary assembly (`1, 2, …` naming), SnpEff `GRCh38.115`, ClinVar
  (NCBI), and population AF from GATK's `af-only-gnomad.hg38` (gnomAD v2) reduced to the padded targets (since
  replaced by the gnomAD v4.1 extract; the GIAB validation below was rerun with it).

Capture QC (Picard HsMetrics): mean target coverage 23×, median 20×, 77 % of target bases ≥ 10×, 50 % ≥ 20×,
fold-80 2.88, 6.1 % duplicates.

Accuracy against the GIAB v4.2.1 benchmark, inside targets ∩ benchmark regions (56.9 Mb). Calls and truth were
normalised with `bcftools norm` and matched on exact alleles (stricter than hap.py for complex indels):

| Call set | Type | Recall | Precision | F1 | Genotype match |
|---|---|--:|--:|--:|--:|
| `annotated_vcf`, all evaluated targets | SNV | 0.728 | 0.993 | 0.840 | 0.996 |
| `annotated_vcf`, all evaluated targets | indel | 0.504 | 0.828 | 0.627 | 0.925 |
| `annotated_vcf`, targets fully ≥ 10× (13.8 Mb) | SNV | 0.956 | 0.990 | 0.973 | 0.997 |
| `annotated_vcf`, targets fully ≥ 10× (13.8 Mb) | indel | 0.889 | 0.707 | 0.787 | 0.960 |

The low overall recall comes from coverage: with one lane, 23 % of the target bases stay below the `FMT/DP>=10`
quality filter. On well-covered targets SNV recall is 95.6 %. 85 % of the 55,655 reported variants received a
gnomAD AF and 21 % a ClinVar record. The candidate table had 1,357 variants, 104 of them ClinVar (likely)
pathogenic or with conflicting classifications.

Lessons from the run, reflected in section 2:
- Fetch inputs into Galaxy by `https://` URLs. The same NCBI files fetched by `ftp://` URL arrived corrupted.
- SnpEff 5.4 has no `GRCh38.105`; use a name that the **SnpEff databases** tool lists, such as `GRCh38.115`.

### Clinical report on HG001

The same run with the clinical report (26 jobs, all finished). Every non-pharmacogenomic finding in the report
was checked against the GIAB v4.2.1 truth set for HG001:

| Report finding | Report section | QC flags in the report | GIAB v4.2.1 |
|---|---|---|---|
| *PAH* c.1241A>G p.Tyr414Cys het (ClinVar P, expert panel) | Carrier | none | confirmed (0/1) |
| *TPP1* c.509-1G>A het (ClinVar P) | Carrier | low depth (10×) | confirmed (0/1) |
| *PLG* c.112A>G p.Lys38Glu het (ClinVar P/LP) | Review: gene with dominant and recessive disorders | none | confirmed (0/1) |
| *MNS1* c.605delA het (ClinVar LP, 1 star) | Carrier | homopolymer (T×10) | **false positive** (inside benchmark, no variant) |
| *PRSS1* p.Ala16Val het (ClinVar P/LP, gnomAD AF 0.26) | Review: too common for a penetrant disorder | segmental duplication | outside benchmark |
| *PIGQ* c.1680delT hom (not in ClinVar) | Candidate (recessive model) | indel | outside benchmark; GIAB has a nearby SNV only |
| *PRPF3* c.-48-4_-48-3dupTT het | Candidate | low depth, allele balance, homopolymer (T×12) | wrong allele (truth: hom 4-bp insertion) |

All three findings inside the benchmark without a QC flag were confirmed; the false positive and the wrong
allele carry repeat flags. Pharmacogenomics (read depth from the BAM): CYP2C9 \*1/\*2 (intermediate
metabolizer), SLCO1B1 c.521TC (decreased function), CYP2C19 \*2 heterozygous (indeterminate: \*17 position not
covered), CYP4F2, NUDT15, CYP2B6 normal; DPYD, TPMT, VKORC1, CYP3A5 and UGT1A1 indeterminate (a defining position
below 15 reads or off-target). These agree with the CDC GeT-RM consensus for NA12878 (CYP2C9 \*1/\*2 or \*2/\*61,
CYP2C19 \*1/\*2, SLCO1B1 \*1/\*15, which carries c.521T>C). Genetic sex: female; runs of homozygosity 18 Mb.

| Step | Wall time | Peak memory |
|---|--:|--:|
| Clinical report (HG001, all inputs incl. BAM and reference) | 9 s | 1.3 GB |
| bcftools roh | 0.4 s | 14 MB |
| mosdepth over 193,535 coding exons | 41 s | not recorded |

The coverage of this one-lane test (59 % of target bases ≥ 20×) is far below a clinical exome; the report says
so and lists the 79 ACMG SF genes below 95 % at 20×.

### Validation on recent GIAB exomes (NovaSeq 6000, IDT capture)

Recently sequenced exomes of the GIAB reference samples from
[PRJNA734598](https://www.ebi.ac.uk/ena/browser/view/PRJNA734598) (Baid et al., 2021; NovaSeq 6000, IDT
capture, 2×151 bp, about 18.5 M read pairs each), each run through v1 in Galaxy with the IDT target BED of
[DeepVariant's exome case study](https://storage.googleapis.com/deepvariant/exome-case-study-testdata/idt_capture_novogene.grch38.bed).
Truth: NIST GIAB v4.2.1 benchmark (an independent, multi-technology consensus) and, for HG001/NA12878, the CDC
GeT-RM consensus pharmacogenetic genotypes (inter-laboratory testing of Coriell reference DNA).

**Clinical concordance.** The same v2 workflow was also run on each sample's GIAB truth VCF (restricted to the
IDT targets × benchmark regions, quality filter disabled). The findings of the two reports were then matched
variant by variant inside those regions with [`tests/compare_findings.py`](tests/compare_findings.py) (CHROM/POS/REF/ALT). A finding "missed"
is one the truth genotypes support but the sequenced sample's report lacks; an "extra" one is reported from the
sequenced sample but absent from the truth.

**Coverage, accuracy and cost.** Accuracy is allele-exact after normalisation (`bcftools isec`; rtg/hap.py do
not run on this aarch64 server, so complex indels are judged more strictly than with haplotype matching).
Runtime is the sum of the 26 job runtimes; BWA-MEM takes about 61 min per sample, of which roughly 50 min are
spent building the GRCh38 index inside the job (a cached index removes that). Peak memory per job was not
recorded by this Galaxy's job metrics.

| Sample (run) | Mean target cov. / ≥20× | Region | SNV recall / precision | Indel recall / precision | Job runtime | Peak memory |
|---|---|---|--:|--:|--:|--:|
| HG001 (SRR14724473) | 54.6× / 93.0 % | IDT ∩ benchmark, 35.5 Mb | 0.981 / 0.988 | 0.900 / 0.885 | 1.5 h | not recorded |
| | | targets fully ≥ 20×, 32.0 Mb | 0.995 / 0.988 | 0.948 / 0.886 | | |
| HG002 (SRR14724472) | 54.2× / 92.1 % | IDT ∩ benchmark, 36.1 Mb | 0.979 / 0.984 | 0.900 / 0.878 | 1.5 h | not recorded |
| | | targets fully ≥ 20×, 32.5 Mb | 0.994 / 0.984 | 0.946 / 0.880 | | |
| HG005 (SRR14724469) | 53.4× / 91.9 % | IDT ∩ benchmark, 35.6 Mb | 0.979 / 0.989 | 0.918 / 0.916 | 1.5 h | not recorded |
| | | targets fully ≥ 20×, 32.1 Mb | 0.995 / 0.989 | 0.945 / 0.913 | | |

Genotype agreement among matched variants: 99.6–99.7 % for SNVs, 92–94 % for indels. Missed SNVs are mostly
in targets with stretches below 20×.

**Clinical concordance** with report 0.2.0 (ACMG SF v3.3, gnomAD v4.1 joint extract), matched / in truth report;
extra findings in brackets:

| Finding group | HG001 | HG002 | HG005 |
|---|---|---|---|
| ClinVar pathogenic / likely pathogenic | 2/2 (0) | 1/1 (+2) | 4/4 (0) |
| ACMG SF v3.3 (reportable or review) | 0/0 (0) | 0/0 (0) | 1/1 (0) |
| Candidate rare loss-of-function (report 0.2.0; with the predictors of 0.4.0: 6/6, 5/6, 5/5, below) | 1/1 (0) | 1/1 (0) | 3/3 (0) |
| ClinVar drug response | 26/26 (0) | 24/24 (+1) | 20/20 (0) |
| ClinVar risk factor / association | 27/27 (0) | 30/30 (0) | 28/28 (0) |
| Pharmacogenomic genotypes (8 assessable genes + APOE; first run, with the BAM) | all equal | all equal | all equal |

No finding supported by the truth genotypes was missed. All three extra findings, all in HG002 (*MNS1*
c.605delA, alt fraction 0.17 in a T×10 homopolymer; *CNGA1* p.Glu117\*, alt fraction 0.24, QUAL 20; *CYP2D6*
c.941G>A drug-response variant, alt fraction 0.16 in the CYP2D6/CYP2D7 paralog region), are marked **low
confidence** by the report (section 4a), so precision of the high-confidence findings was 100 % in all three
samples. *MNS1* and *CNGA1* now also carry gnomAD's own warnings (site filtered, exome/genome frequencies
discordant). The fourth extra finding of the gnomAD v2 run, *PRSS1* p.Gly177\* in the PRSS1/PRSS2 segmental
duplication, is gone: gnomAD v4.1 has it at AF 0.24 (FAF95 0.39 in African-ancestry samples, filtered in both
exomes and genomes). Findings in the samples' reports that the truth
confirms include *PAH* p.Tyr414Cys and *TPP1* c.509-1G>A (HG001 carrier), *F11* p.Glu135\* (HG002, the
Ashkenazi factor XI deficiency allele), *BBS5* p.Arg89\* (HG005 carrier) and a rare *RYR2* frameshift (HG005,
ACMG SF review), and NUDT15 \*1/\*3 (HG005, intermediate metabolizer for thiopurines).

**gnomAD v2 → v4.1.** The VCF workflow was rerun on each sample's called VCF and truth VCF with the v4.1 extract
and report 0.2.0 (same calls, so the accuracy above is unchanged). The candidate table lost 36–48 % of its rows,
nearly all of them alleles the gnomAD v2 AF-only file had as rare or absent but that are ≥ 1 % overall or in one
genetic ancestry group in v4.1. The only report findings that changed were in the "inheritance not met" group:
*MYO5A*, *IFIH1*, *VPS13B*, *PITRM1*, *HOXD12*, *IFNAR2* and *PRSS1* alleles left (v4.1 AF or FAF95 1.1–24 %),
and an *HLA-DRB1* allele (v4.1 AF 0.14 %, flagged discordant and difficult region) entered.

| Sample | Candidate table rows, called VCF (v2 → v4.1) | Report findings, called VCF (v2 → v4.1) | Clinical concordance (above) | VCF workflow, 13 jobs (called VCF) | gnomAD annotation step | Peak memory, gnomAD annotation |
|---|--:|--:|---|--:|--:|--:|
| HG001 | 1,151 → 722 | 88 → 84 | no finding missed | 5.6 min | 22 s (v2: 8 s) | 18 MB |
| HG002 | 1,298 → 797 | 83 → 81 | no finding missed | 5.5 min | 22 s (v2: 7 s) | not measured |
| HG005 | 1,542 → 797 | 79 → 78 | no finding missed | 5.5 min | 22 s (v2: 8 s) | not measured |

98.9 % of HG001's 54,975 reported variants received a gnomAD v4.1 record. Peak memory was measured outside Galaxy
(`/usr/bin/time`, same input), as this Galaxy does not record job memory.

**Pharmacogenomics against the CDC GeT-RM consensus (HG001 = NA12878).** HG002 and HG005 are not GeT-RM samples.

| Gene | GeT-RM consensus | Report (HG001, SRR14724473) | |
|---|---|---|---|
| CYP2C9 | \*1/\*2 | \*1/\*2, intermediate metabolizer | agrees |
| SLCO1B1 | \*1/\*15 (\*15 carries c.521T>C) | c.521TC, decreased function | agrees |
| CYP3A5 | \*3/\*3 | \*3/\*3, poor metabolizer | agrees |
| CYP4F2 | \*1/\*1 | \*1/\*1 | agrees |
| TPMT | \*1/\*1 | \*1/\*1 | agrees |
| CYP2B6 | \*1/\*1 | c.516GG | agrees |
| DPYD | \*1/(\*4) (\*4 = normal function) | activity score 2, normal | agrees |
| CYP2C19 | \*1/\*2 | \*2 heterozygous; indeterminate (\*17 position off-target) | not assessable |
| VKORC1 | −1639 GA | indeterminate (promoter, off-target) | not assessable |
| UGT1A1 | \*60/\*28 | indeterminate (\*80 promoter proxy off-target) | not assessable |

Seven genes agree, none disagree, three are reported as not assessable rather than guessed.

**Recurrent artefacts.** *PIGQ* c.1680delT (16:582232 CT>C) is called homozygous in all three samples at 11–22×
in a GC-rich target (> 65 % GC). It is in neither GIAB truth set nor gnomAD v4.1 and lies outside the GIAB
benchmark. Without an artefact list it was reported as a high-confidence candidate. The GIAB difficult-region BEDs
do not help: only the broad "all difficult regions" union (625 Mb, 20 % of the genome) covers it, and the narrower
low-mappability/segmental-duplication union also covers confirmed findings (*PLG* p.Lys38Glu, CYP2C19 and UGT2B7
drug-response alleles), so neither can mark calls low confidence.

The recurrent-artefact list was tested leave-one-out: for each sample, the list was built from the other two
samples' calls (`MIN_SAMPLES=2`) and given to the VCF workflow with the sample's own calls.

| Sample | List built from | Alleles on the list | Calls of the sample on the list | … inside IDT ∩ benchmark: false / real | Report findings changed | Clinical concordance | Build time, peak memory (3 samples) | Report step |
|---|---|--:|--:|---|---|---|--:|--:|
| HG001 | HG002, HG005 | 13,925 | 243 | 20 / 6 | *PIGQ* candidate → low confidence | unchanged, no finding missed | 33 s, 51 MB | 22 s |
| HG002 | HG001, HG005 | 14,205 | 243 | 21 / 6 | *PIGQ* candidate → low confidence | unchanged, no finding missed | | 21 s |
| HG005 | HG001, HG002 | 15,263 | 251 | 16 / 6 | *PIGQ* candidate → low confidence | unchanged, no finding missed | | 22 s |

Every listed call outside HLA was a false call (57 of 57 in the benchmark regions, 4–5 % of each sample's false
calls there). The six real alleles on the list are in *HLA-C* and *HLA-B*, whose alleles gnomAD undercounts; the report
already marks the HLA genes low confidence. No other finding changed category or confidence. With only two
samples behind each list this is a lower bound: the list catches more as a laboratory's samples accumulate.

**In-silico predictors and PharmCAT (release 0.4, report 0.4.0).** The VCF workflow was rerun on each sample's
called VCF (with its leave-one-out artefact list) and on its truth VCF, now with AlphaMissense, SpliceAI and PharmCAT.

| Sample | Candidates (inheritance-consistent), 0.3 → 0.4 | of them via AlphaMissense PP3 ≥ Moderate / SpliceAI PP3 | Candidates in IDT ∩ benchmark: matched / in truth report (extra) | Workflow, 16 jobs (called VCF) | SpliceAI step | PharmCAT step | Peak memory |
|---|--:|--:|---|--:|--:|--:|---|
| HG001 | 4 (1) → 15 (4) | 5 / 6 | 6/6 (0) | 20.1 min | 13.9 min (2,464 variants) | 14 s | SpliceAI about 0.7 GB per process × 32; PharmCAT 1.5 GB |
| HG002 | 8 (3) → 15 (5) | 3 / 3 | 5/6 (0) | 20.1 min | 13.8 min | 28 s | as HG001 |
| HG005 | 9 (3) → 15 (3) | 2 / 3 | 5/5 (0) | 16.6 min | 10.4 min | 27 s | as HG001 |

The six runs shared the server, so the SpliceAI times are upper bounds; alone, 32 processes scored HG001's 2,464 rare
variants in 6.2 min. The one candidate missed in HG002 is an 18-bp in-frame duplication in *NBEAL2* that FreeBayes
did not call, not a report rule. All other findings groups were unchanged (ClinVar P/LP, SF, drug response, risk:
as above). New candidates include rare missense variants at PP3_Strong (*DNMT3A* in HG002, AlphaMissense 0.993;
*PKP1* and *PRKRA* in HG005) and splice-region or deep intronic variants with SpliceAI 0.5–0.99 (e.g. *TECPR2*, *GALC*,
*CLCNKB*). Where AlphaMissense is benign but SpliceAI is high (*CR2*, 0.125 / 0.72) the splice PP3 is used. No
reported variant met the BP4 conditions, because rare missense variants with AlphaMissense ≤ 0.169 are not
candidates themselves.

PharmCAT 3.4.0 gave the same diplotype from the sequenced calls as from the GIAB truth VCF for 19 of the 20 genes
both call in the three samples; for HG001 CYP2C19 the exome call adds \*2/\*38 as an alternative to \*1/\*2 because
of positions missing from the exome VCF. From the sequenced exome it also calls HG001's UGT1A1 \*28 (promoter
TA repeat) and CYP4F2, which the truth VCF lacks. Against the GeT-RM consensus for HG001 it agrees for CYP2C19
(\*1/\*2; the built-in genotyper said *indeterminate*), CYP2C9 (\*1/\*2 among PharmCAT's ambiguous set), SLCO1B1 (\*1/\*15
among its set), CYP3A5 \*3/\*3 and DPYD (\*4; c.1627A>G is a normal-function variant), and calls UGT1A1 \*1/\*28 where
GeT-RM has \*60/\*28. For HG005 it calls NUDT15 \*3/\*6, where the built-in table tests only \*3.

**Clinical knowledge (report 0.5.0).** ClinGen dosage sensitivity, the ClinVar residue table (PS1/PM5), HPO
disease annotations and read-backed phase were added. On the three GIAB exomes (rerun of the report on each sample's
release-0.4 annotated VCF, and the full VCF workflow for HG001) every finding kept its category, evidence codes and
confidence: these healthy individuals carry no rare missense variant at a pathogenic ClinVar residue and no pair of
variants in a recessive gene, and no PVS1 changed. The workflow run gave the same findings as the standalone report.
The additions were therefore measured where they act:

*Candidate ranking (spike-in).* One ClinVar pathogenic missense variant (≥ 2 stars, SNV) of a dominant disease gene
was inserted into the HG001 exome as an unclassified variant (no ClinVar fields, absent from gnomAD), with phenotype
terms of an HPO-annotated disease of that gene, 300 cases with one gene each. For PS1/PM5 the variant's own ClinVar
entry was left out. Rank of the inserted variant among all candidates:

| Phenotype given | Report | Rank 1 | Top 5 | Top 10 | Median rank | Time per report (median) | Peak memory |
|---|---|--:|--:|--:|--:|--:|--:|
| 5 disease terms + 1 random term | 0.4.0 | 32/300 | 182/300 | 272/300 | 4.5 | 3.4 s | 0.42 GiB |
| | 0.5.0 without ClinVar residues | 110/300 | 299/300 | 299/300 | 2 | 4.4 s | 0.43 GiB |
| | 0.5.0 | **162/300** | **299/300** | **299/300** | **1** | 5.6 s | 0.55 GiB |
| 3 disease terms + 3 random terms | 0.4.0 | 45/300 | 206/300 | 284/300 | 4 | 3.1 s | |
| | 0.5.0 without ClinVar residues | 56/300 | 215/300 | 286/300 | 3 | 4.1 s | |
| | 0.5.0 | **110/300** | **244/300** | **291/300** | **2** | 5.1 s | |

The inserted variant got PM5 in 110 and 116 of the 300 cases and PS1 in 8 and 19. The reports ran on the rare and
ClinVar-annotated part of the HG001 VCF, 16 in parallel; peak memory was measured on one case (the same for both
phenotype settings), most of the 0.5.0 increase being the ClinVar residue table. Disease-level matching helps most when the phenotype is precise; ClinVar residue evidence helps
regardless of the phenotype.

*PS1/PM5 specificity (ClinVar, leave-one-out).* Each GRCh38 missense SNV with ≥ 2 stars tested against all other
ClinVar entries (≥ 1 star) at its residue:

| ClinVar classification | Variants | PS1 | PM5 | PS1 or PM5 | Time, memory (build of the table) |
|---|--:|--:|--:|--:|---|
| Pathogenic / likely pathogenic | 26,023 | 1,269 (4.9 %) | 14,193 (54.5 %) | 14,650 (56.3 %) | 49 s, < 0.1 GiB |
| Benign / likely benign | 41,609 | 1 (0.0 %) | 126 (0.3 %) | 127 (0.3 %) | |

*Read-backed phase (HG002).* Every pair of heterozygous calls within 1 kb in the HG002 exome
(SRR14724472, realigned with BWA-MEM) whose phase the assembly-based GIAB v5.0q truth gives, 3,000 pairs sampled:

| Distance between the variants | Pairs | Phased (cis or trans) | Wrong phase | Time per pair | Peak memory |
|---|--:|--:|--:|--:|--:|
| ≤ 150 bp | 1,271 | 95.2 % | 0 | 29 ms | 1.09 GiB (the whole check, mostly the VCF) |
| 151–400 bp | 781 | 73.5 % | 0 | | |
| 401–1,000 bp | 948 | 0.9 % | 0 | | |
| SNV pairs / pairs with an indel | 2,629 / 371 | 62.6 % / 39.6 % | 0 / 0 | | |
| **All** | **3,000** | **59.8 % (1,793)** | **0** | | |

All 1,793 phase calls agreed with the truth (1,444 cis, 349 trans); the rest were reported as unknown. Beyond the
read-pair span (about 400 bp for this library) phase stays unknown, so most compound heterozygotes in different
exons still need parental testing.

*Full report on HG001* (all inputs, BAM and reference, median of 3 runs on 4 cores): 0.4.0 7.3 s and 1.09 GiB, 0.5.0
9.2 s and 1.22 GiB.


### Trio analysis (report 0.6.0, GIAB trios)

The parents' IDT exomes from the same NovaSeq 6000 study (HG003 SRR14724471, HG004 SRR14724470, HG006
SRR14724468, HG007 SRR14724467, NA12891 SRR14724466, NA12892 SRR14724464) were run through the FASTQ workflow; the
report was given the proband's release-0.4 annotated VCF, the parents' *Normalised variants* and *Aligned reads
(duplicates marked)*, and *Affected parents* `none`. Parental origin and de novo status were scored against the
GIAB v4.2.1 truth of all three members (HG002's parents with the `noinconsistent` benchmark regions); the
Utah/CEPH parents of HG001 have no GIAB truth.

*Parentage check* (common SNVs, gnomAD AF ≥ 5 %: share of the proband's homozygous sites where the parent carries
the allele, and of its heterozygous sites where at least one parent does; fails below 0.9):

| Proband | Mother given | Father given | Mother | Father | Het sites | Result | Time | Memory |
|---|---|---|--:|--:|--:|---|--:|--:|
| HG002 | HG004 | HG003 | 0.999 | 0.999 | 0.999 | pass | 8 s | in the report |
| HG005 | HG007 | HG006 | 0.999 | 0.999 | 0.998 | pass | 8 s | |
| HG001 | NA12892 | NA12891 | 0.999 | 0.998 | 0.998 | pass | 8 s | |
| HG002 | HG006 (unrelated) | NA12892 (unrelated) | 0.815 | 0.851 | 0.756 | **fail** | 8 s | |
| HG001 | HG006 (unrelated) | NA12892 (the true mother) | 0.805 | 0.999 | 0.835 | **fail** | 8 s | |

16,000–17,500 homozygous and 21,000–25,000 heterozygous sites per proband. With the swapped parents the report
listed 87 "de novo" variants under a parentage warning and applied no PS2/PM6.

*Parental origin* of the proband's true heterozygous autosomal calls inside all three benchmarks (rare = gnomAD
< 1 %, checked in the parents' reads; common variants from the parents' calls only):

| Proband | Variants | Origin assigned | Wrong origin | Truth de novo labelled de novo | Time (rare, with reads) | Memory |
|---|--:|--:|--:|--:|--:|--:|
| HG002, rare | 954 | 954 (100 %) | 1 (0.10 %) | 13 / 13 | 26 s | in the report |
| HG002, all | 15,598 | 15,592 | 29 (0.19 %) | 13 / 13 | | |
| HG005, rare | 1,738 | 1,738 (100 %) | 1 (0.06 %) | 1 / 1 | 44 s | |
| HG005, all | 14,486 | 14,485 | 16 (0.11 %) | 1 / 1 | | |

The wrong origins are mostly "both" for a variant only one parent carries, or a single parent for one both carry
(a parental call present in one family member's call set and not the other's); 7 inherited common variants
missed in a parent's calls (6 in HG002, 1 in HG005) came out as *possibly de novo*. For HG001 the 1,288 rare true heterozygous calls were
assigned maternal 599, paternal 551, both 133 and de novo 5.

*De novo table* (protein- or splice-altering variants; each row checked in the truth of all three members):

| Proband | Truth de novo in the benchmarks | Listed as confirmed de novo | of which true | Listed as possibly de novo | Report time | Memory |
|---|--:|--:|--:|--:|--:|--:|
| HG002 | 13 (11 protein-altering) | 11 | **11** | 2 (no reference reads: TBP CAG repeat, NAALADL2; both outside the benchmarks) | 10.3 s | 1.23 GiB |
| HG005 | 1 (RYR2 frameshift) | 1 | **1** | 6 (two clusters of 3 changes in ANKRD36 / ANKRD36C, outside the benchmarks) | 10.4 s | 1.23 GiB |
| HG001 | no parental truth | 12 | 11 outside HG001's benchmark, 1 inside and absent from its truth | 7 (6 ANKRD36 / ANKRD36C cluster changes) | 11.0 s | 1.24 GiB |
| Without parents | | | | | 5.7–5.9 s | 1.09 GiB |

Every confirmed row in HG002 and HG005 was a true de novo variant, and the true de novo variants not listed were a
synonymous change and a non-coding splice-region change. Before the paralogue-cluster and no-reference rules seven
of these eight artefacts were listed as confirmed de novo. HG002's 13 coding-region de novo variants are more than the
1–2 germline variants expected per exome; DNMT3A p.Cys911Tyr at an alt fraction of 0.27 (marked possibly
mosaic) is the kind acquired in blood or cell culture. HG001's confirmed rows have alt fractions of 0.26–0.38
(4 marked possibly mosaic); without parental truth they cannot be told apart from variants acquired by the cell
line and missing from the GIAB truth. Times are the whole report (GenCC only, 4 cores) including the parents'
read checks.

*End-to-end run in Galaxy:* the VCF workflow (release 0.6) on HG005 with both parents' calls and reads gave the
same de novo result (RYR2 confirmed from 27 and 44 reference reads, the ANKRD36 cluster as possibly de novo) and
passed the parentage check at 99.9 %.

## 8. Scope and limitations

Germline SNVs and short indels in the exome targets only. Not covered: CNV/SV calling, mitochondrial variants
and repeat expansions. The ACMG/AMP evidence codes are an automated subset for triage, not a classification, and
the HPO ranking depends on the terms given. Trio analysis needs the parents' calls on the same reference and
works on the variants the proband's caller found; de novo status at repeats, segmental duplications and other
loci where calls are written differently stays *possibly de novo*, and a parent's low-level mosaicism (a single
alt read among ≥ 15) is read as reference. The reports support expert review; they are not a diagnosis. This is a research workflow,
not a validated clinical pipeline.
