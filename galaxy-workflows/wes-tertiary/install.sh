#!/usr/bin/env bash
# Install everything the WES workflows need on a Galaxy server: command-line programs, Tool Shed tools, the local
# tools and their dependencies, the workflows, and the reference and clinical knowledge files.
#
#   ./install.sh [options] STAGE [STAGE ...]
#
# Stages (run in the order given; "all" = system venv tools local-tools workflows data upload):
#   check        report what is installed and what is missing; changes nothing
#   system       command-line programs for the data scripts (bcftools, bedtools, samtools, tabix, curl, gawk,
#                python3-venv) from apt, or into a conda environment with --conda
#   venv         Python virtual environment with Ephemeris and BioBlend ($VENV)
#   tools        Tool Shed repositories from tools.yml (Ephemeris shed-tools)
#   local-tools  copy tools/* to the local tool directory, check the local tool config and the HTML allowlist
#                (changed only with --configure), then let Galaxy install the tools' Conda dependencies
#   workflows    import the two .ga workflows (with --update: new version of workflows already present)
#   data         download and build the reference and knowledge files into --data-dir, with MANIFEST.tsv
#   upload       upload the --data-dir files into a new Galaxy history and fetch the SnpEff database there
#
# Options:
#   --galaxy URL            Galaxy server (default http://localhost)
#   --data-dir DIR          reference data directory (default /data/galaxy/wes-ref/release-<YYYY-MM>)
#   --local-tools-dir DIR   tool_path of the local tool config (default /srv/galaxy/local_tools)
#   --galaxy-config DIR     Galaxy config directory (default /etc/galaxy)
#   --galaxy-user USER      Unix user Galaxy runs as (default galaxy)
#   --regions BED           regions for the gnomAD extract (capture kit BEDs, Ensembl names); without it the
#                           gnomAD extract (about 2 h, 1.6 GB) is skipped
#   --jobs N                parallel gnomAD readers (default 32)
#   --ensembl N             Ensembl release of the FASTA, the coding-exon BED and SnpEff (default 115)
#   --history NAME          history name for "upload" (default "WES reference data <data-dir name>")
#   --previous-history NAME "upload": copy the files that did not change from this history (the previous release)
#                           instead of uploading them again (Galaxy shares their storage)
#   --conda                 "system": install the programs with conda/mamba into $PREFIX instead of apt
#   --update                "workflows": update workflows that already exist
#   --configure             "local-tools": write a missing local_tool_conf.xml and add the report tools to the
#                           HTML allowlist; without it the Galaxy configuration is only checked, never changed
#   --dry-run               print the commands that change something instead of running them
#
# The Galaxy admin API key is read from $GALAXY_API_KEY or ~/.galaxy_admin_key, and handed to the Python helper in
# its environment, never on a command line. Every stage can be re-run: what is already in place is left alone.
# Needs bash, sudo (for apt and the Galaxy directories) and an admin account on the Galaxy server.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GALAXY_URL="http://localhost"
DATA_DIR="/data/galaxy/wes-ref/release-$(date +%Y-%m)"
LOCAL_TOOLS="/srv/galaxy/local_tools"
GX_CONF="/etc/galaxy"
GX_USER="galaxy"
REGIONS=""
JOBS=32
ENSEMBL=115
HISTORY=""
PREV_HISTORY=""
USE_CONDA=0
UPDATE=0
CONFIGURE=0
DRY=0
VENV="${VENV:-$HOME/wes-venv}"
PREFIX="${PREFIX:-$HOME/wes-env}"
LOCAL_TOOL_IDS=(wes_clinical_report spliceai_score pharmcat_pgx wes_cnv wes_mito wes_couple wes_hla wes_yhaplo wes_kir)

die() { echo "install.sh: $*" >&2; exit 1; }
say() { printf '\n== %s\n' "$*"; }
run() { if ((DRY)); then printf '[dry-run] %s\n' "$*"; else "$@"; fi; }
have() { command -v "$1" >/dev/null 2>&1; }
as_root() { if [[ $EUID -eq 0 ]]; then run "$@"; else run sudo "$@"; fi; }

STAGES=()
while (($#)); do
    case "$1" in
        --galaxy) GALAXY_URL="$2"; shift 2 ;;
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        --local-tools-dir) LOCAL_TOOLS="$2"; shift 2 ;;
        --galaxy-config) GX_CONF="$2"; shift 2 ;;
        --galaxy-user) GX_USER="$2"; shift 2 ;;
        --regions) REGIONS="$(realpath "$2")"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --ensembl) ENSEMBL="$2"; shift 2 ;;
        --history) HISTORY="$2"; shift 2 ;;
        --previous-history) PREV_HISTORY="$2"; shift 2 ;;
        --conda) USE_CONDA=1; shift ;;
        --update) UPDATE=1; shift ;;
        --configure) CONFIGURE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0 ;;
        all) STAGES+=(system venv tools local-tools workflows data upload); shift ;;
        check|system|venv|tools|local-tools|workflows|data|upload) STAGES+=("$1"); shift ;;
        *) die "unknown argument '$1' (see --help)" ;;
    esac
done
((${#STAGES[@]})) || die "no stage given (see --help)"
[[ -d "$PREFIX/bin" ]] && export PATH="$PREFIX/bin:$PATH"

api_key() {
    if [[ -n "${GALAXY_API_KEY:-}" ]]; then return; fi
    [[ -r "$HOME/.galaxy_admin_key" ]] || die "no API key: set GALAXY_API_KEY or create ~/.galaxy_admin_key"
    GALAXY_API_KEY="$(tr -d ' \n\r\t' < "$HOME/.galaxy_admin_key")"
    export GALAXY_API_KEY
}
helper() {  # Python helper with the API key in its environment
    [[ -x "$VENV/bin/python" ]] || die "no virtual environment at $VENV: run the 'venv' stage first"
    api_key
    run env GALAXY_URL="$GALAXY_URL" "$VENV/bin/python" "$HERE/install/galaxy_setup.py" "$@"
}

PROGRAMS=(bcftools bedtools samtools tabix bgzip curl gawk python3)
stage_check() {
    say "programs"
    local p
    for p in "${PROGRAMS[@]}"; do
        if have "$p"; then printf '  %-9s %s\n' "$p" "$(command -v "$p")"; else printf '  %-9s MISSING\n' "$p"; fi
    done
    if have bcftools; then
        local v; v="$(bcftools --version | head -1 | awk '{print $2}')"
        printf '  bcftools version %s (data scripts need >= 1.17)\n' "$v"
    fi
    say "Python environment $VENV"
    if [[ -x "$VENV/bin/python" ]]; then "$VENV/bin/pip" list 2>/dev/null | grep -iE '^(ephemeris|bioblend) ' | sed 's/^/  /'
    else echo "  missing"; fi
    say "local tools in $LOCAL_TOOLS"
    local t
    for t in "${LOCAL_TOOL_IDS[@]}"; do
        if [[ -f "$LOCAL_TOOLS/$t/$t.xml" ]]; then
            if diff --strip-trailing-cr -rq "$HERE/tools/$t" "$LOCAL_TOOLS/$t" >/dev/null 2>&1; then echo "  $t: installed, same as the repository"
            else echo "  $t: installed, DIFFERS from the repository"; fi
        else echo "  $t: not installed"; fi
    done
    say "data directory $DATA_DIR"
    if [[ -f "$DATA_DIR/MANIFEST.tsv" ]]; then column -t -s $'\t' "$DATA_DIR/MANIFEST.tsv" | sed 's/^/  /'
    else echo "  no MANIFEST.tsv (run the 'data' stage)"; fi
    if [[ -x "$VENV/bin/python" ]] && { [[ -n "${GALAXY_API_KEY:-}" ]] || [[ -r "$HOME/.galaxy_admin_key" ]]; }; then
        say "Galaxy $GALAXY_URL"
        helper check "$HERE/tools.yml" "${LOCAL_TOOL_IDS[@]}"
    fi
}

stage_system() {
    say "command-line programs"
    if ((USE_CONDA)); then
        local c; c="$(command -v micromamba || command -v mamba || command -v conda || true)"
        [[ -n "$c" ]] || die "--conda: no micromamba, mamba or conda on PATH"
        if [[ -d "$PREFIX/conda-meta" ]]; then echo "  $PREFIX exists"
        else run "$c" create -y -p "$PREFIX" -c conda-forge -c bioconda \
                 "bcftools>=1.17" bedtools samtools htslib curl gawk python=3.12; fi
        return
    fi
    have apt-get || die "no apt-get: use --conda, or install ${PROGRAMS[*]} with the system's package manager"
    local missing=() p
    for p in bcftools bedtools samtools tabix curl gawk python3-venv ca-certificates; do
        dpkg -s "$p" >/dev/null 2>&1 || missing+=("$p")
    done
    if ((${#missing[@]} == 0)); then echo "  all present"; return; fi
    echo "  installing: ${missing[*]}"
    as_root apt-get update -q
    as_root env DEBIAN_FRONTEND=noninteractive apt-get install -y -q "${missing[@]}"
}

stage_venv() {
    say "Python environment $VENV"
    [[ -x "$VENV/bin/python" ]] || run python3 -m venv "$VENV" \
        || die "python3 -m venv failed: install python3-venv (stage 'system')"
    run "$VENV/bin/pip" install -q --upgrade pip
    run "$VENV/bin/pip" install -q "ephemeris>=0.10" "bioblend>=1.2"
}

stage_tools() {
    say "Tool Shed repositories (tools.yml)"
    helper tools "$HERE/tools.yml"
}

stage_local_tools() {
    say "local tools -> $LOCAL_TOOLS"
    local t owner restart=0
    owner="$(stat -c %U "$LOCAL_TOOLS" 2>/dev/null || echo "$GX_USER")"
    [[ -d "$LOCAL_TOOLS" ]] || as_root install -d -o "$GX_USER" -g "$GX_USER" "$LOCAL_TOOLS"
    for t in "${LOCAL_TOOL_IDS[@]}"; do
        if diff --strip-trailing-cr -rq "$HERE/tools/$t" "$LOCAL_TOOLS/$t" >/dev/null 2>&1; then echo "  $t: up to date"; continue; fi
        echo "  $t: copying"
        as_root mkdir -p "$LOCAL_TOOLS/$t"
        as_root cp -r "$HERE/tools/$t/." "$LOCAL_TOOLS/$t/"
        as_root chown -R "$owner:" "$LOCAL_TOOLS/$t"
        restart=1
    done

    say "Galaxy configuration in $GX_CONF"
    local conf="$GX_CONF/local_tool_conf.xml" allow="$GX_CONF/sanitize_allowlist.txt" yml="$GX_CONF/galaxy.yml"
    local todo=0
    if as_root test -f "$conf"; then
        for t in "${LOCAL_TOOL_IDS[@]}"; do
            if as_root grep -q "$t/$t.xml" "$conf"; then echo "  $conf lists $t"
            else echo "  ADD to $conf:  <tool file=\"$t/$t.xml\"/>"; todo=1; fi
        done
    elif ((! CONFIGURE)); then
        echo "  MISSING $conf (re-run with --configure to write it, listing the local tools)"; todo=1
    else
        echo "  writing $conf"
        local tmp; tmp="$(mktemp)"
        {
            echo '<?xml version="1.0"?>'
            echo "<toolbox tool_path=\"$LOCAL_TOOLS\">"
            echo '    <section id="wes_analysis_local" name="WES analysis">'
            for t in "${LOCAL_TOOL_IDS[@]}"; do echo "        <tool file=\"$t/$t.xml\"/>"; done
            echo '    </section>'
            echo '</toolbox>'
        } > "$tmp"
        as_root install -m 0644 "$tmp" "$conf"; rm -f "$tmp"
        todo=1
    fi
    if as_root test -f "$yml"; then
        if as_root grep -q "local_tool_conf.xml" "$yml"; then echo "  galaxy.yml loads local_tool_conf.xml"
        else echo "  ADD $conf to tool_config_file in $yml"; todo=1; fi
        if ! as_root grep -q "sanitize_allowlist_file" "$yml"; then
            echo "  ADD to $yml:  sanitize_allowlist_file: $allow"; todo=1; fi
    fi
    for t in wes_clinical_report pharmcat_pgx; do  # optional: without it the report opens without its styling
        if as_root test -f "$allow" && as_root grep -qx "$t" "$allow"; then echo "  allowlist has $t"
        elif ((CONFIGURE)); then echo "  adding $t to $allow"; as_root sh -c "echo $t >> '$allow'"; restart=1
        else echo "  optional: add the line $t to $allow (or Admin -> Manage Allowlist) for the report's styling"; fi
    done
    if ((todo)); then
        echo "  Galaxy does not load the local tools yet: make the changes above (or re-run with --configure),"
        echo "  restart Galaxy, then re-run 'local-tools'. Galaxy configured by Ansible: put them in its templates."
        return
    fi
    if ((restart)); then echo "  restart Galaxy to load the changes above (copied tools, allowlist)"; fi
    say "dependencies of the local tools (Conda, through Galaxy; SpliceAI's TensorFlow is about 2 GB)"
    helper local-deps "${LOCAL_TOOL_IDS[@]}"
}

stage_workflows() {
    say "workflows"
    local upd=(); ((UPDATE)) && upd=(--update)
    helper workflows "${upd[@]}" "$HERE/wes-tertiary-from-vcf.ga" "$HERE/wes-tertiary-from-fastq.ga"
}

stage_data() {
    say "reference and knowledge files -> $DATA_DIR"
    local p
    for p in bcftools samtools curl gawk python3; do have "$p" || die "$p missing: run the 'system' stage"; done
    ((DRY)) && { echo "[dry-run] would download into $DATA_DIR"; return; }
    mkdir -p "$DATA_DIR"; cd "$DATA_DIR"
    local M=MANIFEST.tsv
    [[ -f $M ]] || printf 'file\tgalaxy_datatype\tsource\tdate\tbytes\n' > $M
    manifest() {  # file datatype source; an unchanged file (same size, same source) keeps its date
        local f="$1" bytes d
        bytes="$(stat -c %s "$f")"
        d="$(awk -F'\t' -v f="$f" -v b="$bytes" -v s="$3" '$1 == f && $5 == b && $3 == s {print $4}' $M)"
        awk -F'\t' -v f="$f" '$1 != f' $M > $M.tmp && mv $M.tmp $M
        printf '%s\t%s\t%s\t%s\t%s\n' "$f" "$2" "$3" "${d:-$(date +%F)}" "$bytes" >> $M
    }
    fetch() {  # url out: resumable download, kept only when complete
        [[ -s "$2" ]] && { echo "  $2: present"; return; }
        echo "  $2 <- $1"
        curl -fL --no-progress-meter --retry 8 --retry-all-errors -C - -o "$2.part" "$1"
        mv "$2.part" "$2"
    }
    local ENS="https://ftp.ensembl.org/pub/release-$ENSEMBL"
    local FA="GRCh38.ensembl$ENSEMBL.primary_assembly.fa"
    if [[ ! -s $FA ]]; then
        fetch "$ENS/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" "$FA.gz"
        gzip -dc "$FA.gz" | sed 's/^>\([^ ]*\) .*/>\1/' > "$FA.tmp" && mv "$FA.tmp" "$FA" && rm -f "$FA.gz"
    fi
    [[ -s $FA.fai ]] || samtools faidx "$FA"
    manifest "$FA" fasta "$ENS/fasta/homo_sapiens/dna/ (Ensembl names, headers trimmed)"

    local CDS="canonical_cds.GRCh38.$ENSEMBL.bed"
    if [[ ! -s $CDS ]]; then
        echo "  $CDS <- Ensembl $ENSEMBL GTF (canonical protein-coding CDS)"
        curl -fsL --retry 8 "$ENS/gtf/homo_sapiens/Homo_sapiens.GRCh38.$ENSEMBL.gtf.gz" | gzip -dc |
          gawk -F'\t' '$3=="CDS" && $9 ~ /tag "Ensembl_canonical"/ && $9 ~ /gene_biotype "protein_coding"/ {
            match($9, /gene_name "[^"]+"/); g = substr($9, RSTART+11, RLENGTH-12); print $1"\t"$4-1"\t"$5"\t"g }' |
          sort -k1,1 -k2,2n -u > "$CDS.tmp" && mv "$CDS.tmp" "$CDS"
    fi
    manifest "$CDS" bed "$ENS/gtf/homo_sapiens/"

    local CV=https://ftp.ncbi.nlm.nih.gov/pub/clinvar
    fetch "$CV/vcf_GRCh38/clinvar.vcf.gz" clinvar.GRCh38.vcf.gz
    manifest clinvar.GRCh38.vcf.gz vcf_bgzip "$CV/vcf_GRCh38/clinvar.vcf.gz ($(gzip -dc clinvar.GRCh38.vcf.gz | grep -m1 '^##fileDate' | cut -d= -f2))"
    if [[ ! -s clinvar_protein.tsv.gz ]]; then
        fetch "$CV/tab_delimited/variant_summary.txt.gz" variant_summary.txt.gz
        python3 "$HERE/data/build_clinvar_protein.py" variant_summary.txt.gz clinvar_protein.tsv.gz
        rm -f variant_summary.txt.gz
    fi
    manifest clinvar_protein.tsv.gz tabular "$CV/tab_delimited/variant_summary.txt.gz via data/build_clinvar_protein.py"

    fetch "https://search.thegencc.org/download/action/submissions-export-tsv?format=new" gencc-submissions.tsv
    manifest gencc-submissions.tsv tabular "https://search.thegencc.org/download (CC0)"
    local GN=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1
    fetch "$GN/constraint/gnomad.v4.1.constraint_metrics.tsv" gnomad.v4.1.constraint_metrics.tsv
    manifest gnomad.v4.1.constraint_metrics.tsv tabular "$GN/constraint/ (CC0)"
    local HPO=https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download
    fetch "$HPO/phenotype_to_genes.txt" phenotype_to_genes.txt
    local hpo_v; hpo_v="$(curl -fsI "$HPO/phenotype.hpoa" | grep -io '/download/v[0-9-]*' | head -1 | cut -d/ -f3 || true)"
    manifest phenotype_to_genes.txt tabular "$HPO/ (${hpo_v:-release unknown})"
    fetch "$HPO/phenotype.hpoa" phenotype.hpoa
    manifest phenotype.hpoa tabular "$HPO/ (${hpo_v:-release unknown})"
    fetch https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv ClinGen_gene_curation_list_GRCh38.tsv
    manifest ClinGen_gene_curation_list_GRCh38.tsv tabular "https://ftp.clinicalgenome.org/"
    fetch https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv ClinGen_region_curation_list_GRCh38.tsv
    manifest ClinGen_region_curation_list_GRCh38.tsv tabular "https://ftp.clinicalgenome.org/ (recurrent CNV and population regions)"
    local MANE_V=1.5 MN=https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human
    fetch "$MN/release_$MANE_V/MANE.GRCh38.v$MANE_V.summary.txt.gz" "MANE.GRCh38.v$MANE_V.summary.txt.gz"
    manifest "MANE.GRCh38.v$MANE_V.summary.txt.gz" tabular "$MN/release_$MANE_V/ (public domain)"
    if [[ ! -s gene_disease.tsv ]]; then
        local G2P=https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads g2p_dir
        g2p_dir="$(curl -fsL "$G2P/" | grep -o 'href="20[0-9_]*/"' | cut -d'"' -f2 | sort | tail -1)"
        fetch https://search.clinicalgenome.org/kb/gene-validity/download clingen_gene_validity.csv
        fetch "$G2P/$g2p_dir" g2p.listing.html
        local g2p_file; g2p_file="$(grep -o 'href="allG2P_[^"]*\.csv\.gz"' g2p.listing.html | cut -d'"' -f2 | head -1)"
        fetch "$G2P/$g2p_dir$g2p_file" "$g2p_file"
        fetch https://www.orphadata.com/data/xml/en_product6.xml en_product6.xml
        fetch https://www.orphadata.com/data/xml/en_product9_ages.xml en_product9_ages.xml
        python3 "$HERE/data/build_gene_disease.py" clingen_gene_validity.csv "$g2p_file" en_product6.xml \
            en_product9_ages.xml gene_disease.tsv
        rm -f g2p.listing.html clingen_gene_validity.csv allG2P_*.csv.gz en_product6.xml en_product9_ages.xml
    fi
    if [[ ! -s vcep_rules.tsv ]]; then
        fetch "https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed" erepo.tsv
        python3 "$HERE/data/build_vcep.py" --fetch cspec erepo.tsv vcep_rules.tsv erepo_variants.tsv
        rm -rf cspec erepo.tsv
    fi
    manifest vcep_rules.tsv tabular "ClinGen CSpec registry via data/build_vcep.py (ClinGen terms of use)"
    manifest erepo_variants.tsv tabular "ClinGen Evidence Repository via data/build_vcep.py (ClinGen terms of use)"
    fetch https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz         Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
    manifest Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz tabular "https://zenodo.org/records/6347673 (CC BY 4.0)"
    if [[ ! -s mave_evidence.tsv.gz ]]; then
        local MV; MV="$(curl -fsL https://zenodo.org/api/records/20840937 | python3 -c 'import json,sys; print([f["links"]["self"] for f in json.load(sys.stdin)["files"] if f["key"].endswith(".tar.gz")][0])')"
        fetch "$MV" mavedb.tar.gz
        python3 "$HERE/data/build_mavedb.py" mavedb.tar.gz mave_evidence.tsv && gzip -f mave_evidence.tsv
        rm -f mavedb.tar.gz
    fi
    manifest mave_evidence.tsv.gz tabular "MaveDB bulk download (CC0, Zenodo 20840937) via data/build_mavedb.py"
    if [[ ! -s gnomad.v3.1.chrM.vcf.gz ]]; then  # chrM -> the reference's mitochondrial name (Ensembl: MT)
        fetch https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz gnomad.chrM.vcf.bgz
        local mt; mt="$(cut -f1 "$FA.fai" | grep -m1 -xE 'MT|chrM')"
        gzip -dc gnomad.chrM.vcf.bgz | sed -e "s/^chrM\t/$mt\t/" -e "s/^##contig=<ID=chrM,/##contig=<ID=$mt,/" | bgzip > gnomad.v3.1.chrM.vcf.gz
        rm -f gnomad.chrM.vcf.bgz
    fi
    manifest gnomad.v3.1.chrM.vcf.gz vcf_bgzip "gnomAD v3.1 chrM sites (CC0), renamed to the reference's mitochondrial contig"
    fetch https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz gnomad.v4.1.cnv.all.vcf.gz
    manifest gnomad.v4.1.cnv.all.vcf.gz vcf_bgzip "gnomAD v4.1 exome CNVs (CC0)"
    if [[ ! -s uniprot_features.tsv.gz ]]; then
        fetch "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=tsv&fields=accession,gene_primary,length,xref_ensembl,ft_domain,ft_region,ft_motif,ft_zn_fing,ft_dna_bind,ft_act_site,ft_binding,ft_repeat,ft_coiled" uniprot_human.tsv
        python3 "$HERE/data/build_uniprot_features.py" uniprot_human.tsv uniprot_features.tsv && gzip -f uniprot_features.tsv
        rm -f uniprot_human.tsv
    fi
    manifest uniprot_features.tsv.gz tabular "UniProtKB/Swiss-Prot human (CC BY 4.0) via data/build_uniprot_features.py"
    fetch https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat hla.dat
    local hla_v; hla_v="$(grep -m1 -o 'Release Version [0-9.]*' hla.dat | cut -d' ' -f3 || true)"
    manifest hla.dat txt "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/ (IPD-IMGT/HLA ${hla_v:-release unknown}, CC BY-ND 4.0; HLA typing with wes_hla)"
    fetch https://ftp.ebi.ac.uk/pub/databases/ipd/kir/kir.dat kir.dat
    local kir_v; kir_v="$(grep -m1 -o 'IPD-KIR Release Version [0-9.]*' kir.dat | cut -d' ' -f4 || true)"
    manifest kir.dat txt "https://ftp.ebi.ac.uk/pub/databases/ipd/kir/ (IPD-KIR ${kir_v:-release unknown}, CC BY-ND 4.0; KIR typing with wes_kir)"
    manifest gene_disease.tsv tabular "ClinGen gene-disease validity (CC0), G2P (CC0), Orphanet (CC BY 4.0) via data/build_gene_disease.py"
    local GS=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
    if [[ ! -s GRCh38_segdups.bed ]]; then  # Ensembl names, like the reference
        fetch "$GS/SegmentalDuplications/GRCh38_segdups.bed.gz" GRCh38_segdups.bed.gz
        gzip -dc GRCh38_segdups.bed.gz | sed "s/^chr//" > GRCh38_segdups.bed && rm -f GRCh38_segdups.bed.gz
    fi
    manifest GRCh38_segdups.bed bed "$GS/SegmentalDuplications/ (GIAB v3.6, public domain; CNV confidence)"

    if [[ ! -s alphamissense.hg38.bcf ]]; then
        echo "  alphamissense.hg38.bcf <- Zenodo 10813168 (about 4 min)"
        bash "$HERE/data/build_alphamissense.sh" "$FA" alphamissense.hg38.bcf
    fi
    manifest alphamissense.hg38.bcf bcf "https://zenodo.org/records/10813168 (CC BY 4.0) via data/build_alphamissense.sh"
    if [[ ! -s missense_scores.hg38.bcf ]]; then
        echo "  missense_scores.hg38.bcf <- AlphaMissense + REVEL v1.3 (Zenodo 7072866; about 10 min)"
        bash "$HERE/data/build_revel.sh" "$FA" alphamissense.hg38.bcf missense_scores.hg38.bcf
    fi
    manifest missense_scores.hg38.bcf bcf "AlphaMissense (CC BY 4.0) + REVEL v1.3 (non-commercial) via data/build_revel.sh; use as the missense predictions"

    if [[ -n "$REGIONS" ]]; then
        if [[ ! -s gnomad.v4.1.joint.targets.GRCh38.bcf ]]; then
            have bedtools || die "bedtools missing: run the 'system' stage"
            cat "$REGIONS" "$CDS" | cut -f1-3 > regions.bed
            echo "  gnomad.v4.1.joint.targets.GRCh38.bcf <- gnomAD v4.1 joint VCFs (about 2 h with 64 readers)"
            bash "$HERE/data/build_gnomad_extract.sh" regions.bed "$FA" gnomad.v4.1.joint.targets.GRCh38.bcf "$JOBS"
        fi
        manifest gnomad.v4.1.joint.targets.GRCh38.bcf bcf "$GN/ joint sites VCFs via data/build_gnomad_extract.sh (CC0)"
        local GB=https://storage.googleapis.com/gcp-public-data--gnomad
        [[ -s regions.bed ]] || cat "$REGIONS" "$CDS" | cut -f1-3 > regions.bed
        if [[ ! -s gnomad_exome_coverage.tsv ]]; then
            echo "  gnomad_exome_coverage.tsv <- gnomAD v4.0 exome coverage summary (4 GB streamed; about 6 min)"
            curl -sfL "$GB/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz" | gzip -dc |
              python3 "$HERE/data/build_gnomad_sites.py" coverage regions.bed gnomad_exome_coverage.tsv
        fi
        manifest gnomad_exome_coverage.tsv tabular "gnomAD v4.0 exome coverage (CC0) via data/build_gnomad_sites.py"
        if [[ ! -s gnomad_pext.tsv ]]; then
            echo "  gnomad_pext.tsv <- gnomAD GTEx v10 base-level pext (198 MB streamed; about 7 min)"
            curl -sfL "$GB/release/4.1/pext/gnomad.pext.gtex_v10.base_level.tsv.gz" | gzip -dc |
              python3 "$HERE/data/build_gnomad_sites.py" pext regions.bed gnomad_pext.tsv
        fi
        manifest gnomad_pext.tsv tabular "gnomAD GTEx v10 pext (CC0) via data/build_gnomad_sites.py"
        if [[ ! -s carrier_genes.tsv ]]; then
            echo "  carrier_genes.tsv <- ClinVar P/LP x gnomAD v4.1 extract, curated genes, ClinGen pediatric actionability"
            fetch "https://actionability.clinicalgenome.org/ac/Pediatric/api/summ?flavor=flat" actionability_pediatric.json
            [[ -s gnomad.v4.1.joint.targets.GRCh38.bcf.csi ]] || bcftools index -f gnomad.v4.1.joint.targets.GRCh38.bcf
            bcftools view -i 'INFO/CLNSIG~"athogenic" && INFO/CLNSIG!~"onflicting"' clinvar.GRCh38.vcf.gz -Oz -o clinvar_plp.vcf.gz
            tabix -f -p vcf clinvar_plp.vcf.gz
            bcftools annotate -a gnomad.v4.1.joint.targets.GRCh38.bcf                 -c INFO/gnomAD_AF,INFO/gnomAD_grpmax_AF,INFO/gnomAD_grpmax,INFO/gnomAD_filter clinvar_plp.vcf.gz |
              python3 "$HERE/data/build_carrier_genes.py" gene_disease.tsv actionability_pediatric.json carrier_genes.tsv
            rm -f clinvar_plp.vcf.gz clinvar_plp.vcf.gz.tbi
        fi
        manifest carrier_genes.tsv tabular "ClinVar, gnomAD v4.1, ClinGen/G2P/Orphanet, ClinGen actionability via data/build_carrier_genes.py"
    else
        echo "  gnomAD extract skipped (give --regions kit.bed; see README section 2)"
    fi
    say "done: $DATA_DIR ($(du -sh . | cut -f1))"
    column -t -s $'\t' $M | sed 's/^/  /'
    cd - >/dev/null
}

stage_upload() {
    say "upload $DATA_DIR into Galaxy"
    [[ -f "$DATA_DIR/MANIFEST.tsv" ]] || die "no $DATA_DIR/MANIFEST.tsv: run the 'data' stage first"
    helper upload "$DATA_DIR" "${HISTORY:-WES reference data $(basename "$DATA_DIR")}" "GRCh38.$ENSEMBL" "$PREV_HISTORY"
}

for s in "${STAGES[@]}"; do "stage_${s//-/_}"; done
