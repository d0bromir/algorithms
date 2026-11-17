#!/usr/bin/env bash
# Build the missense-predictor file of the WES workflows from the AlphaMissense hg38 predictions.
#
#   build_alphamissense.sh reference.fa out.bcf
#
# reference.fa  the GRCh38 FASTA of the workflows (Ensembl names 1..22, X, Y, MT; .fai next to it).
#
# AlphaMissense (Cheng et al., Science 2023) scores every possible single-nucleotide missense change of the canonical
# transcripts (71 M); the predictions are CC BY 4.0 (https://zenodo.org/records/10813168). The 643 MB table is
# streamed, not stored. The output keeps one INFO field per allele, AM_score (0..1, higher = more likely pathogenic),
# which the clinical report maps to the ClinGen-calibrated PP3/BP4 strengths of Bergquist et al., Genet Med 2025.
# Needs bcftools >= 1.17, curl and awk.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
REF=$1 OUT=$2
URL=${AM_URL:-https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz}
TMP=$(mktemp -d "${OUT%.bcf}.tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT

{
    echo "##fileformat=VCFv4.2"
    awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "$REF.fai"
    echo '##INFO=<ID=AM_score,Number=A,Type=Float,Description="AlphaMissense pathogenicity (0-1), canonical transcript; CC BY 4.0, Cheng et al. 2023">'
    echo "##alphamissense=$(basename "$URL"), $(date +%F)"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    # input columns: CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
    curl -sfL --retry 5 "$URL" | gzip -dc \
      | awk 'BEGIN {OFS="\t"} /^#/ {next} {c = $1; sub(/^chr/, "", c); if (c == "M") c = "MT";
             print c, $2, ".", $3, $4, ".", ".", "AM_score=" $9}'
} | bcftools sort -m 2G -T "$TMP" -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "$OUT: $(bcftools index -n "$OUT") alleles, $(du -h "$OUT" | cut -f1)"
