#!/usr/bin/env bash
# Add REVEL to the missense-predictor file of the WES workflows.
#
#   build_revel.sh reference.fa alphamissense.hg38.bcf missense_scores.hg38.bcf
#
# REVEL v1.3 (Ioannidis et al., AJHG 2016; https://zenodo.org/records/7072866, free for non-commercial use) scores every
# possible single-nucleotide missense change (82 M rows, GRCh38 positions since 2020). The 667 MB archive is streamed;
# one score per allele is kept (the highest where transcripts disagree). The output joins it with the AlphaMissense
# BCF of build_alphamissense.sh: INFO/AM_score and INFO/REVEL, which the clinical report maps to the ClinGen-calibrated
# PP3/BP4 strengths (Pejaver et al., AJHG 2022 for REVEL; Bergquist et al., Genet Med 2025 for AlphaMissense).
# Needs bcftools >= 1.17, curl, unzip and awk.
set -euo pipefail
REF=$1 AM=$2 OUT=$3
URL=${REVEL_URL:-https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip}
TMP=$(mktemp -d "${OUT%.bcf}.tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT

curl -sfL --retry 5 -o "$TMP/revel.zip" "$URL"
{
    echo "##fileformat=VCFv4.2"
    awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "$REF.fai"
    echo '##INFO=<ID=REVEL,Number=A,Type=Float,Description="REVEL v1.3 missense pathogenicity (0-1); Ioannidis et al. 2016, non-commercial use">'
    echo "##revel=$(basename "$URL"), $(date +%F)"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    # columns: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
    unzip -p "$TMP/revel.zip" | awk -F, 'BEGIN {OFS="\t"} NR > 1 && $3 != "." && $3 != "" {
        k = $1 SUBSEP $3 SUBSEP $4 SUBSEP $5
        if (k == last) { if ($8 > best) best = $8; next }
        if (last != "") { split(last, a, SUBSEP); print a[1], a[2], ".", a[3], a[4], ".", ".", "REVEL=" best }
        last = k; best = $8 }
      END { if (last != "") { split(last, a, SUBSEP); print a[1], a[2], ".", a[3], a[4], ".", ".", "REVEL=" best } }'
} | bcftools sort -m 2G -T "$TMP" -Ob -o "$TMP/revel.bcf"
bcftools index -f "$TMP/revel.bcf"
bcftools merge -m none "$AM" "$TMP/revel.bcf" -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "$OUT: $(bcftools index -n "$OUT") alleles, $(du -h "$OUT" | cut -f1)"
