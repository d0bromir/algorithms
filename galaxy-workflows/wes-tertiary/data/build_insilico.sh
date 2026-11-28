#!/usr/bin/env bash
# Join the missense predictions (AlphaMissense + REVEL, data/build_revel.sh) and the CADD extract
# (data/build_cadd.sh) into one file for the workflows' "Missense predictions" input, whose annotation step copies
# every INFO field: INFO/AM_score, INFO/REVEL (missense SNVs) and INFO/CADD (all SNVs and gnomAD indels of the
# coding exons +/- 50 bp).
#
#   build_insilico.sh missense_scores.hg38.bcf cadd.v1.7.cds.bcf insilico_scores.hg38.bcf
#
# Sites of the missense file outside the CADD regions (missense variants of non-canonical transcripts) are kept.
# Licences: AlphaMissense CC BY 4.0; REVEL and CADD free for non-commercial use.
set -euo pipefail
MS=$1 CADD=$2 OUT=$3
TMP=$(mktemp -d "${OUT%.bcf}.tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT
[[ -s $MS.csi || -s $MS.tbi ]] || bcftools index -f "$MS"
[[ -s $CADD.csi ]] || bcftools index -f "$CADD"
# CADD sites with AlphaMissense/REVEL added where they exist
bcftools annotate --no-version -a "$MS" -c INFO/AM_score,INFO/REVEL "$CADD" -Ob -o "$TMP/a.bcf"
# missense sites CADD does not have (exact allele match)
bcftools isec --no-version -C -w1 "$MS" "$CADD" -Ob -o "$TMP/b.bcf"
bcftools index -f "$TMP/a.bcf"; bcftools index -f "$TMP/b.bcf"
bcftools concat --no-version -a "$TMP/a.bcf" "$TMP/b.bcf" -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "$(bcftools index -n "$OUT") alleles -> $OUT ($(bcftools index -n "$TMP/b.bcf") missense-only)"
