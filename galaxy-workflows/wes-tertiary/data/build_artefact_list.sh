#!/usr/bin/env bash
# Build the recurrent-artefact list of the clinical report from the laboratory's own exome calls.
#
#   build_artefact_list.sh reference.fa gnomad.bcf out.vcf.gz calls1.vcf.gz [calls2.vcf.gz ...]
#
# calls*.vcf.gz  germline calls of earlier samples from the same capture kit, sequencer and pipeline (the "Raw variant
#                calls" of the FASTQ workflow, or the caller's VCF given to the VCF workflow); single- or multi-sample.
# gnomad.bcf     the population-frequency extract of the workflows (data/build_gnomad_extract.sh).
#
# An allele is listed when it is called in at least MIN_SAMPLES samples (default 2) and in at least MIN_FRACTION of
# them (default 0.1), while its gnomAD v4.1 frequency (overall AF and grpmax FAF95) is below MAX_GNOMAD (default
# 0.001) or absent. A real allele carried by 10 % of unrelated people cannot be that rare in gnomAD, so such alleles
# are systematic errors of the assay (GC-rich or repetitive targets, alignment or reference problems). The report
# marks matching findings low confidence; it does not remove them. Use samples of unrelated individuals only:
# families and a cohort from one founder population share real rare alleles. Needs bcftools >= 1.17.
#
# Calls are normalised as in the workflows (split, atomised, left-aligned) and must pass the same quality filter
# (QUALITY, a bcftools expression), so the list holds the alleles the report can see.
#
# Output INFO: ARTEFACT_N (samples with the allele), ARTEFACT_HOM (of them homozygous), ARTEFACT_TOTAL (samples
# examined), gnomAD_AF/gnomAD_FAF95 (when gnomAD has the allele).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
REF=$1 GNOMAD=$2 OUT=$3; shift 3
MIN_SAMPLES=${MIN_SAMPLES:-2} MIN_FRACTION=${MIN_FRACTION:-0.1} MAX_GNOMAD=${MAX_GNOMAD:-0.001}
QUALITY=${QUALITY:-'(FILTER="PASS" || FILTER=".") && QUAL>=20 && FMT/DP>=10'}  # the workflows' default quality filter
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# one line per (allele, sample) with a non-reference genotype: CHROM POS REF ALT sample hom
total=0
for vcf in "$@"; do total=$((total + $(bcftools query -l "$vcf" 2> /dev/null | wc -l))); done
[ "$total" -gt 0 ] || { echo "no samples in the input VCFs" >&2; exit 1; }
f=0
for vcf in "$@"; do
    f=$((f + 1))
    bcftools norm --no-version -m -any --atomize -c x -f "$REF" -Ou "$vcf" 2>> "$TMP/norm.log" \
      | bcftools view -i "$QUALITY" -Ou \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
      | awk -v f="$f" 'BEGIN {OFS="\t"} {for (i = 5; i <= NF; i++) if ($i ~ /[1-9]/)
            print $1, $2, $3, $4, f ":" i, ($i ~ /^[1-9][\/|][1-9]$/ && substr($i, 1, 1) == substr($i, 3, 1))}'
done | sort -u -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 > "$TMP/carriers.tsv"

{
    echo "##fileformat=VCFv4.2"
    awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "$REF.fai"
    echo '##INFO=<ID=ARTEFACT_N,Number=1,Type=Integer,Description="Samples of the laboratory set with the allele">'
    echo '##INFO=<ID=ARTEFACT_HOM,Number=1,Type=Integer,Description="Of them, homozygous">'
    echo '##INFO=<ID=ARTEFACT_TOTAL,Number=1,Type=Integer,Description="Samples in the laboratory set">'
    echo "##artefact_list=recurrent in >= $MIN_SAMPLES and >= $MIN_FRACTION of $total samples, gnomAD AF and FAF95 < $MAX_GNOMAD or absent; $(date +%F)"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    awk -v N="$total" -v m="$MIN_SAMPLES" -v fr="$MIN_FRACTION" 'BEGIN {OFS="\t"}
        function out() {if (n >= m && n >= fr * N) print c, p, ".", r, a, ".", ".", "ARTEFACT_N=" n ";ARTEFACT_HOM=" h ";ARTEFACT_TOTAL=" N}
        {k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
         if (k != prev) {if (NR > 1) out(); prev = k; c = $1; p = $2; r = $3; a = $4; n = 0; h = 0}
         n++; h += $6}
        END {if (NR) out()}' "$TMP/carriers.tsv"
} | bcftools sort -Ob -T "$TMP" -o "$TMP/recurrent.bcf" 2> /dev/null
bcftools index "$TMP/recurrent.bcf"
bcftools annotate --no-version -a "$GNOMAD" -c INFO/gnomAD_AF,INFO/gnomAD_FAF95 -Ou "$TMP/recurrent.bcf" \
  | bcftools view --no-version -e "INFO/gnomAD_AF>=$MAX_GNOMAD || INFO/gnomAD_FAF95>=$MAX_GNOMAD" -Oz -o "$OUT"
bcftools index -t -f "$OUT"
echo "$OUT: $(bcftools index -n "$OUT") recurrent alleles from $total samples in $# files"
