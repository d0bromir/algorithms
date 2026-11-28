#!/usr/bin/env bash
# Build the CADD file of the WES workflows: CADD v1.7 PHRED scores for every possible SNV, and for the indels seen in
# gnomAD v4.0 genomes, in the coding exons (plus padding).
#
#   build_cadd.sh regions.bed reference.fa out.bcf [parallel_jobs]
#
# regions.bed  regions to keep, Ensembl chromosome names, e.g. the canonical coding exons (canonical_cds.*.bed);
#              50 bp of padding is added here (splice regions).
# reference.fa the GRCh38 FASTA of the workflows (Ensembl names, .fai next to it).
#
# The CADD files (87 GB of SNVs, 1.3 GB of indels) are not downloaded: tabix reads only the indexed blocks that
# overlap the regions from the Kircher lab server. The output keeps INFO/CADD (PHRED) per allele.
# CADD (Rentzsch et al., Genome Med 2021; Schubach et al., NAR 2024) is free for non-commercial use; commercial use
# needs a licence from its authors (https://cadd.bihealth.org).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
REGIONS=$1 REF=$2 OUT=$3 JOBS=${4:-32}
BASE=${CADD_URL:-https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38}
TMP=$(mktemp -d "${OUT%.bcf}.tmp.XXXX"); TMP=$(cd "$TMP" && pwd)  # absolute: the chunks run inside it
trap 'rm -rf "$TMP"' EXIT

# padded, merged regions (CADD uses Ensembl names: 1..22, X, Y)
cut -f1-3 "$REGIONS" | awk 'NR==FNR {len[$1]=$2; next} ($1 in len) {s=$2-50; if (s<0) s=0; e=$3+50;
    if (e>len[$1]) e=len[$1]; print $1"\t"s"\t"e}' "$REF.fai" - \
  | sort -k1,1V -k2,2n | bedtools merge -i - > "$TMP/regions.bed"

# Remote reads are latency-bound: regions less than 1 kb apart are fetched together and the fetch regions are cut
# into chunks of $CHUNK consecutive regions that are read in parallel; rows outside the padded regions are dropped.
bedtools merge -d "${GAP:-1000}" -i "$TMP/regions.bed" \
  | awk -v n="${CHUNK:-300}" -v d="$TMP" '{if ($1 != c || k == n) {i++; k = 0; c = $1}
    k++; f = sprintf("%s/part%05d.bed", d, i); print >> f; close(f)}'

one_chunk() {  # partNNNNN.bed -> partNNNNN.tsv (SNVs and indels; retried: remote reads can drop)
    local bed=$1 n f ok
    set -o pipefail
    for n in 1 2 3 4 5 6 7 8; do
        ok=1
        : > "${bed%.bed}.tsv"
        for f in whole_genome_SNVs.tsv.gz gnomad.genomes.r4.0.indel.tsv.gz; do
            # chrom pos ref alt raw phred -> BED-like (chrom, pos-1, pos, ref, alt, phred) for the region filter
            if ! tabix -R "$bed" "$BASE/$f" 2> "${bed%.bed}.log" \
                | awk 'BEGIN {OFS="\t"} !/^#/ {print $1, $2 - 1, $2, $3, $4, $6}' \
                | bedtools intersect -u -a - -b "$TMP/regions.bed" 2>> "${bed%.bed}.log" \
                | cut -f1,3-6 >> "${bed%.bed}.tsv"; then ok=0; break; fi
        done
        if ((ok)); then echo "$(date +%T) $(basename "${bed%.bed}"): $(wc -l < "${bed%.bed}.tsv") alleles"; return 0; fi
        echo "$(date +%T) $(basename "$bed"): attempt $n failed ($(tail -1 "${bed%.bed}.log")), retrying" >&2
        sleep $((15 * n))
    done
    return 1
}
export -f one_chunk; export TMP BASE REF
cd "$TMP"  # tabix writes the downloaded .tbi index into the working directory
curl -fsSL -o whole_genome_SNVs.tsv.gz.tbi "$BASE/whole_genome_SNVs.tsv.gz.tbi"
curl -fsSL -o gnomad.genomes.r4.0.indel.tsv.gz.tbi "$BASE/gnomad.genomes.r4.0.indel.tsv.gz.tbi"
ls part*.bed | xargs -P "$JOBS" -I{} bash -c 'one_chunk {}'
cd - >/dev/null

# TSV (chrom pos ref alt phred) -> sorted VCF -> BCF with INFO/CADD
{
    printf '##fileformat=VCFv4.2\n'
    printf '##source=CADD GRCh38-v1.7 (free for non-commercial use; https://cadd.bihealth.org)\n'
    printf '##INFO=<ID=CADD,Number=A,Type=Float,Description="CADD v1.7 PHRED score">\n'
    awk '{print "##contig=<ID="$1",length="$2">"}' "$REF.fai"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    cat "$TMP"/part*.tsv | awk 'BEGIN {OFS="\t"} {print $1, $2, ".", $3, $4, ".", ".", "CADD=" $5}' \
        | sort -k1,1V -k2,2n -k4,4 -k5,5 -u -S 2G -T "$TMP"
} | bcftools view --no-version -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "$(bcftools index -n "$OUT") alleles -> $OUT"
