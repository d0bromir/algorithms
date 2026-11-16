#!/usr/bin/env bash
# Build the population-frequency file of the WES workflows from the gnomAD v4.1 joint (exomes + genomes) sites VCFs.
#
#   build_gnomad_extract.sh regions.bed reference.fa out.bcf [parallel_jobs]
#
# regions.bed  regions to keep, Ensembl chromosome names (1..22, X, Y), e.g. the capture kit BED(s) plus the
#              canonical coding exons; 100 bp of padding is added here.
# reference.fa the GRCh38 FASTA used by the workflows (Ensembl names, .fai next to it).
#
# The gnomAD VCFs (877 GB in total) are not downloaded: bcftools reads only the indexed blocks that overlap the
# regions, straight from the public bucket (two kits plus all coding exons, 110 Mb padded,
# took 1 h 51 min with 64 parallel readers). Needs bcftools >= 1.17,
# bedtools and awk. The output (1.6 GB, 65 M alleles for those regions) keeps, per allele:
#   gnomAD_AF, gnomAD_AC, gnomAD_AN, gnomAD_nhomalt   joint exome+genome allele frequency, counts, homozygotes
#   gnomAD_grpmax, gnomAD_grpmax_AF, gnomAD_grpmax_AN genetic ancestry group with the highest AF (non-bottlenecked)
#   gnomAD_FAF95, gnomAD_FAF95_grp                    highest filtering allele frequency (95% CI lower bound,
#                                                     the statistic of the ACMG BA1/BS1 frequency rules)
#   gnomAD_filter                                     gnomAD FILTER: PASS, EXOMES_FILTERED, GENOMES_FILTERED, BOTH_FILTERED
#   gnomAD_discordance_p                              exome-vs-genome frequency discordance p-value (< 1e-4: unreliable)
# gnomAD is released under CC0; the in-silico scores of the per-dataset VCFs (CADD, SpliceAI) are not in the joint
# VCF, so the extract carries no non-commercial terms.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1  # bcftools links OpenBLAS; one thread per process
REGIONS=$1 REF=$2 OUT=$3 JOBS=${4:-64}
# Plain http on purpose: over https libcurl negotiates HTTP/2 with Google storage, and long range reads then fail
# with "Error in the HTTP2 framing layer". The files are public and the BGZF/CRC checks catch corruption.
URL=${GNOMAD_URL:-http://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr%s.vcf.bgz}
TMP=$(mktemp -d "${OUT%.bcf}.tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT

# padded, merged regions with gnomAD (chr) names
cut -f1-3 "$REGIONS" | awk 'NR==FNR {len[$1]=$2; next} ($1 in len) {s=$2-100; if (s<0) s=0; e=$3+100;
    if (e>len[$1]) e=len[$1]; print "chr"$1"\t"s"\t"e}' "$REF.fai" - \
  | sort -k1,1V -k2,2n | bedtools merge -i - > "$TMP/regions.chr.bed"
for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > "$TMP/chrs.txt"

cat > "$TMP/rename.txt" <<'EOF'
INFO/AF_joint gnomAD_AF
INFO/AC_joint gnomAD_AC
INFO/AN_joint gnomAD_AN
INFO/nhomalt_joint gnomAD_nhomalt
INFO/grpmax_joint gnomAD_grpmax
INFO/AF_grpmax_joint gnomAD_grpmax_AF
INFO/AN_grpmax_joint gnomAD_grpmax_AN
INFO/fafmax_faf95_max_joint gnomAD_FAF95
INFO/fafmax_faf95_max_gen_anc_joint gnomAD_FAF95_grp
INFO/stat_union_p_value gnomAD_discordance_p
EOF
KEEP=INFO/AF_joint,INFO/AC_joint,INFO/AN_joint,INFO/nhomalt_joint,INFO/grpmax_joint,INFO/AF_grpmax_joint
KEEP=$KEEP,INFO/AN_grpmax_joint,INFO/fafmax_faf95_max_joint,INFO/fafmax_faf95_max_gen_anc_joint,INFO/stat_union_p_value

# Remote reads are latency-bound (about one second per region), so regions less than 1 kb apart are fetched
# together (-R, then -T keeps only the padded regions themselves), and the fetch regions are cut into chunks of
# $CHUNK consecutive regions on one chromosome that are read in parallel.
bedtools merge -d "${GAP:-1000}" -i "$TMP/regions.chr.bed" \
  | awk -v n="${CHUNK:-200}" -v d="$TMP" '{if ($1 != c || k == n) {i++; k = 0; c = $1}
    k++; f = sprintf("%s/part%05d.%s.bed", d, i, substr($1, 4)); print >> f; close(f)}'

one_chunk() {  # partNNNNN.CHR.bed -> partNNNNN.CHR.bcf (retried: remote reads can drop)
    local bed=$1 c n
    set -o pipefail  # a dropped connection must fail the chunk, not truncate it
    c=${bed%.bed}; c=${c##*.}
    for n in 1 2 3 4 5 6 7 8; do
        if bcftools view --no-version -R "$bed" --regions-overlap pos -T "$TMP/regions.chr.bed" -i 'AC_joint>0' -Ou \
                "$(printf "$URL" "$c")" 2> "${bed%.bed}.view.log" \
            | bcftools annotate --no-version -x "^$KEEP" --rename-annots "$TMP/rename.txt" --rename-chrs "$TMP/chrs.txt" -Ov \
            | awk 'BEGIN {OFS="\t"}
                   /^##INFO=<ID=gnomAD_AF,/ {print "##INFO=<ID=gnomAD_filter,Number=1,Type=String,Description=\"gnomAD v4.1 joint FILTER (PASS, EXOMES_FILTERED, GENOMES_FILTERED, BOTH_FILTERED)\">"}
                   /^#/ {print; next} {$8 = $8 ";gnomAD_filter=" $7; $7 = "."; print}' \
            | bcftools norm --no-version -m -any -c x -f "$REF" -Ob -o "${bed%.bed}.bcf" 2> "${bed%.bed}.log"; then
            echo "$(date +%T) $(basename "${bed%.bed}"): $(bcftools view -H "${bed%.bed}.bcf" | wc -l) records"
            return 0
        fi
        echo "$(date +%T) $(basename "$bed"): attempt $n failed ($(tail -1 "${bed%.bed}.view.log")), retrying" >&2
        rm -f "${bed%.bed}.bcf"; sleep $((15 * n))
    done
    return 1
}
export -f one_chunk; export TMP URL KEEP REF
ls "$TMP"/part*.bed | xargs -P "$JOBS" -I{} bash -c 'one_chunk {}'

ls "$TMP"/part*.bcf | sort > "$TMP/parts.txt"
[ "$(wc -l < "$TMP/parts.txt")" -eq "$(ls "$TMP"/part*.bed | wc -l)" ] || { echo "missing chunks" >&2; exit 1; }
bcftools concat --no-version -f "$TMP/parts.txt" -Ou \
  | bcftools annotate --no-version -h <(echo "##gnomAD_extract=gnomAD v4.1 joint sites, $(date +%F), regions $(basename "$REGIONS") +100 bp") \
      -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "$(date +%T) $OUT: $(bcftools index -n "$OUT") records, $(du -h "$OUT" | cut -f1)"
