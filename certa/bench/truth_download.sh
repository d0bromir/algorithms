#!/bin/bash
# truth_download.sh: phase 1 truth data for VALIDATION_PLAN.md, downloaded to ~/certa-data/truth and
# checked: GIAB/strat files against their published md5s, the S3 assemblies
# against their multipart ETags (8 MiB parts) and gzip -t.
set -u
T=~/certa-data/truth; mkdir -p $T/asm $T/giab $T/strat; cd $T
get() { curl -fsS --retry 8 --retry-delay 10 -C - -o "$2" "$1" || { echo "FAILED $1"; return 1; }; }

S3=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies
for f in hg002v1.1.mat_Y_EBV_MT.fasta.gz hg002v1.1.pat_X_EBV_MT.fasta.gz; do get $S3/$f asm/$f & done

G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest
get $G/checksum.md5 giab/checksum.md5
for f in dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap1.bam dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap1.bam.bai \
         dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap2.bam dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.hap2.bam.bai \
         dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed \
         HG002_GRCh38_v5.0q_smvar.vcf.gz HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi HG002_GRCh38_v5.0q_smvar.benchmark.bed \
         HG002_GRCh38_v5.0q_stvar.vcf.gz HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi HG002_GRCh38_v5.0q_stvar.benchmark.bed \
         NIST_HG002_v5.0q_variant-benchmarksets_README.md; do
  mkdir -p giab/$(dirname $f); get $G/$f giab/$f & done

V=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
get "$V/GRCh38-genome-stratifications-md5s.txt" strat/md5s.txt
for f in SegmentalDuplications/GRCh38_segdups.bed.gz Mappability/GRCh38_lowmappabilityall.bed.gz \
         LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz Union/GRCh38_alldifficultregions.bed.gz \
         Union/GRCh38_alllowmapandsegdupregions.bed.gz Union/GRCh38_notinalldifficultregions.bed.gz \
         OtherDifficult/GRCh38_MHC.bed.gz FunctionalTechnicallyDifficult/GRCh38_CMRGv1.00_falselyduplicatedgenes.bed.gz; do
  mkdir -p strat/$(dirname $f); get "$V/$f" strat/$f & done
wait

echo "== giab md5"; (cd giab && grep -E 'GRCh38_HG2|HG002_GRCh38_v5.0q_(sm|st)var\.(vcf|bench)' checksum.md5 | md5sum -c - 2>&1 | grep -v ': OK$'; echo "giab checked")
echo "== strat md5"; (cd strat && for f in */*.bed.gz; do m=$(grep -F "$(basename $f)" md5s.txt | awk '{print $1}' | head -1); [ "$m" = "$(md5sum < $f | cut -d' ' -f1)" ] && echo "OK $f" || echo "BAD $f (published ${m:-none})"; done)
echo "== assemblies"
etag() { python3 - "$1" <<'EOF'
import hashlib, sys
d = b''; n = 0
with open(sys.argv[1], 'rb') as f:
    while True:
        b = f.read(8 << 20)
        if not b: break
        d += hashlib.md5(b).digest(); n += 1
print(f"{hashlib.md5(d).hexdigest()}-{n}")
EOF
}
for p in "hg002v1.1.mat_Y_EBV_MT.fasta.gz 1ab2fd465678bdc66892b236202065ff-110" "hg002v1.1.pat_X_EBV_MT.fasta.gz 7d212d38df97f42917cce7e959db7b66-110"; do
  set -- $p; e=$(etag asm/$1); [ "$e" = "$2" ] && echo "OK etag $1" || echo "BAD etag $1 ($e)"
  gzip -t asm/$1 && echo "OK gzip $1" || echo "BAD gzip $1"
done
du -sh asm giab strat
echo finished
