#!/bin/bash
# a2: vc_call.sh GENOME NAME...: DeepVariant 1.10.0 (WGS model) on
# ~/vc/GENOME/NAME.cram (or .bam) over the target regions, the commands run_deepvariant
# runs (64 make_examples shards, call_variants, postprocess_variants), then
# hap.py 0.3.15 (vcfeval engine) against the GIAB truth: HG002 against v5.0q
# and CMRG v1.00, HG005 against v4.2.1; stratified (segdups, low mappability,
# all difficult, CMRG genes, chr20), PASS calls.
set -euo pipefail
G=$1; shift
V=~/vc; O=$V/$G; T=$V/truth
REF=~/certa-bench/idx/ref.fa; TG=$V/targets.bed; M=$V/dv/rootfs/opt/models/wgs; D=$V/scripts/dvrun.sh
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
happy() {  # happy TRUTH BED OUT [extra]
  local truth=$1 bed=$2 out=$3; shift 3
  PATH=$V/mamba/envs/happy/bin:$PATH HGREF=$REF $V/mamba/envs/happy/bin/python $V/mamba/envs/happy/bin/hap.py \
    $truth $O/$n.vcf.gz -f $bed -r $REF -o $out --engine vcfeval --engine-vcfeval-template $V/ref.sdf \
    --stratification $T/strat.tsv --threads 32 "$@" > $out.log 2>&1
}
for n in "$@"; do
  echo "== $(date +%T) $G $n"
  if [ ! -s $O/$n.vcf.gz ]; then
    X=$O/$n.dv; rm -rf $X; mkdir -p $X
    IN=$O/$n.cram; [ -s $IN ] || IN=$O/$n.bam
    seq 0 63 | $TV -o $O/$n.dv_me.time xargs -P 64 -I{} $D make_examples --mode calling --ref $REF --reads $IN \
      --examples $X/me.tfrecord@64.gz --checkpoint $M --nocall_small_model_examples --regions $TG --sample_name $G \
      --task {} > $O/$n.dv_me.log 2>&1
    OMP_NUM_THREADS=64 $TV -o $O/$n.dv_cv.time $D call_variants --outfile $X/cv.tfrecord.gz --examples $X/me.tfrecord@64.gz \
      --checkpoint $M > $O/$n.dv_cv.log 2>&1
    $TV -o $O/$n.dv_pp.time $D postprocess_variants --ref $REF --infile $X/cv.tfrecord.gz --outfile $O/$n.vcf.gz.tmp.vcf.gz \
      --cpus 64 --regions $TG --sample_name $G > $O/$n.dv_pp.log 2>&1
    mv $O/$n.vcf.gz.tmp.vcf.gz $O/$n.vcf.gz
    [ -e $O/$n.vcf.gz.tmp.vcf.gz.tbi ] && mv $O/$n.vcf.gz.tmp.vcf.gz.tbi $O/$n.vcf.gz.tbi
    rm -rf $X
  fi
  if [ $G = HG002 ]; then
    happy $T/HG002_GRCh38_v5.0q_smvar.vcf.gz $T/HG002_v5.0q.bed $O/$n.v5q -T $TG
    happy $T/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz $T/HG002_CMRG_v1.00.bed $O/$n.cmrg
  else
    happy $T/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $T/HG005_v4.2.1.bed $O/$n.v421 -T $TG
  fi
done
echo "== $(date +%T) done"
