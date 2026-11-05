#!/bin/bash
# a2: the GIAB truth sets and stratifications renamed to the reference's
# Ensembl contig names (1..22, X, Y, MT), the stratification list for hap.py,
# and RTG's SDF of the reference for the vcfeval engine.
set -euo pipefail
V=~/vc; T=$V/truth; mkdir -p $T
E=$V/mamba/envs/happy/bin; export PATH=$E:$PATH
REF=~/certa-bench/idx/ref.fa
map=$T/chr_to_ensembl.txt
{ for c in $(seq 1 22) X Y; do echo "chr$c $c"; done; echo "chrM MT"; } > $map
bed() { sed -e 's/^chrM\t/MT\t/' -e 's/^chr//' "$1"; }
for v in HG002_GRCh38_v5.0q_smvar HG002_GRCh38_CMRG_smallvar_v1.00 HG005_GRCh38_1_22_v4.2.1_benchmark; do
  [ -s $T/$v.vcf.gz.tbi ] || { bcftools annotate --rename-chrs $map -Oz -o $T/$v.vcf.gz $V/dl/truth/$v.vcf.gz && tabix -p vcf $T/$v.vcf.gz; }
done
bed $V/dl/truth/HG002_GRCh38_v5.0q_smvar.benchmark.bed > $T/HG002_v5.0q.bed
bed $V/dl/truth/HG002_GRCh38_CMRG_smallvar_v1.00.bed > $T/HG002_CMRG_v1.00.bed
bed $V/dl/truth/HG005_GRCh38_1_22_v4.2.1_benchmark.bed > $T/HG005_v4.2.1.bed
S=$T/strat; mkdir -p $S
for f in GRCh38_segdups GRCh38_notinsegdups GRCh38_alllowmapandsegdupregions GRCh38_notinalllowmapandsegdupregions \
         GRCh38_alldifficultregions GRCh38_notinalldifficultregions; do
  zcat $V/dl/strat/$f.bed.gz | bed /dev/stdin | grep -v -E '^#|^track|^browser' | gzip > $S/$f.bed.gz
done
awk -v OFS='\t' '{print $1, $2, $3}' $T/HG002_CMRG_v1.00.bed | gzip > $S/CMRG_genes.bed.gz
grep -P '^20\t' $V/targets.bed | gzip > $S/chr20.bed.gz
{ echo -e "segdups\t$S/GRCh38_segdups.bed.gz"; echo -e "notinsegdups\t$S/GRCh38_notinsegdups.bed.gz";
  echo -e "lowmap_or_segdup\t$S/GRCh38_alllowmapandsegdupregions.bed.gz";
  echo -e "not_lowmap_or_segdup\t$S/GRCh38_notinalllowmapandsegdupregions.bed.gz";
  echo -e "alldifficult\t$S/GRCh38_alldifficultregions.bed.gz"; echo -e "notinalldifficult\t$S/GRCh38_notinalldifficultregions.bed.gz";
  echo -e "CMRG_genes\t$S/CMRG_genes.bed.gz"; echo -e "chr20\t$S/chr20.bed.gz"; } > $T/strat.tsv
[ -d $V/ref.sdf ] || rtg format -o $V/ref.sdf $REF > $V/sdf.log 2>&1
for v in HG002_GRCh38_v5.0q_smvar HG002_GRCh38_CMRG_smallvar_v1.00 HG005_GRCh38_1_22_v4.2.1_benchmark; do
  echo "$v: $(bcftools index -s $T/$v.vcf.gz | cut -f1 | tr '\n' ' ' | cut -c1-80)"
done
wc -l $T/*.bed; du -sh $V/ref.sdf; echo "truth ready"
