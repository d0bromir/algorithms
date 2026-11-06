#!/bin/bash
# a2 (x86-64): tools and truth for the variant-calling test of CERTA's MAPQ
# options. DeepVariant 1.10.0 and hap.py 0.3.15 from bioconda via micromamba
# (no containers on a2); GIAB truth sets and stratifications. Downloads go to
# their own directories under ~/vc/dl, scripts live in ~/vc/scripts.
set -euo pipefail
V=~/vc
mkdir -p $V/bin $V/dl/truth $V/dl/strat
if [ ! -x $V/bin/micromamba ]; then
  curl -fsSL -o $V/bin/micromamba https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
  chmod +x $V/bin/micromamba
fi
$V/bin/micromamba --version
export MAMBA_ROOT_PREFIX=$V/mamba
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
cd $V/dl/truth
for f in AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz \
         AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi \
         AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed \
         AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md \
         AshkenazimTrio/HG002_NA24385_son/v5.0q/checksum.md5 \
         AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz \
         AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi \
         AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed \
         ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
         ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
         ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed; do
  [ -s $(basename $f) ] || curl -fsS -O $G/$f
done
grep smvar checksum.md5 | grep GRCh38 | md5sum -c - || true
S=$G/genome-stratifications/v3.6/GRCh38@all
cd $V/dl/strat
for f in SegmentalDuplications/GRCh38_segdups.bed.gz SegmentalDuplications/GRCh38_notinsegdups.bed.gz \
         Union/GRCh38_alllowmapandsegdupregions.bed.gz Union/GRCh38_notinalllowmapandsegdupregions.bed.gz \
         Union/GRCh38_alldifficultregions.bed.gz Union/GRCh38_notinalldifficultregions.bed.gz \
         GRCh38-genome-stratifications-md5s.txt; do
  [ -s $(basename $f) ] || curl -fsS -O $S/$f
done
echo "== $(date +%T) envs"
[ -d $MAMBA_ROOT_PREFIX/envs/dv ] || $V/bin/micromamba create -y -n dv -c conda-forge -c bioconda \
  deepvariant=1.10.0 samtools bcftools htslib > $V/dv_env.log 2>&1
[ -d $MAMBA_ROOT_PREFIX/envs/happy ] || $V/bin/micromamba create -y -n happy -c conda-forge -c bioconda \
  hap.py=0.3.15 rtg-tools > $V/happy_env.log 2>&1
$V/bin/micromamba list -n dv | grep -E "deepvariant|tensorflow|samtools" || true
$V/bin/micromamba list -n happy | grep -E "hap.py|rtg" || true
du -sh $V/mamba $V/dl
echo "== $(date +%T) setup finished"
