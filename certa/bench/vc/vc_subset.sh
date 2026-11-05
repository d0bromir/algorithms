#!/bin/bash
# a2: vc_subset.sh GENOME MODE...: how CERTA's MAPQ affects DeepVariant (a
# diagnostic). The read pairs in GENOME's five target BAMs (their union;
# primary records; both mates present) are mapped again by the frozen PR 19
# build with each MAPQ mode (the placements are the same), and the records
# overlapping the targets go to DeepVariant and hap.py (call.sh) as
# GENOME/S_<mode>.bam. Modes: point, cert, cal, pair. Waits for GENOME's own
# calls (GENOME.call.log) to finish first.
set -euo pipefail
G=$1; shift
V=~/vc; O=$V/$G; S=$V/subset_$G; mkdir -p $S
ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; T=$V/targets.bed; P=~/certa-pm; C=$P/build/certa
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
step() { echo "== $(date +%T) $*"; }
while ! grep -q "^== .* done" $V/$G.call.log 2>/dev/null; do sleep 60; done
step "build $(cd $P && find include src -type f | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -c1-16)"
if [ ! -x $C ]; then
  mkdir -p $P/build/obj
  for f in index seqio mapper_cpu check chain main cmd_map map_config map_pipeline map_stats sam_output gpu_stub; do
    echo "g++ -O3 -DNDEBUG -march=native -std=c++17 -I$P/include -I$P/src -c $P/src/$f.cpp -o $P/build/obj/$f.o"
  done | xargs -P 12 -I{} sh -c {}
  g++ -O3 -march=native $P/build/obj/*.o -lpthread -o $C
fi
step "the reads"
if [ ! -s $S/R2.fq ]; then
  { $ST view -H $O/minibwa.cram
    for n in ivcal pair pairt2t minibwa bwa-mem2; do $ST view -F 0x900 $O/$n.cram; done \
      | awk 'BEGIN{FS=OFS="\t"} { k = $1 "\t" int($2 / 64) % 4; if (!(k in s)) { s[k] = 1; print } }'
  } | $ST collate -u -O - $S/col | $ST fastq -n -1 $S/R1.fq -2 $S/R2.fq -s $S/single.fq -0 /dev/null - 2> $S/fastq.log
fi
wc -l $S/R1.fq $S/single.fq
IV="--intervals --iv-pair all -t 64 -o /dev/stdout"
names=""
for n in "$@"; do
  names="$names S_$n"
  [ -s $O/S_$n.bam.bai ] && continue
  step "map S_$n"
  case $n in
    point) extra="--iv-mapq point" ;;
    cert) extra="--iv-mapq cert" ;;
    cal) extra="--iv-mapq cert --mapq cal" ;;
    pair) extra="--iv-mapq cert --mapq pair" ;;
    *) echo "unknown mode $n"; exit 2 ;;
  esac
  $TV -o $S/S_$n.map.time $C map $IDX $S/R1.fq $S/R2.fq $IV $extra 2> $S/S_$n.log \
    | $ST view -u -L $T - | $ST sort -@ 8 -m 1500M -T $S/sort.$n -o $O/S_$n.bam -
  $ST index $O/S_$n.bam
done
step "DeepVariant"
cd $V && scripts/call.sh $G $names
step "all finished"
