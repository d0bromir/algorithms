#!/bin/bash
# galaxy: t2t_feat.sh BUILD: for HG002 and HG005, the vc subset of read pairs
# (vc_subset.sh's extraction, from galaxy's copies of the five target CRAMs)
# mapped by BUILD/certa with --iv-rescue local, --t2t and --feature-tags;
# per primary record: qname, mate, contig, pos, MAPQ, XT, XC -> zstd TSV
# (~/vc/t2tfeat_G.tsv.zst; t2t_gate.py reads it). The chain is the
# decompressed hs1ToHg38.over.chain (a build without zlib reads no .gz).
set -euo pipefail
C=$1/certa
V=~/vc; ST=~/micromamba/envs/certa-bench/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; HS1=~/certa-t2t/hs1.cidx; CHAIN=~/vc/hs1ToHg38.over.chain
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
step() { echo "== $(date +%T) $*"; }
for G in HG002 HG005; do
  O=$V/$G; S=$V/subset_$G; mkdir -p $S
  out=$V/t2tfeat_$G.tsv.zst
  [ -s $out ] && continue
  if [ ! -s $S/R2.fq ]; then
    step "$G reads"
    { $ST view -H $O/minibwa.cram
      for n in ivcal pair pairt2t minibwa bwa-mem2; do taskset -c 64-127 $ST view -T $REF -F 0x900 $O/$n.cram; done \
        | awk 'BEGIN{FS=OFS="\t"} { k = $1 "\t" int($2 / 64) % 4; if (!(k in s)) { s[k] = 1; print } }'
    } | taskset -c 64-127 $ST collate -u -O - $S/col | taskset -c 64-127 $ST fastq -n -1 $S/R1.fq -2 $S/R2.fq \
        -s $S/single.fq -0 /dev/null - 2> $S/fastq.log
  fi
  ls -la $S/R1.fq $S/R2.fq $S/single.fq
  step "$G map with --t2t --feature-tags"
  taskset -c 0-63 $TV -o $S/t2tfeat.time $C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all --iv-rescue local \
      --t2t $HS1 --t2t-chain $CHAIN --feature-tags -t 64 -o /dev/stdout 2> $S/t2tfeat.log \
    | taskset -c 64-127 awk 'BEGIN{FS=OFS="\t"} $1 !~ /^@/ && int($2 / 256) % 8 == 0 {
        xt = "-"; xc = "-"
        for (i = 12; i <= NF; i++) { t = substr($i, 1, 5); if (t == "XT:Z:") xt = substr($i, 6); else if (t == "XC:Z:") xc = substr($i, 6) }
        print $1, int($2 / 64) % 4, $3, $4, $5, xt, xc }' \
    | taskset -c 64-127 zstd -q -T8 -o $out.tmp
  mv $out.tmp $out
  cat $S/t2tfeat.time
  rm -f $S/R1.fq $S/R2.fq $S/single.fq
done
step "done"
