#!/bin/bash
# a2: vote_sweep.sh: HG002's lost read pairs mapped by ~/certa-fix with seed
# votes off and at several settings; the reads still lost (lost_diag.py) and
# the CPU time.
set -uo pipefail
cd ~/vc/diag/lost
ST=~/vc/dv/rootfs/opt/conda/envs/bio/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; C=~/certa-fix/build/certa
# lost_run.sh's lost read pairs (names.txt) from the subset's reads
[ -s lost_R2.fq ] || for r in 1 2; do
  awk 'NR == FNR { w["@" $1] = 1; next } FNR % 4 == 1 { k = ($1 in w) } k' names.txt ~/vc/subset_HG002/R$r.fq > lost_R$r.fq
done
for s in "0 8" "64 8" "256 8" "1024 8" "256 16" "4096 8"; do
  set -- $s
  n=v$1_$2
  /usr/bin/time -f "%e %U %S %M" -o $n.time $C map $IDX lost_R1.fq lost_R2.fq --intervals --iv-pair all --iv-rescue local \
    --feature-tags --witness --iv-vote-occ $1 --iv-vote-top $2 -t 64 -o $n.sam 2> $n.log
  lost=$(python3 ~/vc/scripts/lost_diag.py $ST $REF mb.sam $n.sam 2>&1 > $n.md | head -1)
  echo "occ $1 top $2: $(awk '{printf "%.1f s CPU, %.2f GiB", $2 + $3, $4 / 1048576}' $n.time) | $lost | $(grep -c '^| search' $n.md) search rows"
done
