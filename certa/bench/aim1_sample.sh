#!/bin/bash
# Aim 1 measurement for one sample: q22 (GPU) baseline and best config, then
# the q15 index on each baseline's uncertified reads (CPU, 16 threads).
#   aim1_sample.sh NAME READS.fq
set -u
cd ~/certa-bench
C=./certa-sl/certa/build/certa
N=$1; R=$2; O=aim1/$N; mkdir -p $O
for base in "def:--local 2" "best:--budget3 4096 --local 2"; do
  b=${base%%:*}; ba=${base#*:}
  $C map idx/grch38.s1v3.cidx $R -t 16 --gpu --device 1 $ba -o /dev/null -u $O/unc_$b.fq 2> $O/q22_$b.log
  nu=$(( $(wc -l < $O/unc_$b.fq) / 4 ))
  echo "== $N q22 $b ($ba): $(grep -oE 'reads: certified [0-9.]+%' $O/q22_$b.log); uncertified $nu"
  for cfg in "l2:--local 2" "l2b:--local 2 --budget3 4096"; do
    n=${cfg%%:*}; a=${cfg#*:}
    /usr/bin/time -f "%e %U" -o $O/t_${b}_$n.txt $C map idx/grch38.s1q15.cidx $O/unc_$b.fq -t 16 $a \
      -o $O/q15_${b}_$n.sam --stats $O/q15_${b}_$n.json 2> $O/q15_${b}_$n.log
    read wall user < $O/t_${b}_$n.txt
    cert=$(grep -oE 'certified [0-9.]+%' $O/q15_${b}_$n.log | head -1)
    cand=$(grep -oE 'last pass\): [0-9.]+' $O/q15_${b}_$n.log | grep -oE '[0-9.]+$')
    ms=$(python3 -c "print(round($user / $nu * 1000, 3))")
    echo "   q15 $n ($a): $cert of residual; wall $wall s, CPU $user s = $ms ms/read; candidates/read $cand"
  done
done
