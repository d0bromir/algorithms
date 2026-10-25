#!/bin/bash
# Faster pass V (batched lookups, prefetched cluster windows): every output
# byte against the c9ee1c7 build, then the time of --full-cert.
set -u
cd /home/mpiuser/certa-bench
D=/home/mpiuser/certa-bench/cpu/pv2/out
mkdir -p $D
R1=reads/hg002_2M.fq
R2=reads/hg002_2M_2.fq
IDX=idx/grch38.s1v3.cidx
P="taskset -c 0-63"
NEW=/home/mpiuser/certa-fc/build/certa
OLD=/home/mpiuser/certa-bench/cpu/pv2/certa.c9ee1c7
n=0; d=0
for b in old new; do
  C=$OLD; [ $b = new ] && C=$NEW
  for x in "" "--full-cert"; do
    t=${x:+fc}; t=${t:-def}
    $P $C map $IDX $R1 -t 64 $x -o $D/se_$t.$b.sam -u $D/se_$t.$b.fq 2>/dev/null
    $P $C map $IDX $R1 $R2 -t 64 $x -o $D/pe_$t.$b.sam -u $D/pe_$t.$b.u1 -U $D/pe_$t.$b.u2 2>/dev/null
    $P $C map $IDX $R1 -t 64 $x --intervals -o $D/sei_$t.$b.sam 2>/dev/null
    $P $C map $IDX $R1 $R2 -t 64 $x --intervals -o $D/pei_$t.$b.sam 2>/dev/null
  done
  $P $C map $IDX $R1 $R2 -t 64 --full-cert --gpu --device 1 --intervals -o $D/peig_fc.$b.sam 2>/dev/null
done
for f in $D/*.old.*; do
  g=${f/.old./.new.}
  n=$((n + 1))
  cmp -s $f $g || { d=$((d + 1)); echo "DIFFERENT $(basename $f)"; }
done
echo "new vs c9ee1c7: $((n - d)) of $n files identical (default and --full-cert; SE/PE; fallback/interval; and PE --full-cert on the A100)"
rm -f $D/*
echo "== time, interval mode, 64 threads, two rounds"
for rep in 1 2; do
  for b in old new; do
    C=$OLD; [ $b = new ] && C=$NEW
    for m in "SE:$R1" "PE:$R1 $R2"; do
      /usr/bin/time -f "%e %U %S %M" -o $D/t $P $C map $IDX ${m#*:} -t 64 --full-cert --intervals -o /dev/null 2>/dev/null
      echo "$b ${m%%:*} --full-cert: wall/user/sys/maxrss $(cat $D/t)"
    done
  done
done
echo ALL_DONE
