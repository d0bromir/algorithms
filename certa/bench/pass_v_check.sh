#!/bin/bash
# Pass V (adaptive parts): default output against 77b1df2, pass V off against
# 77b1df2's equivalent flags, CPU = GPU and serial = fused with pass V on;
# then certification and time on 2 M HG002 reads and pairs.
set -u
cd /home/mpiuser/certa-bench
D=/home/mpiuser/certa-bench/cpu/pv/out
mkdir -p $D
R1=reads/hg002_2M.fq
R2=reads/hg002_2M_2.fq
IDX=idx/grch38.s1v3.cidx
IDX2=idx/grch38.s1q15.cidx
P="taskset -c 0-63"
NEW=/home/mpiuser/certa-fc/build/certa
OLD=/home/mpiuser/certa-final/build/certa
G="--gpu --device 1"
nopg() { grep -v "^@PG" "$1"; }
same() {  # label, files a b (SAM compared without @PG; others byte for byte)
  local lab=$1 a=$2 b=$3
  if [[ "$a" == *.sam ]]; then cmp -s <(nopg $a) <(nopg $b) && echo "identical (but @PG): $lab" || echo "DIFFERENT: $lab"
  else cmp -s $a $b && echo "identical: $lab" || echo "DIFFERENT: $lab"; fi
}
echo "== $(date +%T) default output: new vs 77b1df2 (byte for byte, @PG included)"
n=0; d=0
for dev in cpu gpu; do
  X=""; [ $dev = gpu ] && X=$G
  for b in old new; do
    C=$OLD; [ $b = new ] && C=$NEW
    $P $C map $IDX $R1 -t 64 $X -o $D/se_$dev.$b.sam -u $D/se_$dev.$b.fq 2>/dev/null
    $P $C map $IDX $R1 $R2 -t 64 $X -o $D/pe_$dev.$b.sam -u $D/pe_$dev.$b.u1 -U $D/pe_$dev.$b.u2 2>/dev/null
    $P $C map $IDX $R1 -t 64 $X --intervals -o $D/sei_$dev.$b.sam 2>/dev/null
    $P $C map $IDX $R1 $R2 -t 64 $X --intervals --witness -o $D/pei_$dev.$b.sam 2>/dev/null
  done
  for f in $D/*_$dev.old.*; do
    g=${f/.old./.new.}
    n=$((n + 1))
    cmp -s $f $g || { d=$((d + 1)); echo "DIFFERENT $(basename $f)"; }
  done
done
echo "default output: $((n - d)) of $n files identical"
rm -f $D/*.old.* $D/*.new.*
echo "== $(date +%T) pass V off: new --full-cert --adaptive-parts 0 vs 77b1df2 --budget3 4096 --local 2"
$P $OLD map $IDX $R1 $R2 -t 64 --budget3 4096 --local 2 -o $D/off.old.sam -u $D/off.old.u1 -U $D/off.old.u2 2>/dev/null
$P $NEW map $IDX $R1 $R2 -t 64 --full-cert --adaptive-parts 0 -o $D/off.new.sam -u $D/off.new.u1 -U $D/off.new.u2 2>/dev/null
same "paired fallback SAM" $D/off.old.sam $D/off.new.sam
same "paired unmapped mates 1" $D/off.old.u1 $D/off.new.u1
same "paired unmapped mates 2" $D/off.old.u2 $D/off.new.u2
$P $OLD map $IDX $R1 -t 64 --budget3 4096 --local 2 --intervals -o $D/offi.old.sam 2>/dev/null
$P $NEW map $IDX $R1 -t 64 --full-cert --adaptive-parts 0 --intervals -o $D/offi.new.sam 2>/dev/null
same "single-end interval SAM" $D/offi.old.sam $D/offi.new.sam
rm -f $D/off*
echo "== $(date +%T) pass V on: CPU = GPU, fused = serial stages"
$P $NEW map $IDX $R1 $R2 -t 64 --full-cert -o $D/v_cpu.sam -u $D/v_cpu.u1 -U $D/v_cpu.u2 2>/dev/null
$P $NEW map $IDX $R1 $R2 -t 64 --full-cert $G -o $D/v_gpu.sam -u $D/v_gpu.u1 -U $D/v_gpu.u2 2>/dev/null
same "paired --full-cert CPU vs GPU" $D/v_cpu.sam $D/v_gpu.sam
same "paired --full-cert CPU vs GPU, mates 1" $D/v_cpu.u1 $D/v_gpu.u1
$P $NEW map $IDX $R1 -t 64 --full-cert -o $D/vs_f.sam -u $D/vs_f.fq 2>/dev/null
$P $NEW map $IDX $R1 -t 64 --full-cert --serial-stages -o $D/vs_s.sam -u $D/vs_s.fq 2>/dev/null
same "single-end --full-cert fused vs serial" $D/vs_f.sam $D/vs_s.sam
same "single-end --full-cert fused vs serial, unmapped" $D/vs_f.fq $D/vs_s.fq
$P $NEW map $IDX $R1 -t 64 --full-cert --intervals -o $D/vi_cpu.sam 2>/dev/null
$P $NEW map $IDX $R1 -t 64 --full-cert --intervals $G -o $D/vi_gpu.sam 2>/dev/null
same "single-end --full-cert interval CPU vs GPU" $D/vi_cpu.sam $D/vi_gpu.sam
rm -f $D/v*
echo "== $(date +%T) certification and time, 2 M reads / pairs, interval mode, 64 threads"
run() {  # label, certa arguments
  local label=$1
  shift
  /usr/bin/time -f "%e %U %S %M" -o $D/t $P $NEW map $IDX "$@" -t 64 --intervals -o /dev/null 2> $D/log
  echo "-- $label: wall/user/sys/maxrss $(cat $D/t)"
  grep -E "reads: certified|pairs .*both mates|pass V|second index" $D/log
}
for rep in 1 2; do
  run "SE --full-cert --adaptive-parts 0" $R1 --full-cert --adaptive-parts 0
  run "SE --full-cert --adaptive-parts 1" $R1 --full-cert --adaptive-parts 1
  run "SE --full-cert" $R1 --full-cert
  run "SE --full-cert --index2" $R1 --full-cert --index2 $IDX2
  run "SE --full-cert --adaptive-parts 0 --index2" $R1 --full-cert --adaptive-parts 0 --index2 $IDX2
  run "PE --full-cert --adaptive-parts 0" $R1 $R2 --full-cert --adaptive-parts 0
  run "PE --full-cert" $R1 $R2 --full-cert
  run "PE --full-cert --index2" $R1 $R2 --full-cert --index2 $IDX2
  run "SE --full-cert, A100" $R1 --full-cert $G
  run "PE --full-cert, A100" $R1 $R2 --full-cert $G
done
echo ALL_DONE
