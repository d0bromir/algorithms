#!/bin/bash
# Memory changes (input pages dropped once parsed; the directory held as 32
# bits): every output byte against the build they start from, then peak
# memory and time at 2 M reads / pairs.
#   memory_identity.sh OLD_CERTA NEW_CERTA OUT_DIR [GPU_DEVICE]
set -u
OLD=$1; NEW=$2; D=$3; DEV=${4:-}
cd ~/certa-bench
mkdir -p $D
R1=reads/hg002_2M.fq
R2=reads/hg002_2M_2.fq
IDX=idx/grch38.s1v3.cidx
if [ "$(hostname -s)" = a2 ]; then P="numactl --interleave=all"; else P="taskset -c 0-63"; fi
n=0; d=0
for b in old new; do
  C=$OLD; [ $b = new ] && C=$NEW
  for x in "" "--full-cert"; do
    t=${x:+fc}; t=${t:-def}
    $P $C map $IDX $R1 -t 64 $x -o $D/se_$t.$b.sam -u $D/se_$t.$b.fq 2>$D/se_$t.$b.log
    $P $C map $IDX $R1 $R2 -t 64 $x -o $D/pe_$t.$b.sam -u $D/pe_$t.$b.u1 -U $D/pe_$t.$b.u2 2>$D/pe_$t.$b.log
    $P $C map $IDX $R1 -t 64 $x --intervals -o $D/sei_$t.$b.sam 2>$D/sei_$t.$b.log
    $P $C map $IDX $R1 $R2 -t 64 $x --intervals -o $D/pei_$t.$b.sam 2>$D/pei_$t.$b.log
    if [ -n "$DEV" ]; then
      $P $C map $IDX $R1 -t 64 $x --gpu --device $DEV --intervals -o $D/seig_$t.$b.sam 2>$D/seig_$t.$b.log
      $P $C map $IDX $R1 $R2 -t 64 $x --gpu --device $DEV --intervals -o $D/peig_$t.$b.sam 2>$D/peig_$t.$b.log
    fi
  done
done
for f in $D/*.old.*; do
  g=${f/.old./.new.}
  case $f in
    *.log)  # the summary lines, without timings
      diff -q <(grep -E "certified|pass V \(|pairs " $f) <(grep -E "certified|pass V \(|pairs " $g) > /dev/null ;;
    *) cmp -s $f $g ;;
  esac || { d=$((d + 1)); echo "DIFFERENT $(basename $f)"; }
  n=$((n + 1))
done
echo "new vs old: $((n - d)) of $n files identical (SAM, uncertified FASTQ and summary lines; default and --full-cert; single-end and paired; fallback and interval mode${DEV:+; interval mode on the GPU})"
rm -f $D/*.sam $D/*.fq $D/*.u1 $D/*.u2 $D/*.log
echo "== peak memory and time, interval mode, 64 threads, two rounds"
printf "build\tmode\tdevice\tconfig\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > $D/timings.tsv
for rep in 1 2; do
  for b in old new; do
    C=$OLD; [ $b = new ] && C=$NEW
    for x in default --full-cert; do
      y=""; [ $x != default ] && y=$x
      for m in "single-end:$R1" "paired:$R1 $R2"; do
        for dev in cpu ${DEV:+gpu}; do
          g=""; [ $dev = gpu ] && g="--gpu --device $DEV"
          /usr/bin/time -f "%e\t%U\t%S\t%M" -o $D/t $P $C map $IDX ${m#*:} -t 64 $y $g --intervals -o /dev/null 2>/dev/null
          printf "%s\t%s\t%s\t%s\t%s\t%s\n" $b ${m%%:*} $dev $x $rep "$(cat $D/t)" | tee -a $D/timings.tsv
        done
      done
    done
  done
done
rm -f $D/t
echo ALL_DONE
