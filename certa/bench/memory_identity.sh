#!/bin/bash
# Memory changes (input pages dropped once parsed; the directory held as 32
# bits): every output byte against the build they start from, then peak
# memory and time at 2 M reads / pairs. Outputs are hashed as they are
# written (SAM through a pipe, -u/-U through FIFOs): nothing large on disk.
#   memory_identity.sh OLD_CERTA NEW_CERTA OUT_DIR [GPU_DEVICE]
set -u
OLD=$1; NEW=$2; D=$3; DEV=${4:-}
cd ~/certa-bench
mkdir -p $D
R1=reads/hg002_2M.fq
R2=reads/hg002_2M_2.fq
IDX=idx/grch38.s1v3.cidx
if [ "$(hostname -s)" = a2 ]; then
  P="numactl --interleave=all"
  R1=/data/certa-bench/reads/hg002_2M.fq; R2=/data/certa-bench/reads/hg002_2M_2.fq
  [ -f $R1 ] || R1=reads/hg002_2M.fq
  [ -f $R2 ] || R2=reads/hg002_2M_2.fq
else
  P="taskset -c 0-63"
fi
fifo() {  # a FIFO whose contents are hashed into $1.md5
  rm -f $1; mkfifo $1
  md5sum < $1 | cut -c1-32 > $1.md5 &
}
run() {  # tag, certa map arguments (SAM to stdout)
  local t=$1; shift
  $P $C map $IDX "$@" 2> $D/$t.log | md5sum | cut -c1-32 > $D/$t.sam.md5
  wait
  rm -f $D/$t.u1 $D/$t.u2
}
for b in old new; do
  C=$OLD; [ $b = new ] && C=$NEW
  for x in "" "--full-cert"; do
    t=${x:+fc}; t=${t:-def}
    fifo $D/se_$t.$b.u1
    run se_$t.$b $R1 -t 64 $x -u $D/se_$t.$b.u1
    fifo $D/pe_$t.$b.u1; fifo $D/pe_$t.$b.u2
    run pe_$t.$b $R1 $R2 -t 64 $x -u $D/pe_$t.$b.u1 -U $D/pe_$t.$b.u2
    run sei_$t.$b $R1 -t 64 $x --intervals
    run pei_$t.$b $R1 $R2 -t 64 $x --intervals
    if [ -n "$DEV" ]; then
      run seig_$t.$b $R1 -t 64 $x --gpu --device $DEV --intervals
      run peig_$t.$b $R1 $R2 -t 64 $x --gpu --device $DEV --intervals
    fi
  done
done
n=0; d=0
for f in $D/*.old.*; do
  g=${f/.old./.new.}
  case $f in
    *.log)  # the summary lines, without timings
      diff -q <(grep -E "certified|pass V \(|pairs " $f) <(grep -E "certified|pass V \(|pairs " $g) > /dev/null ;;
    *) cmp -s $f $g ;;
  esac || { d=$((d + 1)); echo "DIFFERENT $(basename $f)"; }
  n=$((n + 1))
done
echo "new vs old: $((n - d)) of $n outputs identical (SAM, uncertified FASTQ and summary lines; default and --full-cert; single-end and paired; fallback and interval mode${DEV:+; interval mode on the GPU})"
grep -h "certified" $D/sei_*.new.log $D/pei_*.new.log
rm -f $D/*.md5 $D/*.log
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
