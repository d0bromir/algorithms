#!/bin/bash
# CPU cost per read of each read class (class_*.fq from read_classes.py),
# default and --full-cert, interval mode, 64 threads, three rounds; a 1,000-
# read file gives the fixed cost (index load) to subtract.
#   class_cost.sh CERTA CLASS_DIR OUT_TSV
set -u
C=$1; D=$2; O=$3
cd ~/certa-bench
IDX=idx/grch38.s1v3.cidx
head -4000 $D/random.fq > $D/base.fq
printf "class\tconfig\trep\treads\twall_s\tcpu_s\n" > $O
for rep in 1 2 3; do
  for f in $D/*.fq; do
    c=$(basename $f .fq); n=$(( $(wc -l < $f) / 4 ))
    for x in default --full-cert; do
      y=""; [ $x != default ] && y=$x
      /usr/bin/time -f "%e %U %S" -o /tmp/certa_cc_t taskset -c 0-63 $C map $IDX $f -t 64 --intervals $y -o /dev/null 2>/dev/null
      read w u s < /tmp/certa_cc_t
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" $c $x $rep $n $w $(echo "$u + $s" | bc) >> $O
    done
  done
done
rm -f /tmp/certa_cc_t
echo ALL_DONE
