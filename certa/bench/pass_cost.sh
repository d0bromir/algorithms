#!/bin/bash
# CPU per read of each pass, by differences between configurations, on the
# class files given (100,000 reads each): fallback mode (no tier R) and
# interval mode, default, --full-cert, and --full-cert without pass V or L.
#   pass_cost.sh CERTA OUT_TSV CLASS.fq...
set -u
C=$1; O=$2; shift 2
cd ~/certa-bench
IDX=idx/grch38.s1v3.cidx
printf "class\tconfig\trep\treads\twall_s\tcpu_s\n" > $O
CFG=("" "--intervals" "--full-cert" "--full-cert --intervals" "--full-cert --adaptive-parts 0 --intervals"
     "--full-cert --local 0 --intervals" "--full-cert --budget3 0 --intervals")
for rep in 1 2 3; do
  for f in "$@"; do
    c=$(basename $f .fq); n=$(( $(wc -l < $f) / 4 ))
    for x in "${CFG[@]}"; do
      /usr/bin/time -f "%e %U %S" -o /tmp/certa_pc_t taskset -c 0-63 $C map $IDX $f -t 64 $x -o /dev/null 2>/dev/null
      read w u s < /tmp/certa_pc_t
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" $c "${x:-fallback}" $rep $n $w $(echo "$u + $s" | bc) >> $O
    done
  done
done
rm -f /tmp/certa_pc_t
echo ALL_DONE
