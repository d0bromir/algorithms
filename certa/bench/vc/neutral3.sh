#!/bin/bash
# galaxy: N3, a --mapq pair model from the locus and the mate only (an
# experiment): the read's own tier, MAPQ, XE, clusters, XQ, NM, clips and PR
# margin all depend on whether it shows the sample's variant (an alt read is
# never S0), so N3 may not split on them. Cross-validated on HG002 as
# neutral.sh's models.
set -uo pipefail
O=~/certa-trace/runs/neutral; cd $O
F=~/certa-trace/runs/cal; M=$F/models; TD=~/certa-trace/runs/t2t; W=~/certa-truth
P="python3 $O/truth_mapq_pair.py"
N3="tier,mq,xe,xn,xq,nm,clip,margin"
echo "== $(date +%T) fit N3"
taskset -c 125 $P fit --cv --drop $N3 $F/A.f4.npz $F/B.f4.npz > cv_N3.h 2> cv_N3.log
echo "== $(date +%T) scoring N3"
taskset -c 124 $P eval --half --shares 0.9413,0.8310 --model $M/cv_pair.h --model cv_N3.h \
  $W/A/truth.npz $F/A.f4.npz "minibwa=$TD/A/minibwa.npz" > cv3_A.md 2> cv3_A.log &
taskset -c 125 $P eval --half --shares 0.9975,0.9944,0.9558 --model $M/cv_pair.h --model cv_N3.h \
  $W/B/truth.npz $F/B.f4.npz "minibwa=$TD/B/minibwa.npz" > cv3_B.md 2> cv3_B.log &
wait
cat cv3_A.md cv3_B.md
echo "== $(date +%T) finished"
