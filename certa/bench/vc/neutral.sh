#!/bin/bash
# galaxy: allele-neutral --mapq pair models (an experiment). DeepVariant calls
# fewer variants on CERTA's fitted MAPQ because the trees lower it for every
# edit, the sample's own variants included. Fitted here without the inputs
# that count the read's own edits (N1: XE, XQ, NM, soft clips, both mates;
# N2: also the record MAPQ), cross-validated on HG002 (even pairs; the odd
# ones scored at matched N) next to the committed inputs (cv_pair.h).
set -uo pipefail
O=~/certa-trace/runs/neutral; mkdir -p $O; cd $O
F=~/certa-trace/runs/cal; M=$F/models; TD=~/certa-trace/runs/t2t; W=~/certa-truth
P="python3 $O/truth_mapq_pair.py"
N1="xe,xq,nm,clip,p_xe,p_xq,p_nm,p_clip"; N2="$N1,mq,p_mq"
echo "== $(date +%T) fits"
taskset -c 126 $P fit --cv --drop $N1 $F/A.f4.npz $F/B.f4.npz > cv_N1.h 2> cv_N1.log &
taskset -c 127 $P fit --cv --drop $N2 $F/A.f4.npz $F/B.f4.npz > cv_N2.h 2> cv_N2.log &
wait
echo "== $(date +%T) scoring"
for X in A B; do
  sh="0.9413,0.8310"; core=126
  [ $X = B ] && { sh="0.9975,0.9944,0.9558"; core=127; }
  taskset -c $core $P eval --half --shares $sh --model $M/cv_pair.h --model cv_N1.h --model cv_N2.h \
    $W/$X/truth.npz $F/$X.f4.npz "minibwa=$TD/$X/minibwa.npz" > cv_$X.md 2> cv_$X.log &
done
wait
cat cv_A.md cv_B.md
echo "== $(date +%T) finished"
