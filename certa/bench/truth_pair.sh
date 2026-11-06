#!/bin/bash
# truth_pair.sh: the pair-aware MAPQ (--mapq pair) and T2T-CHM13 in the same
# pass (--t2t), fitted on HG002 tests A and B and held out on HG005
# (bench/VALIDATION_RESULTS.md, "A MAPQ from both mates"). Aligners on cores
# 0-63, the rest on 64-127.
#   truth_pair.sh features BUILD_DIR HS1_CIDX CHAIN OUT_DIR
#   truth_pair.sh fit OUT_DIR T2T_DIR
#   truth_pair.sh check BUILD_DIR HS1_CIDX CHAIN OUT_DIR
#   truth_pair.sh holdout BUILD_DIR HS1_CIDX CHAIN OUT_DIR H_PLACE_DIR
# HS1_CIDX is `certa index hs1.fa -q 22 -s 1`; CHAIN is UCSC's
# hs1ToHg38.over.chain.gz; T2T_DIR is truth_t2t.sh's OUT_DIR (the other
# aligners' placements, X/*.npz); H_PLACE_DIR is truth_hg005.sh's
# OUT_DIR/place.
# - features: tests A, B and H with --t2t --feature-tags (X.f4.npz);
# - fit: the two models, cross-validated on HG002 (even pairs, then the odd
#   pairs at matched N: cv_A.md, cv_B.md) and on every pair
#   (mapq_model.h, mapq_model_t2t.h, to copy to include/certa/);
# - check (a build with those headers): --mapq pair, with and without
#   --t2t, on tests A and B: CERTA's MAPQ against the fit's replay;
# - holdout (the same build, run once): HG005 with --mapq pair, with and
#   without --t2t, against the other aligners and BWA-MEM2's T2T veto
#   (holdout_H.md).
set -euo pipefail
step=$1; shift
T=$(cd "$(dirname "$0")" && pwd)
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai; IDX=~/certa-bench/idx/grch38.s1v3.cidx
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
IV="--intervals --iv-mapq cert --iv-pair all -t 64 -o /dev/stdout"
P="$AUX python3 $T/truth_mapq_pair.py"
reads() { echo "$W/$1/${1}_1.fq $W/$1/${1}_2.fq"; }
case $step in
features)
  C=$1/certa; T2T="--t2t $2 --t2t-chain $3"; O=$(mkdir -p "$4" && cd "$4" && pwd)
  for X in A B H; do
    echo "== $(date +%T) $X: --t2t --feature-tags"
    [ -s $O/$X.f4.npz ] || $PIN $TV -o $O/$X.f4.time $C map $IDX $(reads $X) $IV --mapq-limit 10 --feature-tags $T2T \
      2> $O/$X.f4.log | $P feat $W/$X/truth.npz $FAI $O/$X.f4.npz
  done ;;
fit)
  O=$(cd "$1" && pwd); TD=$(cd "$2" && pwd)
  # Depth 4 (300 trees) and 6 (400) cross-validated; the rule fixed beforehand
  # (depth 6 if it lowers test B's errors by >= 5 % at two of three N and
  # raises none by > 5 %) chose 6 for the T2T model only, as fitted below.
  echo "== $(date +%T) fits"
  D6="--depth 6 --rounds 400"
  $P fit --cv $O/A.f4.npz $O/B.f4.npz > $O/cv_pair.h 2> $O/cv_pair.log &
  $P fit --cv $D6 $O/A.f4.npz $O/B.f4.npz > $O/cv_pair_d6.h 2> $O/cv_pair_d6.log &
  $P fit --cv --t2t $O/A.f4.npz $O/B.f4.npz > $O/cv_t2t.h 2> $O/cv_t2t.log &
  $P fit --cv --t2t $D6 $O/A.f4.npz $O/B.f4.npz > $O/cv_t2t_d6.h 2> $O/cv_t2t_d6.log &
  $P fit $O/A.f4.npz $O/B.f4.npz > $O/mapq_model.h 2> $O/mapq_model.log &
  $P fit --t2t $D6 $O/A.f4.npz $O/B.f4.npz > $O/mapq_model_t2t.h 2> $O/mapq_model_t2t.log &
  wait
  for X in A B; do
    sh="0.9413,0.8310"; [ $X = B ] && sh="0.9975,0.9944,0.9558"
    $P eval --half --shares $sh --model $O/cv_pair.h --model $O/cv_t2t.h --model $O/cv_pair_d6.h \
      --model $O/cv_t2t_d6.h $W/$X/truth.npz $O/$X.f4.npz \
      "certa --mapq cal=$TD/$X/ivcal.npz" "minibwa=$TD/$X/minibwa.npz" "bwa-mem2=$TD/$X/bwa-mem2.npz" \
      "strobealign=$TD/$X/strobealign.npz" "minimap2=$TD/$X/minimap2.npz" > $O/cv_$X.md 2> $O/cv_$X.log &
  done
  wait ;;
check)
  C=$1/certa; T2T="--t2t $2 --t2t-chain $3"; O=$(cd "$4" && pwd)
  for X in A B; do
    echo "== $(date +%T) $X: replay"
    $PIN $C map $IDX $(reads $X) $IV --mapq pair 2> $O/$X.pair.log \
      | $P replay $T/../include/certa/mapq_model.h $O/$X.f4.npz | tee $O/$X.pair.replay
    $PIN $C map $IDX $(reads $X) $IV --mapq pair $T2T 2> $O/$X.pairt2t.log \
      | $P replay $T/../include/certa/mapq_model_t2t.h $O/$X.f4.npz | tee $O/$X.pairt2t.replay
  done ;;
holdout)
  C=$1/certa; T2T="--t2t $2 --t2t-chain $3"; O=$(cd "$4" && pwd); HP=$(cd "$5" && pwd)
  for run in pair pairt2t; do
    extra=""; [ $run = pairt2t ] && extra=$T2T
    echo "== $(date +%T) H: --mapq pair $extra"
    if [ ! -s $O/H.$run.npz ]; then
      rm -f $O/H.$run.placed
      $PIN $TV -o $O/H.$run.time $C map $IDX $(reads H) $IV --mapq pair $extra 2> $O/H.$run.log \
        | tee >($AUX python3 $T/truth_t2t.py place $W/H/truth.npz $FAI $O/H.$run.npz; touch $O/H.$run.placed) \
        | $P replay $T/../include/certa/mapq_model$([ $run = pairt2t ] && echo _t2t).h $O/H.f4.npz | tee $O/H.$run.replay
      while [ ! -e $O/H.$run.placed ]; do sleep 5; done  # bash does not wait for the process substitution
    fi
  done
  for a in minibwa bwa-mem2; do
    $AUX python3 $T/truth_veto.py $HP/$a.npz $HP/t2t.npz $O/H.$a.veto.npz
  done
  $P eval --shares 0.9975,0.9929,0.98,0.95 $W/H/truth.npz $O/H.f4.npz \
    "certa --mapq cal=$HP/ivcal.npz" "certa --mapq pair=$O/H.pair.npz" "certa --mapq pair --t2t=$O/H.pairt2t.npz" \
    "minibwa=$HP/minibwa.npz" "bwa-mem2=$HP/bwa-mem2.npz" "strobealign=$HP/strobealign.npz" \
    "minimap2=$HP/minimap2.npz" "minibwa, BWA-MEM2's T2T veto=$O/H.minibwa.veto.npz" \
    "bwa-mem2, BWA-MEM2's T2T veto=$O/H.bwa-mem2.veto.npz" > $O/holdout_H.md ;;
*) echo "unknown step $step" >&2; exit 2 ;;
esac
echo "== $(date +%T) $step finished"
