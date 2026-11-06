#!/bin/bash
# truth_rescue.sh BUILD_DIR OUT_DIR [BEFORE_DIR]: --iv-rescue local against
# the truth on HG002 tests A and B (bench/VALIDATION_RESULTS.md, "A local
# mate rescue"). Interval mode with --iv-pair all, MAPQ cert (as the pairing
# fix was checked) and point (the default), each without the local rescue
# (e2e) and with it (local): placements (truth_t2t.py place), errors at
# matched N with truth_pair.sh's shares, and the mates the rescue moves
# (truth_rescue.py). With BEFORE_DIR, the e2e runs are read from there
# instead of mapped (the e2e path does not change). Aligner on cores 0-63.
set -euo pipefail
T=$(cd "$(dirname "$0")" && pwd)
C=$1/certa; O=$(mkdir -p "$2" && cd "$2" && pwd); E=${3:+$(cd "$3" && pwd)}; E=${E:-$O}
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai; IDX=~/certa-bench/idx/grch38.s1v3.cidx
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
for X in A B; do
  for q in cert point; do
    for m in e2e local; do
      n=$X.$q.$m
      [ $m = e2e ] && [ "$E" != "$O" ] && continue
      [ -s $O/$n.npz ] && continue
      echo "== $(date +%T) $n"
      $PIN $TV -o $O/$n.time $C map $IDX $W/$X/${X}_1.fq $W/$X/${X}_2.fq --intervals --iv-mapq $q --iv-pair all \
        --iv-rescue $m -t 64 -o /dev/stdout 2> $O/$n.log | $AUX python3 $T/truth_t2t.py place $W/$X/truth.npz $FAI $O/$n.npz
    done
  done
  sh="0.9413,0.8310"; [ $X = B ] && sh="0.9975,0.9944,0.9558"
  for q in cert point; do
    { $AUX python3 $T/truth_rescue.py matched $W/$X/truth.npz $sh "without (e2e)=$E/$X.$q.e2e.npz" \
        "--iv-rescue local=$O/$X.$q.local.npz"
      echo
      $AUX python3 $T/truth_rescue.py moved $W/$X/truth.npz $E/$X.$q.e2e.npz $O/$X.$q.local.npz; } > $O/$X.$q.md
  done
done
echo "== $(date +%T) done"
