#!/bin/bash
# truth_prq.sh NEW_BUILD OLD_BUILD OUT_DIR BEFORE_DIR: PR mates' point MAPQ
# from the best other pair found against the truth on HG002 tests A and B
# (bench/VALIDATION_RESULTS.md, "The MAPQ of mates placed by the pair").
# truth_rescue.sh's --iv-rescue local runs (interval mode, --iv-pair all,
# MAPQ cert and point) mapped by NEW_BUILD/certa, placed (truth_t2t.py
# place) and scored at truth_pair.sh's shares against BEFORE_DIR's runs of
# the build before (truth_rescue.py matched). Then --mapq pair on test B by
# both builds: the trees read the MAPQ they were fitted on, so the records
# must be the same (md5 of the SAM body). Aligner on cores 0-63.
set -euo pipefail
T=$(cd "$(dirname "$0")" && pwd)
C=$1/certa; C0=$2/certa; O=$(mkdir -p "$3" && cd "$3" && pwd); E=$(cd "$4" && pwd)
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai; IDX=~/certa-bench/idx/grch38.s1v3.cidx
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
for X in A B; do
  for q in cert point; do
    n=$X.$q
    if [ ! -s $O/$n.npz ]; then
      echo "== $(date +%T) $n"
      $PIN $TV -o $O/$n.time $C map $IDX $W/$X/${X}_1.fq $W/$X/${X}_2.fq --intervals --iv-mapq $q --iv-pair all \
        --iv-rescue local -t 64 -o /dev/stdout 2> $O/$n.log | $AUX python3 $T/truth_t2t.py place $W/$X/truth.npz $FAI $O/$n.npz
    fi
  done
  sh="0.9413,0.8310"; [ $X = B ] && sh="0.9975,0.9944,0.9558"
  for q in cert point; do
    { $AUX python3 $T/truth_rescue.py matched $W/$X/truth.npz $sh "before=$E/$X.$q.local.npz" \
        "PR MAPQ from the pair found=$O/$X.$q.npz"
      echo
      $AUX python3 $T/truth_rescue.py moved $W/$X/truth.npz $E/$X.$q.local.npz $O/$X.$q.npz; } > $O/$X.$q.md
  done
done
for b in old new; do
  c=$C; [ $b = old ] && c=$C0
  echo "== $(date +%T) B --mapq pair, $b build"
  $PIN $TV -o $O/B.pair.$b.time $c map $IDX $W/B/B_1.fq $W/B/B_2.fq --intervals --iv-pair all --mapq pair -t 64 \
    -o /dev/stdout 2> $O/B.pair.$b.log | grep -v "^@" | md5sum > $O/B.pair.$b.md5
done
cmp $O/B.pair.old.md5 $O/B.pair.new.md5 && echo "--mapq pair: the same records" || echo "--mapq pair: RECORDS DIFFER"
echo "== $(date +%T) done"
