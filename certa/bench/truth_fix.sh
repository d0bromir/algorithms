#!/bin/bash
# galaxy: truth_fix.sh SRC_TGZ NAME TAG [EXTRA...]: builds ~/certa-NAME from the
# sources, then HG002 tests A and B (cert and point MAPQ; --iv-pair all
# --iv-rescue local EXTRA) placed by truth_t2t.py and scored at the locus
# against ~/certa-mr3's runs (PR 20's build): truth_rescue.py matched and moved.
# Results in ~/certa-NAME/acc/X.q.TAG.*. Aligner on cores 0-63.
set -euo pipefail
TGZ=$1; NAME=$2; TAG=$3; shift 3
P=~/certa-$NAME; C=$P/build/certa; T=$(cd "$(dirname "$0")" && pwd); E=~/certa-mr3/acc; O=$P/acc
W=~/certa-truth; FAI=~/certa-bench/idx/ref.fa.fai; IDX=~/certa-bench/idx/grch38.s1v3.cidx
export PATH=~/tools/env/truth/bin:$PATH
PIN="taskset -c 0-63"; AUX="taskset -c 64-127"
TV="/usr/bin/time -f %e\t%U\t%S\t%M"
if [ -n "$TGZ" ]; then
  rm -rf $P/include $P/src; mkdir -p $P/build/obj $O; tar xzf $TGZ -C $P
  echo "build $(cd $P && find include src -type f | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -c1-16)"
  for f in index seqio mapper_cpu check chain main cmd_map map_config map_pipeline map_stats sam_output gpu_stub; do
    echo "g++ -O3 -DNDEBUG -mcpu=native -std=c++17 -I$P/include -I$P/src -c $P/src/$f.cpp -o $P/build/obj/$f.o"
  done | xargs -P 12 -I{} sh -c {}
  g++ -O3 -mcpu=native $P/build/obj/*.o -lpthread -o $C
fi
for X in A B; do
  for q in cert point; do
    n=$X.$q.$TAG
    if [ ! -s $O/$n.npz ]; then
      echo "== $(date +%T) $n"
      $PIN $TV -o $O/$n.time $C map $IDX $W/$X/${X}_1.fq $W/$X/${X}_2.fq --intervals --iv-mapq $q --iv-pair all \
        --iv-rescue local "$@" -t 64 -o /dev/stdout 2> $O/$n.log | $AUX python3 $T/truth_t2t.py place $W/$X/truth.npz $FAI $O/$n.npz
    fi
  done
  sh="0.9413,0.8310"; [ $X = B ] && sh="0.9975,0.9944,0.9558"
  for q in cert point; do
    { $AUX python3 $T/truth_rescue.py matched $W/$X/truth.npz $sh "before (PR 20)=$E/$X.$q.local.npz" "$TAG=$O/$X.$q.$TAG.npz"
      echo
      $AUX python3 $T/truth_rescue.py moved $W/$X/truth.npz $E/$X.$q.local.npz $O/$X.$q.$TAG.npz; } > $O/$X.$q.$TAG.md
  done
done
echo "== $(date +%T) done"
