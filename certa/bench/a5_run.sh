#!/bin/bash
# A5: fit the phantom prior (delta, pi) on HG002 by calibration error, freeze,
# validate on the held-out HG005. Disk-safe: everything in /dev/shm (RAM), SAM
# streamed (no BAM). Fetches the a5 branch (which carries fit_a5.py).
set -uo pipefail
mkdir -p /dev/shm/a2/tmp
export TMPDIR=/dev/shm/a2/tmp
cd /dev/shm/a2/src
git fetch --depth 1 origin claude/a5-calibration > /dev/shm/a5_pull.log 2>&1 \
  && git reset --hard FETCH_HEAD >> /dev/shm/a5_pull.log 2>&1 \
  || { echo PULL_FAIL; cat /dev/shm/a5_pull.log; exit 1; }
cd certa
cmake --build build -j 32 --target certa > /dev/shm/a5_cli.log 2>&1 \
  || { echo CLI_BUILD_FAIL; tail -30 /dev/shm/a5_cli.log; exit 1; }
echo CLI_OK
C=/dev/shm/a2/src/certa/build/certa
B=~/certa-bench; W=~/certa-truth/B; H=~/certa-truth/H; T=/dev/shm/a2/src/certa/bench
IDX=$B/idx/grch38.s1v3.cidx; FAI=$B/idx/ref.fa.fai
cd /dev/shm/a2
[ -f segdups.bed ] || zcat ~/certa-data/truth/strat/SegmentalDuplications/GRCh38_segdups.bed.gz > segdups.bed
for D in 0.003 0.006 0.01; do
  for P in 0.003 0.01 0.03; do
    tag="fit_d${D}_p${P}"
    echo "== $(date +%T) HG002 $tag"
    taskset -c 0-63 "$C" map "$IDX" sub_1.fq --bio --bio-phantom "$D" --bio-phantom-bed segdups.bed \
      --bio-phantom-pi "$P" --bio-phantom-pi-unique 0 -t 64 -o /dev/stdout -u /dev/null 2>"${tag}.log" \
      | python3 "$T/truth_eval.py" score "$W/truth.npz" "$FAI" "${tag}.npz" 2>"${tag}.eval.log"
  done
done
echo "== $(date +%T) fit (pick best by calibration error, freeze)"
python3 "$T/fit_a5.py" "$W/truth.npz" /dev/shm/a2
D=$(python3 -c "import json;print(json.load(open('/dev/shm/a2/a5_params.json'))['delta'])")
P=$(python3 -c "import json;print(json.load(open('/dev/shm/a2/a5_params.json'))['pi'])")
head -n 16000000 "$H/H_1.fq" > /dev/shm/a2/hsub_1.fq  # 4M HG005 reads (mate 0)
echo "== $(date +%T) HG005 baseline"
taskset -c 0-63 "$C" map "$IDX" hsub_1.fq --bio -t 64 -o /dev/stdout -u /dev/null 2>h0.log \
  | python3 "$T/truth_eval.py" score "$H/truth.npz" "$FAI" h0.npz 2>h0.eval.log
echo "== $(date +%T) HG005 frozen (delta=$D pi=$P)"
taskset -c 0-63 "$C" map "$IDX" hsub_1.fq --bio --bio-phantom "$D" --bio-phantom-bed segdups.bed \
  --bio-phantom-pi "$P" --bio-phantom-pi-unique 0 -t 64 -o /dev/stdout -u /dev/null 2>h1.log \
  | python3 "$T/truth_eval.py" score "$H/truth.npz" "$FAI" h1.npz 2>h1.eval.log
python3 "$T/truth_eval.py" report "$H/truth.npz" hg005-bio=h0.npz hg005-phantom=h1.npz > a5_holdout.md
echo DONE
echo ---FIT---; cat a5_fit.md
echo ---HOLDOUT---; sed -n '1,16p' a5_holdout.md
