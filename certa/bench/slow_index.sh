#!/bin/bash
# Indexes for the slow but reliable mappers compared with CERTA's full
# certification (bench/slow_bench.sh): Winnowmap's repetitive 15-mers
# (meryl), Yara's GRCh38 index, and minimap2's HiFi index for
# CERTA's long-read pipeline; RazerS 3 needs none. Run on galaxy's cores
# 64-127 so 0-63 stay free. A step whose output exists is skipped.
#   slow_index.sh   (env certa-slow: winnowmap 2.03, yara 1.0.5, razers3 3.5.12;
#                    env certa-bench: minimap2 2.31)
set -u
E=$HOME/micromamba/envs/certa-slow/bin
M=$HOME/micromamba/envs/certa-bench/bin/minimap2
D=$HOME/certa-bench/idx/slow
REF=$HOME/certa-bench/idx/ref.fa
P="taskset -c 64-127 nice"
mkdir -p "$D"
cd "$D" || exit 1
if [ ! -s repetitive_k15.txt ]; then
  echo "== $(date +%T) meryl count"
  /usr/bin/time -v $P "$E/meryl" count k=15 threads=64 memory=100 output merylDB_k15 "$REF" > meryl.log 2>&1
  $P "$E/meryl" print greater-than distinct=0.9998 merylDB_k15 > repetitive_k15.txt 2>> meryl.log
  echo "== $(date +%T) meryl done: $(wc -l < repetitive_k15.txt) repetitive 15-mers"
fi
if [ ! -s yara_grch38.lf.drv ]; then
  echo "== $(date +%T) yara_indexer"
  /usr/bin/time -v $P "$E/yara_indexer" "$REF" -o yara_grch38 > yara_index.log 2>&1
  echo "== $(date +%T) yara_indexer exit $?"
fi
# Winnowmap 2.03 indexes the FASTA in every run: its dumped index (-d) builds
# in 90 s but crashes on load (segmentation fault, aarch64 Bioconda build).
if [ ! -s mm2_hifi.mmi ]; then
  echo "== $(date +%T) minimap2 index (map-hifi)"
  /usr/bin/time -v $P "$M" -x map-hifi -t 16 -d mm2_hifi.mmi "$REF" > mm2_index.log 2>&1
  echo "== $(date +%T) minimap2 index exit $?"
fi
ls -la "$D"
echo ALL_DONE
