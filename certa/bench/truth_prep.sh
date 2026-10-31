#!/bin/bash
# truth_prep.sh: phase 1 preparation (VALIDATION_PLAN.md): haplotype FASTAs and indexes,
# test A simulated pairs, test B trimmed real pairs. Main work on cores
# 64-127; trimming on 0-63.
set -eu
W=~/certa-truth; D=~/certa-data/truth; C=~/certa-bio/certa/build/certa
SIM=~/certa-bio/certa/tools/simulate.py
BIN=~/micromamba/envs/certa-bench/bin
export PATH=~/tools/env/truth/bin:$PATH
TV="/usr/bin/time -f %e\t%M"
mkdir -p $W/hap $W/A $W/B
cd $W
step() { echo "== $(date +%T) $*"; }

step haplotypes
for h in mat_Y_EBV_MT pat_X_EBV_MT; do
  [ -s hap/$h.fa ] || pigz -dc -p 8 $D/asm/hg002v1.1.$h.fasta.gz > hap/$h.fa
  samtools faidx hap/$h.fa
done
# Simulation: a male's genome once. Maternal file: maternal autosomes, X, and
# the paternal Y; paternal file: paternal autosomes only. No chrM / chrEBV.
awk '/^>/{keep = ($1 != ">chrM" && $1 != ">chrEBV")} keep' hap/mat_Y_EBV_MT.fa > hap/sim_mat.fa
awk '/^>/{keep = ($1 ~ /_PATERNAL$/ && $1 != ">chrY_PATERNAL")} keep' hap/pat_X_EBV_MT.fa > hap/sim_pat.fa
grep -c '>' hap/sim_mat.fa hap/sim_pat.fa

step "simulate test A (16 shards, 625k pairs each)"
i=0
for h in mat pat; do
  for s in 0 1 2 3 4 5 6 7; do
    taskset -c $((64 + i)) python3 $SIM reads hap/sim_$h.fa A/sh_${h}_$s.1.fq --mates A/sh_${h}_$s.2.fq \
      --n 625000 --quality novaseqx --snv 0 --indel 0 --insert 429,99 --seed $((1000 + i)) > A/sh_${h}_$s.log 2>&1 &
    i=$((i + 1))
  done
done

step "trim test B (20 M pairs)"
taskset -c 0-63 $TV -o B/trim.time ~/tools/bin/fastp_d0bromir_cpu_v5 \
  -i ~/certa-bench/reads/hg002_20M.fq -I ~/certa-bench/reads/hg002_20M_2.fq -o B/t_1.fq -O B/t_2.fq \
  -2 -g -Q -l 30 -w 32 -j B/trim.json -h B/trim.html > B/trim.log 2>&1
for m in 1 2; do awk 'NR % 4 == 1 {print "@b" (NR - 1) / 4; next} 1' B/t_$m.fq > B/B_$m.fq && rm B/t_$m.fq; done
echo "test B pairs: $(( $(wc -l < B/B_1.fq) / 4 ))"; cat B/trim.time

step "CERTA haplotype indexes"
for h in mat_Y_EBV_MT pat_X_EBV_MT; do
  [ -s hap/$h.cidx ] || taskset -c 80-127 $TV -o hap/$h.cidx.time $C index hap/$h.fa -o hap/$h.cidx -q 22 -s 1 -t 32 > hap/$h.cidx.log 2>&1
  tail -2 hap/$h.cidx.log; cat hap/$h.cidx.time
done

step "wait for the simulation"
wait
cat A/sh_*.log | sort | uniq -c
cat A/sh_mat_?.1.fq A/sh_pat_?.1.fq > A/orig_1.fq
cat A/sh_mat_?.2.fq A/sh_pat_?.2.fq > A/orig_2.fq
rm A/sh_*.fq
for m in 1 2; do awk 'NR % 4 == 1 {print "@a" (NR - 1) / 4; next} 1' A/orig_$m.fq > A/A_$m.fq; done
echo "test A pairs: $(( $(wc -l < A/A_1.fq) / 4 ))"

step "minibwa haplotype indexes"
for h in mat_Y_EBV_MT pat_X_EBV_MT; do
  [ -s hap/$h.mbw ] || taskset -c 64-127 $TV -o hap/$h.mbw.time $BIN/minibwa index -t 32 hap/$h.fa hap/$h > hap/$h.mbw.log 2>&1 || {
    echo "minibwa index failed:"; tail -5 hap/$h.mbw.log; }
  cat hap/$h.mbw.time 2>/dev/null
done
ls -la hap
df -h / | tail -1
step finished
