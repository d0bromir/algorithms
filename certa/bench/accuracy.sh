#!/usr/bin/env bash
# Write SAM from every tool on the same reads, then compare with bwa-mem2.
#   bench/accuracy.sh IDX_DIR READS.fq OUT_DIR
# Environment: THREADS (64), CERTA (binary), BIN (tool dir).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
idx="$1"; reads="$2"; out="$3"
T="${THREADS:-64}"
BIN="${BIN:-$HOME/micromamba/envs/certa-bench/bin}"
mkdir -p "$out"
n=$(( $(wc -l < "$reads") / 4 ))
"$CERTA" map "$idx/grch38.cidx" "$reads" -t "$T" -o "$out/certa.sam" -u "$out/certa.uncertified.fq" 2> "$out/certa.log"
"$BIN/bwa-mem2" mem -t "$T" "$idx/bwa-mem2" "$reads" > "$out/bwa-mem2.sam" 2> "$out/bwa-mem2.log"
"$BIN/minibwa" map -t "$T" "$idx/minibwa" "$reads" > "$out/minibwa.sam" 2> "$out/minibwa.log"
"$BIN/strobealign" -t "$T" --use-index "$idx/ref.fa" "$reads" > "$out/strobealign.sam" 2> "$out/strobealign.log"
"$BIN/minimap2" -ax sr -t "$T" "$idx/minimap2.sr.mmi" "$reads" > "$out/minimap2.sam" 2> "$out/minimap2.log"
"$BIN/bowtie2" -p "$T" -x "$idx/bowtie2" -U "$reads" > "$out/bowtie2.sam" 2> "$out/bowtie2.log"
[[ -f "$idx/ref.fa.fai" ]] || "$BIN/samtools" faidx "$idx/ref.fa"
python3 "$here/compare_sam.py" --reads "$n" --certa "$out/certa.sam" --fasta "$idx/ref.fa" \
  bwa-mem2="$out/bwa-mem2.sam" minibwa="$out/minibwa.sam" strobealign="$out/strobealign.sam" \
  minimap2="$out/minimap2.sam" bowtie2="$out/bowtie2.sam" | tee "$out/accuracy.txt"
