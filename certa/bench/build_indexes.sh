#!/usr/bin/env bash
# Build every tool's index for one reference, timing each build.
#   bench/build_indexes.sh REF.fa IDX_DIR THREADS [tools...]
# Tools: certa minibwa bwa-mem2 strobealign minimap2 bowtie2 (default: all).
# Environment: CERTA (certa binary), BIN (directory with the other tools).
set -euo pipefail
ref="$1"; idx="$2"; t="$3"; shift 3
tools=("$@"); [[ ${#tools[@]} -eq 0 ]] && tools=(certa minibwa bwa-mem2 strobealign minimap2 bowtie2)
BIN="${BIN:-$HOME/micromamba/envs/certa-bench/bin}"
mkdir -p "$idx"
ln -sf "$(readlink -f "$ref")" "$idx/ref.fa"  # strobealign writes its index next to the FASTA
timed() { local name="$1"; shift; /usr/bin/time -v -o "$idx/$name.index.time" "$@" > "$idx/$name.index.log" 2>&1; }
for tool in "${tools[@]}"; do
  echo "== $(date +%T) building $tool index"
  case "$tool" in
    certa)       timed certa "$CERTA" index "$idx/ref.fa" -o "$idx/grch38.cidx" -t "$t" ;;
    minibwa)     timed minibwa "$BIN/minibwa" index -t "$t" "$idx/ref.fa" "$idx/minibwa" ;;
    bwa-mem2)    timed bwa-mem2 "$BIN/bwa-mem2" index -p "$idx/bwa-mem2" "$idx/ref.fa" ;;
    strobealign) timed strobealign "$BIN/strobealign" -t "$t" --create-index -r 150 "$idx/ref.fa" ;;
    minimap2)    timed minimap2 "$BIN/minimap2" -x sr -t "$t" -d "$idx/minimap2.sr.mmi" "$idx/ref.fa" ;;
    bowtie2)     timed bowtie2 "$BIN/bowtie2-build" --threads "$t" "$idx/ref.fa" "$idx/bowtie2" ;;
    *) echo "unknown tool $tool" >&2; exit 1 ;;
  esac
  grep -E "Elapsed|Maximum resident" "$idx/$tool.index.time" | sed 's/^\s*/   /'
done
