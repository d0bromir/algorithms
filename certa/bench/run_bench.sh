#!/usr/bin/env bash
# Timed single-end mapping benchmark: CERTA (alone and CERTA + fallback)
# versus the fastest short-read mappers, same reads, same threads.
#   bench/run_bench.sh IDX_DIR READS.fq OUT_DIR
# Environment: THREADS (64), REPS (3), CERTA (binary), BIN (tool dir),
#   GPU_DEVICE (unset = no GPU runs), FALLBACK (minibwa), TOOLS (list).
# Every tool writes SAM to /dev/null for timing (formatting is included,
# disk writes are not). Index files are warmed into the page cache first.
set -euo pipefail
idx="$1"; reads="$2"; out="$3"
T="${THREADS:-64}"; REPS="${REPS:-3}"; FALLBACK="${FALLBACK:-minibwa}"
BIN="${BIN:-$HOME/micromamba/envs/certa-bench/bin}"
TOOLS="${TOOLS:-certa certa+fallback minibwa strobealign bwa-mem2 minimap2 bowtie2}"
[[ -n "${GPU_DEVICE:-}" ]] && TOOLS="$TOOLS certa-gpu certa-gpu+fallback"
mkdir -p "$out"
tsv="$out/timings.tsv"
[[ -s "$tsv" ]] || printf "host\ttool\trep\treads\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$tsv"
nreads=$(( $(wc -l < "$reads") / 4 ))

cmd_for() {  # tool, reads file -> command array in CMD
  case "$1" in
    certa)       CMD=("$CERTA" map "$idx/grch38.cidx" "$2" -t "$T" -o /dev/null -u "$out/uncertified.fq" --stats "$out/certa.stats.json") ;;
    certa-gpu)   CMD=("$CERTA" map "$idx/grch38.cidx" "$2" -t "$T" --gpu --device "$GPU_DEVICE" -o /dev/null -u "$out/uncertified.fq" --stats "$out/certa-gpu.stats.json") ;;
    minibwa)     CMD=("$BIN/minibwa" map -t "$T" "$idx/minibwa" "$2") ;;
    strobealign) CMD=("$BIN/strobealign" -t "$T" --use-index "$idx/ref.fa" "$2") ;;
    bwa-mem2)    CMD=("$BIN/bwa-mem2" mem -t "$T" "$idx/bwa-mem2" "$2") ;;
    minimap2)    CMD=("$BIN/minimap2" -ax sr -t "$T" "$idx/minimap2.sr.mmi" "$2") ;;
    bowtie2)     CMD=("$BIN/bowtie2" -p "$T" -x "$idx/bowtie2" -U "$2") ;;
    *) echo "unknown tool $1" >&2; exit 1 ;;
  esac
}

run_timed() {  # label rep tool reads-file
  local label="$1" rep="$2" tool="$3" file="$4" tf="$out/.time"
  cmd_for "$tool" "$file"
  /usr/bin/time -v -o "$tf" "${CMD[@]}" > /dev/null 2> "$out/$label.$rep.log"
  local wall user sys rss
  wall=$(awk -F': ' '/Elapsed/ { n = split($2, a, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + a[i]; print s }' "$tf")
  user=$(awk -F': ' '/User time/ { print $2 }' "$tf")
  sys=$(awk -F': ' '/System time/ { print $2 }' "$tf")
  rss=$(awk -F': ' '/Maximum resident/ { print $2 }' "$tf")
  echo "$wall $user $sys $rss"
}

head -n 400000 "$reads" > "$out/warmup.fq"
for tool in $TOOLS; do
  base="${tool%+fallback}"
  echo "== $(date +%T) $tool: warm-up"
  run_timed "$base.warmup" 0 "$base" "$out/warmup.fq" > /dev/null
  for rep in $(seq 1 "$REPS"); do
    read -r wall user sys rss < <(run_timed "$tool" "$rep" "$base" "$reads")
    if [[ "$tool" == *+fallback ]]; then  # map the uncertified reads with the fallback
      read -r w2 u2 s2 r2 < <(run_timed "$tool.fallback" "$rep" "$FALLBACK" "$out/uncertified.fq")
      wall=$(awk -v a="$wall" -v b="$w2" 'BEGIN { print a + b }')
      user=$(awk -v a="$user" -v b="$u2" 'BEGIN { print a + b }')
      sys=$(awk -v a="$sys" -v b="$s2" 'BEGIN { print a + b }')
      rss=$(( rss > r2 ? rss : r2 ))
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(hostname -s)" "$tool" "$rep" "$nreads" "$wall" "$user" "$sys" "$rss" | tee -a "$tsv"
  done
done
rm -f "$out/.time" "$out/warmup.fq"
