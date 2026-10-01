#!/usr/bin/env bash
# Paired-end timing: BWA-MEM2, minibwa, and CERTA (GPU) with minibwa as the
# concurrent paired fallback (two named pipes; the fallback's records of mates
# CERTA already certified are dropped by drop_certified_mates.awk).
#   bench/run_bench_pe.sh IDX_DIR R1.fq R2.fq OUT_DIR
# Env: THREADS (64), REPS (3), GPU_DEVICE (1), PIPE_CERTA_THREADS (16),
#      CERTA (binary), CERTA_INDEX (grch38.s1v3.cidx), BIN, TOOLS, CERTA_ARGS.
set -euo pipefail
idx="$1"; r1="$2"; r2="$3"; out="$4"
T="${THREADS:-64}"; REPS="${REPS:-3}"; DEV="${GPU_DEVICE:-1}"; GT="${PIPE_CERTA_THREADS:-16}"
CERTA="${CERTA:-certa}"; CIDX="${CERTA_INDEX:-grch38.s1v3.cidx}"
BIN="${BIN:-$HOME/micromamba/envs/certa-bench/bin}"
TOOLS="${TOOLS:-certa-gpu-pe-pipe minibwa bwa-mem2}"
here="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$out"
tsv="$out/timings_pe.tsv"
[ -f "$tsv" ] || printf "tool\trep\tpairs\twall_s\tuser_s\tsys_s\tmax_rss_kb\n" > "$tsv"
pairs=$(( $(wc -l < "$r1") / 4 ))
for rep in $(seq 1 "$REPS"); do
  for tool in $TOOLS; do
    case "$tool" in
      certa-gpu-pe-pipe)
        rm -f "$out/u1" "$out/u2"; mkfifo "$out/u1" "$out/u2"
        CMD=(bash -c "set -o pipefail
          '$CERTA' map '$idx/$CIDX' '$r1' '$r2' -t $GT --gpu --device $DEV ${CERTA_ARGS:-} -o '$out/certa_pe.sam' \
            -u '$out/u1' -U '$out/u2' --stats '$out/certa_pe.stats.json' 2> '$out/certa_pe.log' &
          '$BIN/minibwa' map -t $((T - GT)) '$idx/minibwa' '$out/u1' '$out/u2' 2>/dev/null \
            | awk -f '$here/drop_certified_mates.awk' > '$out/fallback_pe.sam'
          wait") ;;
      minibwa)  CMD=(bash -c "'$BIN/minibwa' map -t $T '$idx/minibwa' '$r1' '$r2' > /dev/null 2>&1") ;;
      bwa-mem2) CMD=(bash -c "'$BIN/bwa-mem2' mem -t $T '$idx/bwa-mem2' '$r1' '$r2' > /dev/null 2>&1") ;;
      *) echo "unknown tool $tool" >&2; exit 1 ;;
    esac
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$out/time.txt" "${CMD[@]}"
    printf "%s\t%s\t%s\t%s\n" "$tool" "$rep" "$pairs" "$(cat "$out/time.txt")" >> "$tsv"
    tail -1 "$tsv"
  done
done
