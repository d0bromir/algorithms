#!/bin/bash
# Mutation tests of every bound (Aim 5): for each line of tests/mutants.tsv,
# a copy of the tree with one sed-applied mutation, built and run through the
# oracle test on two seeds: tests/test_certa.cpp, or, for a mutant of
# include/certa/longread.h, tests/test_longread.cpp and the exhaustive lemma
# checks (tests/verify_lemmas.cpp). A mutant expected "killed" must fail a
# check on some run; one expected "equivalent" (proven not to change any
# claim, THEORY.md) must pass them all.
#
#   tests/run_mutants.sh [parallel jobs = 8] [id regex = .] [seeds = "1 2"]
#
# Run from certa/. Prints one line per mutant and a summary; exits 1 if any
# mutant behaves unexpectedly.
set -u
here=$(cd "$(dirname "$0")/.." && pwd)
jobs=${1:-8}
filter=${2:-.}
seeds=${3:-1 2}
cxx=${CXX:-g++}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

run_one() {
  local id=$1 file=$2 expr=$3 expect=$4 what=$5
  local d=$work/$id
  mkdir -p "$d"
  cp -r "$here/include" "$here/src" "$here/tests" "$d/"
  ( cd "$d" || exit
    sed -i "$expr" "$file"
    if cmp -s "$file" "$here/$file"; then echo "$id	NOT_APPLIED	$what"; exit; fi
    local runs
    if [[ "$file" == include/certa/longread.h ]]; then
      if ! "$cxx" -O2 -std=c++17 -Iinclude -Isrc tests/test_longread.cpp src/index.cpp src/seqio.cpp -pthread -o t \
           2> build.log || ! "$cxx" -O2 -std=c++17 -Iinclude -Isrc tests/verify_lemmas.cpp -pthread -o v 2>> build.log; then
        echo "$id	BUILD_FAILED	$what"; exit
      fi
      runs="lemmas $seeds"
    else
      if ! "$cxx" -O2 -std=c++17 -Iinclude -Isrc tests/test_certa.cpp src/check.cpp src/gpu_stub.cpp \
           src/index.cpp src/seqio.cpp src/mapper_cpu.cpp -pthread -o t 2> build.log; then
        echo "$id	BUILD_FAILED	$what"; exit
      fi
      runs="$seeds"
    fi
    local res="" killed=0
    for s in $runs; do
      local r
      if [[ "$s" == lemmas ]]; then
        r=$(./v 2>&1 | grep -oE '[0-9]+ check\(s\) FAILED|all lemma checks passed' | tail -1)
      else
        r=$(./t "$s" 2>&1 | grep -oE '[0-9]+ check\(s\) FAILED|all certificate checks passed|all long-read checks passed' | tail -1)
      fi
      [[ "$r" == *FAILED* ]] && killed=1
      local label="seed $s"
      [[ "$s" == lemmas ]] && label=lemmas
      res="$res${res:+, }$label: ${r:-crashed}"
      [[ -z "$r" ]] && killed=1
    done
    local got=equivalent
    [[ $killed == 1 ]] && got=killed
    local ok=OK
    [[ "$got" != "$expect" ]] && ok=UNEXPECTED
    echo "$id	$ok	$got	$what	($res)" )
  rm -rf "$d"
}

n=0
while IFS=$'\t' read -r id file expr expect what; do
  [[ -z "$id" || "$id" == \#* ]] && continue
  [[ "$id" =~ $filter ]] || continue
  run_one "$id" "$file" "$expr" "$expect" "$what" &
  n=$((n + 1))
  if (( n % jobs == 0 )); then wait; fi
done < "$here/tests/mutants.tsv" | tee "$work/results.txt"
wait
total=$(wc -l < "$work/results.txt")
bad=$(grep -cvP '\tOK\t' "$work/results.txt")
echo "mutants: $total, unexpected: $bad"
[[ $bad == 0 ]]
