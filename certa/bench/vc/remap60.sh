#!/bin/bash
# a2: remap60.sh GENOME NAME...: a diagnostic. CERTA's fitted MAPQ tops out
# near 43-47 (calibrated: 43 is one misplaced record in 20,000), while
# DeepVariant's model was trained on BWA-MEM, which gives unique reads 60.
# Copies NAME.cram with every MAPQ >= 30 set to 60 (the rest unchanged) as
# NAME60.cram, then runs call.sh on the copies once GENOME's own calls are done.
set -euo pipefail
G=$1; shift
V=~/vc; O=$V/$G; ST=$V/dv/rootfs/opt/conda/envs/bio/bin/samtools; REF=~/certa-bench/idx/ref.fa
while ! grep -q "^== .* done" $V/$G.call.log 2>/dev/null; do sleep 60; done
names=""
for n in "$@"; do
  echo "== $(date +%T) $n -> ${n}60"
  $ST view -@8 -h $O/$n.cram \
    | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} {if ($5 >= 30 && $5 != 255) $5 = 60; print}' \
    | $ST view -@8 -C -T $REF -o $O/${n}60.cram -
  $ST index $O/${n}60.cram
  cp $O/$n.time $O/${n}60.time
  names="$names ${n}60"
done
cd $V && scripts/call.sh $G $names
