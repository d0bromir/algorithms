#!/bin/bash
# a2: lost_run.sh: HG002's lost reads (minibwa MAPQ >= 20 in the targets; CERTA
# > 1 kb away, another contig or unmapped, by the --t2t run's placements),
# both mates re-mapped by the same build with --feature-tags --witness, and
# minibwa's records of them; then lost_diag.py.
set -euo pipefail
D=~/vc/diag/lost; mkdir -p $D; cd $D
ST=~/vc/dv/rootfs/opt/conda/envs/bio/bin/samtools; REF=~/certa-bench/idx/ref.fa
IDX=~/certa-bench/idx/grch38.s1v3.cidx; S=~/vc/subset_HG002; C=~/certa-prq/build/certa
step() { echo "== $(date +%T) $*"; }
step "lost names"
$ST view -T $REF -F 0x904 -q 20 ~/vc/HG002/minibwa.cram | cut -f1-4 | awk -F'\t' '{print $1 "\t" int($2 / 64) % 4 "\t" $3 "\t" $4}' > mb_pos.tsv
zstd -dc ~/vc/diag/t2tfeat_HG002.tsv.zst | awk -F'\t' 'NR == FNR { c[$1 "\t" $2] = $3; p[$1 "\t" $2] = $4; next }
  { k = $1 "\t" $2; if (!(k in c)) next; d = $4 - p[k]; if ($3 == "*" || $3 != c[k] || d > 1000 || d < -1000) print $1 }' \
  mb_pos.tsv /dev/stdin | sort -u > names.txt
wc -l < names.txt
step "minibwa's records"
$ST view -T $REF -F 0x900 -N names.txt ~/vc/HG002/minibwa.cram > mb.sam
step "CERTA's records"
$C map $IDX $S/R1.fq $S/R2.fq --intervals --iv-pair all --iv-rescue local --feature-tags --witness -t 64 -o /dev/stdout 2> certa.log \
  | awk -F'\t' 'NR == FNR { w[$1] = 1; next } /^@/ || ($1 in w)' names.txt - > certa.sam
grep -c -v "^@" mb.sam certa.sam
step "diagnosis"
python3 ~/vc/scripts/lost_diag.py $ST $REF mb.sam certa.sam > diag.md
cat diag.md
step "done"
