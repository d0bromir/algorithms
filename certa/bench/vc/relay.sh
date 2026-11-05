#!/bin/bash
# relay.sh GENOME (runs here): waits for galaxy's five filtered BAMs of
# GENOME, converts them to CRAM there (cores 120-127), copies them to a2
# through this machine (the hosts cannot reach each other), checks md5s, and
# starts DeepVariant and hap.py on a2.
set -euo pipefail
G=$1
names="ivcal pair pairt2t minibwa bwa-mem2"
while true; do
  n_ok=$(ssh galaxy "ls ~/vc/$G/*.ok 2>/dev/null | wc -l")
  fin=$(ssh galaxy "grep -c '$G finished' ~/vc/$G.map.log 2>/dev/null || true")
  [ "$n_ok" -ge 5 ] && break
  [ "${fin:-0}" -ge 1 ] && { echo "$G mapping finished with $n_ok of 5 BAMs"; break; }
  sleep 120
done
echo "== $(date +%T) $G: CRAM on galaxy"
ssh galaxy "cd ~/vc/$G && for n in $names; do [ -e \$n.ok ] || continue; [ -s \$n.cram.crai ] || { taskset -c 120-127 ~/micromamba/envs/certa-bench/bin/samtools view -@ 8 -C -T ~/certa-bench/idx/ref.fa -o \$n.cram \$n.bam && taskset -c 120-127 ~/micromamba/envs/certa-bench/bin/samtools index \$n.cram; }; done; md5sum *.cram *.crai > cram.md5; ls -la *.bam *.cram"
ssh a2 "mkdir -p ~/vc/$G"
scp -3 -q galaxy:vc/$G/cram.md5 a2:vc/$G/cram.md5
for n in $names; do
  for f in $n.cram $n.cram.crai $n.time; do
    ssh galaxy "[ -e ~/vc/$G/$f ]" || continue
    echo "$(date +%T) copy $f"
    scp -3 -q galaxy:vc/$G/$f a2:vc/$G/$f
  done
done
ssh a2 "cd ~/vc/$G && md5sum -c cram.md5"
echo "== $(date +%T) $G: DeepVariant on a2"
ssh a2 "cd ~/vc && nohup scripts/call.sh $G $names > $G.call.log 2>&1 < /dev/null &" &
sleep 20
ssh a2 "tail -n 3 ~/vc/$G.call.log"
echo "== $(date +%T) relay of $G finished"
