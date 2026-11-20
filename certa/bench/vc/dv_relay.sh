#!/bin/bash
# dv_relay.sh TAG (here): galaxy's ~/vc/G/S_TAG.bam (+ .bai, .time) to a2,
# md5-checked, then call.sh G S_TAG on a2 for HG002 and HG005 in turn.
set -euo pipefail
TAG=$1; n=S_$TAG
for G in HG002 HG005; do
  until ssh galaxy "test -s ~/vc/$G/$n.bam.md5"; do sleep 60; done
  for f in $n.bam $n.bam.bai $n.time; do scp -3 -q galaxy:vc/$G/$f a2:vc/$G/$f; done
  m1=$(ssh galaxy "cat ~/vc/$G/$n.bam.md5"); m2=$(ssh a2 "md5sum < ~/vc/$G/$n.bam")
  [ "$m1" = "$m2" ] || { echo "$G $n md5 differs"; exit 1; }
  echo "$(date +%T) $G $n on a2; calling"
  ssh a2 "cd ~/vc && scripts/call.sh $G $n > ~/vc/$G.$n.call.log 2>&1" || echo "$G calls failed"
  ssh a2 "tail -2 ~/vc/$G.$n.call.log"
done
