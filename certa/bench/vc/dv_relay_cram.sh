#!/bin/bash
# dv_relay_cram.sh TAG... (here, the machine that reaches both hosts):
# galaxy's ~/vc/G/S_TAG.cram (+ .crai, .time) for every TAG to a2,
# md5-checked; then call.sh on a2, HG002's runs and then HG005's, in the
# background there (~/vc/G.final.call.log).
set -euo pipefail
names=""
for TAG in "$@"; do
  n=S_$TAG; names="$names $n"
  for G in HG002 HG005; do
    until ssh galaxy "test -s ~/vc/$G/$n.cram.md5"; do sleep 60; done
    for f in $n.cram $n.cram.crai $n.time; do scp -3 -q galaxy:vc/$G/$f a2:vc/$G/$f; done
    m1=$(ssh galaxy "cat ~/vc/$G/$n.cram.md5"); m2=$(ssh a2 "md5sum < ~/vc/$G/$n.cram")
    [ "$m1" = "$m2" ] || { echo "$G $n md5 differs"; exit 1; }
    echo "$(date +%T) $G $n on a2"
  done
done
ssh a2 "cd ~/vc && setsid nohup sh -c 'scripts/call.sh HG002$names > HG002.final.call.log 2>&1; scripts/call.sh HG005$names > HG005.final.call.log 2>&1' > /dev/null 2>&1 < /dev/null &"
echo "$(date +%T) calling$names on a2"
