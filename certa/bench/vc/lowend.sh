#!/bin/bash
# a2: the low end of kBwaScale for DeepVariant's MAPQ floor (5): S_wpair (--mapq pair, refit, wide band) with
# Phred 3 (MAPQ 10) below the floor (L1), and Phred 3 and 4 (MAPQ 10, 16) below it (L2); then call.sh.
set -euo pipefail
cd ~/vc
for G in HG002 HG005; do
  scripts/bwscale.sh $G S_wpair S_wpairL1 0,1,2,3,4,5,6,7,8,9,3,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60
  scripts/bwscale.sh $G S_wpair S_wpairL2 0,1,2,3,4,5,6,7,8,9,3,11,12,13,14,15,4,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60
done
scripts/call.sh HG002 S_wpairL1 S_wpairL2 > HG002.lowend.call.log 2>&1
scripts/call.sh HG005 S_wpairL1 S_wpairL2 > HG005.lowend.call.log 2>&1
echo "== $(date +%T) lowend done"
