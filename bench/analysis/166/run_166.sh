#!/bin/bash
# specs/166: run every arm through scripts/run_arm.sh, as worker W of 4 (port 8100+W).
# Usage: bench/analysis/166/run_166.sh W
set -u
W=$1
mapfile -t ARMS < <(ls -d "$HOME"/deltasignal-catalogs/builds/20261118-lit-h12b-* | xargs -n1 basename | sort)
i=0
for a in "${ARMS[@]}"; do
  if [ $((i % 4 + 1)) -eq "$W" ]; then
    echo "=== $a ($(date -Is))"
    scripts/run_arm.sh lit166 --catalog "$a" --axis curator --axis experimental --axis korkut --port $((8100 + W)) \
      || echo "FAILED $a"
  fi
  i=$((i + 1))
done
echo "worker $W complete"
