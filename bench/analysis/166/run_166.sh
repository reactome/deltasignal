#!/bin/bash
# specs/166: run every arm through scripts/run_arm.sh, as worker W of 4 (port 8100+W).
# Usage: bench/analysis/166/run_166.sh W
set -u
W=$1
# primary arms first (without, with, shuffles), descriptive slices last
mapfile -t ARMS < <(ls -d "$HOME"/deltasignal-catalogs/builds/20261118-lit-h12b-* | xargs -n1 basename |
  awk '{k = /-without$/ ? 0 : /-with$/ ? 1 : /-shuf/ ? 2 : 3; print k "\t" $0}' | sort | cut -f2)
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
