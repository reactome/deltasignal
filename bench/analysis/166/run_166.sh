#!/bin/bash
# specs/166: run every arm through scripts/run_arm.sh, as worker W of 4 (port 8100+W).
# Superseded on 2026-11-19 by sched_166.sh (amendment 1); kept as the record of the first runs.
# Usage: bench/analysis/166/run_166.sh W
set -u
W=$1
# primary arms first (without, with, shuffles), descriptive slices last
mapfile -t ARMS < <(ls -d "$HOME"/deltasignal-catalogs/builds/20261118-lit-h12b-* | xargs -n1 basename |
  awk '{k = /-without$/ ? 0 : /-with$/ ? 1 : /-shuf/ ? 2 : 3; print k "\t" $0}' | sort | cut -f2)
i=0
for a in "${ARMS[@]}"; do
  if [ $((i % 4 + 1)) -eq "$W" ]; then
    sha=$(git rev-parse --short HEAD)
    if python3 -c 'import json,sys; sys.exit(json.load(open(sys.argv[1]))["exit"])' \
         "$HOME/deltasignal-catalogs/builds/$a/results/$sha/lit166/ARM.json" 2>/dev/null; then
      echo "=== $a already complete at $sha"; i=$((i + 1)); continue
    fi
    echo "=== $a ($(date -Is))"
    # a shuffled network can hold one slow solve (> the default 600 s, shuf02 Mitotic G1)
    scripts/run_arm.sh lit166 --catalog "$a" --axis curator --axis experimental --axis korkut --port $((8100 + W)) \
      --bench DS_SOLVE_TIMEOUT=3600 \
      || echo "FAILED $a"
  fi
  i=$((i + 1))
done
echo "worker $W complete"
