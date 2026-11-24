#!/bin/bash
# specs/170: run the registered arms through run_arm.sh (run name lit170), at most SLOTS at once (default 3).
# Usage: sched_168.sh ARM...   (build dir names under ~/deltasignal-catalogs/builds; *-shufNN arms run the gate axes only)
cd ~/gitroot/deltasignal || exit 1
LOG=$HOME/deltasignal-catalogs/analysis/170/sched.log
mkdir -p "$(dirname "$LOG")"
CAP=21600
PORTS=(8111 8112 8113 8114)
running() { ps -eo args= | grep -c "[r]un_arm.sh lit170"; }
cap() {
  ps -eo pid=,etimes=,args= | grep "[r]un_arm.sh lit170 --catalog [^ ]*-shuf" | while read -r pid et rest; do
    [ "$et" -gt "$CAP" ] && { echo "$(echo "$rest" | sed 's/.*--catalog \([^ ]*\).*/\1/') stopped at $(date -Is) after ${et}s" >> "$LOG"; kill "$pid"; }
  done
}
for a in "$@"; do
  while [ "$(running)" -ge "${SLOTS:-3}" ]; do cap; sleep 60; done
  for p in "${PORTS[@]}"; do ps -eo args= | grep -q "[r]un_arm.sh lit170.*--port $p " || break; done
  echo "start $a on $p at $(date -Is)" >> "$LOG"
  case "$a" in *-shuf*) AX=(--axis experimental --axis korkut) ;; *) AX=(--axis curator --axis experimental --axis korkut) ;; esac
  ( scripts/run_arm.sh lit170 --catalog "$a" "${AX[@]}" --port "$p" --bench DS_SOLVE_TIMEOUT=3600 --bench DS_PIN_CURATED_ROOT=$HOME/deltasignal-catalogs/builds/20261120-lit-v3-without \
      > "$HOME/deltasignal-catalogs/analysis/170/arm-$a.log" 2>&1
    echo "end $a exit $? at $(date -Is)" >> "$LOG" ) &
  sleep 30
done
while [ "$(running)" -gt 0 ]; do cap; sleep 60; done
echo "all done $(date -Is)" >> "$LOG"
