#!/bin/bash
# specs/166 amendment 1 (rev 06:05: shuffles on the gate's axes only, run name lit166p): slices first, then the remaining shuffles, then retry shuf01-03; at most 4 arms at once;
# a shuffle arm running past 21600 s is stopped and recorded as not finished.
cd ~/gitroot/deltasignal || exit 1
LOG=$HOME/deltasignal-catalogs/analysis/166/run2/sched.log
CAP=21600
QUEUE=(shuf08 shuf11 shuf09 shuf10 shuf12 shuf13 shuf14 shuf15 shuf16 shuf17 shuf18 shuf19 shuf04 shuf05 shuf06 shuf07 shuf01 shuf02 shuf03)
running() { ps -eo args= | grep -c "[r]un_arm.sh lit166"; }
cap() {
  ps -eo pid=,etimes=,args= | grep "[r]un_arm.sh lit166 --catalog 20261118-lit-h12b-shuf" | while read -r pid et rest; do
    if [ "$et" -gt "$CAP" ]; then
      echo "$(echo "$rest" | sed 's/.*--catalog \([^ ]*\).*/\1/') stopped at $(date -Is) after ${et}s" >> "$LOG"; kill "$pid"
    fi
  done
}
for a in "${QUEUE[@]}"; do
  while [ "$(running)" -ge "${SLOTS:-2}" ]; do cap; sleep 120; done
  for p in 8101 8102; do ps -eo args= | grep -q "[r]un_arm.sh lit166.*--port $p " || break; done
  echo "start $a on $p at $(date -Is)" >> "$LOG"
  case "$a" in shuf*) NAME=lit166p; AX=(--axis experimental --axis korkut) ;;
                 *) NAME=lit166; AX=(--axis curator --axis experimental --axis korkut) ;; esac
  ( scripts/run_arm.sh "$NAME" --catalog "20261118-lit-h12b-$a" "${AX[@]}" \
      --port "$p" --bench DS_SOLVE_TIMEOUT=3600 > "$HOME/deltasignal-catalogs/analysis/166/run2/arm-$a.log" 2>&1
    echo "end $a exit $? at $(date -Is)" >> "$LOG" ) &
  sleep 30
done
while [ "$(running)" -gt 0 ]; do cap; sleep 120; done
echo "all done $(date -Is)" >> "$LOG"
