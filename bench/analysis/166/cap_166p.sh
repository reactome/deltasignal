#!/bin/bash
# the 21600 s cap for lit166p shuffle arms (the scheduler's own cap pattern predates the lit166p name)
LOG=$HOME/deltasignal-catalogs/analysis/166/run2/sched.log
while ps -eo args= | grep -q "[s]ched_166.sh"; do
  ps -eo pid=,etimes=,args= | grep "[r]un_arm.sh lit166p --catalog" | while read -r pid et rest; do
    if [ "$et" -gt 21600 ]; then
      echo "$(echo "$rest" | sed 's/.*--catalog \([^ ]*\).*/\1/') stopped at $(date -Is) after ${et}s (cap)" >> "$LOG"; kill "$pid"
    fi
  done
  sleep 120
done
