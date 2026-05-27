#!/usr/bin/env bash
# Re-tune with inhibition_mode grid:
#   - spec: multiplicative Hill with β=0 default → no effect (current behavior)
#   - inversion: invert source state, treat as AND activator (naive-baseline-style)
#
# Locks AND_MODE=signed (winner from Phase 1) and thresholds (0.5, 2.0).
set -eo pipefail

MAX_EDGES=80000
DOWN=0.5
UP=2.0
REPORT_DIR=/tmp/ds_phase1_inhib
mkdir -p "$REPORT_DIR"
SUMMARY="$REPORT_DIR/summary.tsv"
printf "and_mode\tinhibition_mode\tdown\tup\te2e_correct\te2e_total\te2e_acc\tvalid_correct\tvalid_total\tvalid_acc\n" > "$SUMMARY"

restart_api() {
    local and=$1 inhib=$2
    echo "  -- DS_AND_MODE=$and DS_INHIBITION_MODE=$inhib --"
    DS_AND_MODE=$and DS_INHIBITION_MODE=$inhib \
        docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api 2>&1 | tail -2
    until curl -s -o /dev/null -m 2 http://127.0.0.1:8080/api/health; do sleep 3; done
    sleep 2
}

run_one() {
    local and=$1 inhib=$2
    local label="${and}_${inhib}"
    local report="$REPORT_DIR/${label}.tsv"

    echo
    echo "=== AND=${and}  INHIB=${inhib}  thresholds=(${DOWN}, ${UP}) ==="
    DS_DOWN_CUTOFF=$DOWN DS_UP_CUTOFF=$UP \
        python3 bench/benchmark_vs_mpbiopath.py \
            --ground-truth experimental \
            --max-edges "$MAX_EDGES" \
            --report "$report" \
        2>&1 | grep -E "^  [A-Za-z]|DeltaSignal|Valid|Confusion|pred=" || true

    python3 - "$report" "$and" "$inhib" "$DOWN" "$UP" "$SUMMARY" <<'PY'
import csv, sys
report, and_mode, inhib_mode, down, up, summary = sys.argv[1:7]
e2e_c=e2e_t=v_c=v_t=0
with open(report) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row.get('status') != 'ok': continue
        e2e_c += int(row['correct']); e2e_t += int(row['total'])
        v_c += int(row['valid_correct']); v_t += int(row['valid_total'])
acc = e2e_c/e2e_t if e2e_t else 0
vacc = v_c/v_t if v_t else 0
with open(summary,'a') as f:
    f.write(f"{and_mode}\t{inhib_mode}\t{down}\t{up}\t{e2e_c}\t{e2e_t}\t{acc:.4f}\t{v_c}\t{v_t}\t{vacc:.4f}\n")
PY
}

cd /home/awright/gitroot/deltasignal

# Six configs (3 AND × 2 inhibition). Restart between every distinct env combo.
for and in geomean min signed; do
    for inhib in spec inversion; do
        restart_api "$and" "$inhib"
        run_one "$and" "$inhib"
    done
done

echo
echo "===================== SUMMARY ====================="
column -t -s $'\t' "$SUMMARY"
echo "Reports in $REPORT_DIR"
