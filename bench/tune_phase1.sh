#!/usr/bin/env bash
# Phase 1 tuning: grid search across {AND aggregator mode} × {discretization
# thresholds} on the 7 experimental-ground-truth pathways we can serve via
# the HTTP API. The other 3 of the 10 experimental pathways are excluded
# because RAF_MAP_kinase_cascade (1.6M edges) and Mitotic_Prophase (272k)
# exceed the HTTP-benchmark's practical size limit, and the regenerated
# HDR_through_HRR_SSA network is missing from the catalog.
set -eo pipefail

# Path-edges threshold that keeps PIP3, ERBB2, Cell_Cycle_Checkpoints,
# Transcriptional_Regulation_by_TP53, Signaling_by_WNT in and the giants out.
MAX_EDGES=80000

DOWN_TIGHT="0.5";  UP_TIGHT="2.0"     # spec-inspired: half / double baseline
DOWN_LOOSE="0.85"; UP_LOOSE="1.15"    # ±15% sensitivity (slight change counts)
DOWN_WIDE="0.2";   UP_WIDE="5.0"      # only big changes count

REPORT_DIR=/tmp/ds_phase1_tuning
mkdir -p "$REPORT_DIR"
SUMMARY="$REPORT_DIR/summary.tsv"
printf "and_mode\tdown\tup\tend_to_end_correct\tend_to_end_total\tend_to_end_acc\tvalid_correct\tvalid_total\tvalid_acc\n" > "$SUMMARY"

restart_api_with_mode() {
    local mode=$1
    echo "  -- restarting julia-api with DS_AND_MODE=$mode --"
    DS_AND_MODE=$mode docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api 2>&1 | tail -2
    until curl -s -o /dev/null -m 2 http://127.0.0.1:8080/api/health; do
        sleep 3
    done
    sleep 2  # safety margin after health-check passes
}

run_one() {
    local mode=$1 down=$2 up=$3
    local label="${mode}_d${down}_u${up}"
    local report="$REPORT_DIR/${label}.tsv"

    echo
    echo "=== AND=${mode}  DOWN<${down}  UP>=${up} ==="
    DS_DOWN_CUTOFF=$down DS_UP_CUTOFF=$up \
        python3 bench/benchmark_vs_mpbiopath.py \
            --ground-truth experimental \
            --max-edges "$MAX_EDGES" \
            --report "$report" \
        2>&1 | grep -E "^  [A-Z]|DeltaSignal|Valid|Confusion|pred=" || true

    # Pull headline numbers from the report
    python3 - "$report" "$mode" "$down" "$up" "$SUMMARY" <<'PY'
import csv, sys
report, mode, down, up, summary = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
e2e_c=e2e_t=v_c=v_t=0
with open(report) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row.get('status') != 'ok': continue
        e2e_c += int(row['correct']); e2e_t += int(row['total'])
        v_c += int(row['valid_correct']); v_t += int(row['valid_total'])
acc  = e2e_c/e2e_t if e2e_t else 0
vacc = v_c/v_t if v_t else 0
with open(summary,'a') as f:
    f.write(f"{mode}\t{down}\t{up}\t{e2e_c}\t{e2e_t}\t{acc:.4f}\t{v_c}\t{v_t}\t{vacc:.4f}\n")
PY
}

cd /home/awright/gitroot/deltasignal

for mode in geomean min signed; do
    restart_api_with_mode "$mode"
    run_one "$mode" "$DOWN_TIGHT"  "$UP_TIGHT"
    run_one "$mode" "$DOWN_LOOSE"  "$UP_LOOSE"
    run_one "$mode" "$DOWN_WIDE"   "$UP_WIDE"
done

echo
echo "==================== GRID SUMMARY ===================="
column -t "$SUMMARY"
echo
echo "Reports per config in $REPORT_DIR"
