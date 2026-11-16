#!/bin/bash
# specs/164: one server's chain of curve arms. Usage: run_164.sh NAME PORT ARM...
# ARM is cf_off or a bench/analysis/129 table stem (cf_K9_h075, ...). Each arm:
# start a replay server on the served 98 build with exactly that DS_EDGE_CURVES,
# verify edge_curve_rule POSITIVELY from a solve, run korkut, decryptM and CPPA
# benches, stop. Outputs under ~/deltasignal-catalogs/analysis/164/<arm>/.
set -u
NAME=$1; PORT=$2; shift 2
BUILD=$HOME/deltasignal-catalogs/builds/20261113-1311_baf608b_r98fix
A=$HOME/deltasignal-catalogs/analysis
PY=/home/awright/.cache/pypoetry/virtualenvs/logic-network-generator-ndE5cM8Y-py3.14/bin/python
cd ~/gitroot/deltasignal
for ARM in "$@"; do
  OUT=$A/164/$ARM; mkdir -p "$OUT"
  docker stop "$NAME" >/dev/null 2>&1
  if [ "$ARM" = cf_off ]; then EXTRA=(); else EXTRA=(-e "DS_EDGE_CURVES=/app/bench/analysis/129/$ARM.tsv"); fi
  PATHWAY_CATALOG=$BUILD docker compose -f docker-compose.dev.yml run --rm -d --name "$NAME" \
    -p "$PORT:8080" -e DS_SOFT_OBS=off "${EXTRA[@]}" julia-api >/dev/null 2>&1
  curl -s --retry 60 --retry-delay 3 --retry-all-errors -m 5 "http://localhost:$PORT/api/health" >/dev/null \
    || { echo "ABORT $ARM: server never healthy"; exit 1; }
  RULE=$(DELTASIGNAL_BASE=http://localhost:$PORT DS_CATALOG_ROOT=$BUILD $PY - <<PYEOF
import json, sys, urllib.request
sys.path.insert(0, "bench")
import benchmark_vs_mpbiopath as bm
p = bm.parse_pathway_via_ds_api("R-HSA-109606")
u = p["nodes"][0]["uuid"]
rq = urllib.request.Request(bm.DS_BASE + "/api/solve",
     json.dumps({"network_id": p["network_id"], "observations": {u: [0.0, 1.0]}}).encode(),
     {"Content-Type": "application/json"})
print(json.loads(urllib.request.urlopen(rq, timeout=600).read()).get("edge_curve_rule", "?"))
PYEOF
)
  case "$ARM" in
    cf_off) [ "$RULE" = off ] || { echo "ABORT $ARM: edge_curve_rule=$RULE, expected off"; exit 1; } ;;
    *_h05) echo "$RULE" | grep -q "h=0.5" || { echo "ABORT $ARM: rule=$RULE"; exit 1; } ;;
    *_h075) echo "$RULE" | grep -q "h=0.75" || { echo "ABORT $ARM: rule=$RULE"; exit 1; } ;;
    *_h1) echo "$RULE" | grep -q "h=1" || { echo "ABORT $ARM: rule=$RULE"; exit 1; } ;;
  esac
  echo "$RULE" > "$OUT/edge_curve_rule.txt"
  echo "=== $ARM on :$PORT (verified: ${RULE:0:60}) ==="
  export DELTASIGNAL_BASE=http://localhost:$PORT DS_CATALOG_ROOT=$BUILD
  $PY bench/korkut_bench.py "$OUT/korkut" > "$OUT/korkut.log" 2>&1 || { echo "ABORT $ARM korkut"; exit 1; }
  $PY bench/decryptm_bench.py $A/114/cases_frozen92_kinobeads.tsv \
      $HOME/deltasignal-catalogs/external/lincs_p100/analysis/phospho_ewas_r98.csv "$OUT/decryptm" \
      > "$OUT/decryptm.log" 2>&1 || { echo "ABORT $ARM decryptm"; exit 1; }
  $PY bench/cppa_bench.py $A/121/drug_targets.tsv $A/121/acc_gene_r97.tsv $A/121/readouts_cppa.tsv \
      "$OUT/cppa" > "$OUT/cppa.log" 2>&1 || { echo "ABORT $ARM cppa"; exit 1; }
  echo "$ARM done"
done
docker stop "$NAME" >/dev/null 2>&1
echo "chain $NAME complete"
