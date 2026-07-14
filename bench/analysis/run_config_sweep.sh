#!/usr/bin/env bash
# Run the 9-pathway experimental benchmark under several configurations
# and dump per-case TSVs side-by-side. Use compare_configs.py to diff.
#
# Each config restarts the Julia API server (env vars are read at startup),
# waits for /api/health, then runs benchmark_selective_joint.py.
#
# Outputs: /tmp/ds_<label>.tsv
set -euo pipefail
cd "$(dirname "$0")/../.."

api_health() {
  until curl -s http://localhost:8080/api/health 2>/dev/null | grep -q ok; do
    sleep 3
  done
}

run_one() {
  local label="$1"; shift
  echo "============================================================"
  echo "Config: $label"
  for kv in "$@"; do echo "  $kv"; done
  echo "============================================================"
  env "$@" docker compose -f docker-compose.dev.yml up -d julia-api --force-recreate
  api_health
  env DS_DUMP_CASES="/tmp/ds_${label}.tsv" \
    python3 bench/benchmark_selective_joint.py 2>&1 \
    | tee "/tmp/ds_${label}.log" | tail -16
  echo
}

# Baseline (current published best)
run_one or_max \
  DS_AND_MODE=hill_log DS_INHIBITION_MODE=devspec DS_INHIBITOR_BETA=2 \
  DS_OR_MODE=max DS_DEPLETION_H_MAX=10.0 DS_HILL_LOG_ZMAX=10.0

# Mean OR — captures KOs propagating but over-predicts paralog compensation
run_one or_mean \
  DS_AND_MODE=hill_log DS_INHIBITION_MODE=devspec DS_INHIBITOR_BETA=2 \
  DS_OR_MODE=mean DS_DEPLETION_H_MAX=10.0 DS_HILL_LOG_ZMAX=10.0

# Median OR — "majority wins" alternative
run_one or_median \
  DS_AND_MODE=hill_log DS_INHIBITION_MODE=devspec DS_INHIBITOR_BETA=2 \
  DS_OR_MODE=median DS_DEPLETION_H_MAX=10.0 DS_HILL_LOG_ZMAX=10.0

# Lower z_max → continuous outputs (publication figure: GSEA-style ranking)
run_one zmax_4p6 \
  DS_AND_MODE=hill_log DS_INHIBITION_MODE=devspec DS_INHIBITOR_BETA=2 \
  DS_OR_MODE=max DS_DEPLETION_H_MAX=10.0 DS_HILL_LOG_ZMAX=4.605

echo "All configs done. Compare with:"
echo "  python bench/analysis/compare_configs.py /tmp/ds_or_max.tsv /tmp/ds_or_mean.tsv \\"
echo "    --label-a or_max --label-b or_mean"
