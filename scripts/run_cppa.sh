#!/bin/bash
# CPPA single-drug predictions (specs/121; values unread) on the top114 build from a detached worktree at HEAD, as
# run_arm.sh pins an arm: N servers with the dev compose environment, one bench shard each. Behind specs/134 R4/R7/R9.
#   scripts/run_cppa.sh NAME DS_TOTAL_READOUT [N]       output: builds/<top114>/results/<sha>/NAME
# CPPA_PATHWAYS=a,b,... limits the run to those pathways (recorded in RUN.txt). PORT_BASE (default 8094) sets ports.
set -euo pipefail
[ $# -ge 2 ] || { echo "usage: $0 NAME DS_TOTAL_READOUT [N]"; exit 2; }
NAME=$1; MODE=$2; N=${3:-4}
REPO=$HOME/gitroot/deltasignal; cd "$REPO"
[ -z "$(git status --porcelain src bench)" ] || { echo "dirty"; exit 1; }
SHA=$(git rev-parse --short HEAD); BID=20261027-1346_8b6527e_top114; BUILD=$HOME/deltasignal-catalogs/builds/$BID
OUT=$BUILD/results/$SHA/$NAME; [ -e "$OUT" ] && { echo "$OUT exists"; exit 1; }
PY="$(cd $HOME/gitroot/logic-network-generator && poetry env info -p)/bin/python"
export NEO4J_PASSWORD="$(docker exec reactome-neo4j sh -c 'echo ${NEO4J_AUTH#neo4j/}')"
WT=$(mktemp -d "$HOME/gitroot/.arm-$NAME-XXXX"); git worktree add --detach "$WT" "$SHA" -q
ENVF="$WT/.arm.env"
env -i HOME="$HOME" PATH="$PATH" docker compose -f "$WT/docker-compose.dev.yml" config --format json 2>/dev/null |
  python3 -c 'import json,sys; e=json.load(sys.stdin)["services"]["julia-api"].get("environment") or {}
[print(f"{k}={v}") for k,v in sorted(e.items()) if v is not None]' > "$ENVF"
NAMES=()
cleanup() { for c in "${NAMES[@]}"; do docker rm -f "$c" >/dev/null 2>&1; done; git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1; }
trap cleanup EXIT
mkdir -p "$OUT"; A=$HOME/deltasignal-catalogs/analysis/121
sha256sum $A/drug_targets.tsv $A/acc_gene_r97.tsv $A/readouts_cppa.tsv $A/split.tsv > "$OUT/INPUTS.sha256"
echo "commit $SHA build $BID DS_TOTAL_READOUT=$MODE DS_SOLVE_TIMEOUT=3600 shards $N pathways ${CPPA_PATHWAYS:-all}" > "$OUT/RUN.txt"
pids=()
for i in $(seq 0 $((N-1))); do
  PORT=$((${PORT_BASE:-8094}+i)); C="arm-$NAME-$i-$$"; NAMES+=("$C")
  docker run -d --name "$C" -p "$PORT:8080" --env-file "$ENVF" \
    -v "$WT/Project.toml:/app/Project.toml:ro" -v "$WT/src:/app/src:ro" -v "$WT/bench:/app/bench:ro" \
    -v "$WT/cli:/app/cli:ro" -v "$WT/examples:/app/examples:ro" -v "$WT/test:/app/test:ro" \
    -v "$BUILD:/app/pathway_catalog:ro" -v deltasignal_julia-cache:/root/.julia \
    -w /app deltasignal-julia-api julia --project=/app /app/src/api/server.jl --port=8080 --host=0.0.0.0 >/dev/null
  cp "$ENVF" "$OUT/server.env.${i}of$N"
done
for i in $(seq 0 $((N-1))); do
  PORT=$((${PORT_BASE:-8094}+i))
  for t in $(seq 1 120); do curl -sf "localhost:$PORT/api/health" >/dev/null && break; sleep 5; done
  [ "$(curl -s localhost:$PORT/api/health | python3 -c 'import json,sys; print(json.load(sys.stdin)["catalog"]["build_id"])')" = "$BID" ] || { echo "port $PORT wrong build"; exit 1; }
done
for i in $(seq 0 $((N-1))); do
  PORT=$((${PORT_BASE:-8094}+i))
  env DS_TOTAL_READOUT=$MODE DS_SOLVE_TIMEOUT=3600 CPPA_SHARD=$i/$N CPPA_PATHWAYS="${CPPA_PATHWAYS:-}" DELTASIGNAL_BASE="http://localhost:$PORT" DS_CATALOG_ROOT="$BUILD" \
    "$PY" "$WT/bench/cppa_bench.py" $A/drug_targets.tsv $A/acc_gene_r97.tsv $A/readouts_cppa.tsv "$OUT" > "$OUT/cppa_predict.${i}of$N.log" 2>&1 &
  pids+=($!)
done
st=0; for p in "${pids[@]}"; do wait $p || st=1; done
echo "exit $st" >> "$OUT/RUN.txt"; echo "$OUT"; exit $st
