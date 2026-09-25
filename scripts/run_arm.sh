#!/bin/bash
# Run one benchmark ARM, reproducibly, and leave a record of exactly what ran.
#
#   scripts/run_arm.sh NAME [--server K=V]... [--bench K=V]... [--port N] [--limit N]
#
#   --limit N      benchmark only the first N pathways (a smoke test, not an arm)
#
#   --server K=V   a DS_* override for the SOLVER (set in the API container)
#   --bench  K=V   an override for the BENCHMARK script (e.g. DS_PIN_SCOPE=entry)
#
# Why this exists: arms used to run from one-off scratch scripts, so the only
# record of which overrides were live, which container served them and which
# commit they ran was the conversation that launched them. specs/023 found the
# benchmark had silently been pinning every complex containing a gene for two
# months; nobody could have seen it from the results alone.
#
# What it guarantees:
#   * the code is a DETACHED worktree at HEAD (src/ and bench/ must be
#     committed), so editing the repo mid-run cannot change the arm;
#   * the API container runs that worktree, with the dev compose file's
#     environment plus --server overrides, on the current catalog BUILD
#     (resolved once, mounted by its real path, not the moving `current` link);
#   * it aborts unless the container reports that build and each --server
#     override is visible in the container's own environment;
#   * results go to builds/<id>/results/<sha>/<NAME>/ with ARM.json recording
#     sha, build, overrides, the protocol line the benchmark printed, and exits.
#
# The control for an arm is the production scoring of the same build at a
# solver-identical commit, or an arm run with no overrides.
set -u -o pipefail

die() { echo "run_arm: $*" >&2; exit 1; }

[ $# -ge 1 ] || die "usage: scripts/run_arm.sh NAME [--server K=V]... [--bench K=V]... [--port N]"
NAME=$1; shift
[[ "$NAME" =~ ^[A-Za-z0-9._-]+$ ]] || die "NAME must be [A-Za-z0-9._-]+"
SERVER=(); BENCH=(); PORT=8090; LIMIT=()
while [ $# -gt 0 ]; do
  case "$1" in
    --server) [[ "${2:-}" == DS_*=* ]] || die "--server needs DS_NAME=value"; SERVER+=("$2"); shift 2 ;;
    --bench)  [[ "${2:-}" == *=* ]] || die "--bench needs NAME=value"; BENCH+=("$2"); shift 2 ;;
    --port)   PORT=$2; shift 2 ;;
    --limit)  LIMIT=(--limit "$2"); shift 2 ;;
    *) die "unknown argument $1" ;;
  esac
done

REPO=$(git rev-parse --show-toplevel) || die "not in a git repo"
cd "$REPO"
dirty=$(git status --porcelain -- src bench | wc -l)
[ "$dirty" = "0" ] || die "src/ or bench/ has uncommitted changes; an arm must run a commit"
SHA=$(git rev-parse --short HEAD)
BUILD=$(readlink -f "$HOME/deltasignal-catalogs/current") || die "no current catalog"
[ -f "$BUILD/BUILD.json" ] || die "$BUILD has no BUILD.json"
BID=$(basename "$BUILD")
OUT="$BUILD/results/$SHA/$NAME"
[ -e "$OUT" ] && die "$OUT already exists; pick a new NAME or remove it deliberately"
PY="$(cd "$HOME/gitroot/logic-network-generator" && poetry env info -p)/bin/python"
[ -x "$PY" ] || die "generator venv python not found (the benchmark needs py2neo)"

WT=$(mktemp -d "$HOME/gitroot/.arm-$NAME-XXXX")
git worktree add --detach "$WT" "$SHA" -q || die "worktree failed"
# An override the pinned code never reads is set, verified, recorded -- and
# inert. That silently measures the control. Refuse it.
for kv in "${SERVER[@]}"; do
  grep -rqF "\"${kv%%=*}\"" "$WT/src" || { git worktree remove --force "$WT"; die "${kv%%=*} is not read anywhere in src/ at $SHA"; }
done
for kv in "${BENCH[@]}"; do
  grep -rqF "\"${kv%%=*}\"" "$WT/bench" || { git worktree remove --force "$WT"; die "${kv%%=*} is not read anywhere in bench/ at $SHA"; }
done
CNAME="arm-$NAME-$$"
cleanup() { docker rm -f "$CNAME" >/dev/null 2>&1; git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1; }
trap cleanup EXIT

# The solver environment is the dev compose file's, so an arm differs from the
# dev API only by what is named on the command line.
ENVF="$WT/.arm.env"
docker compose -f "$WT/docker-compose.dev.yml" config --format json 2>/dev/null |
  python3 -c 'import json,sys; e=json.load(sys.stdin)["services"]["julia-api"].get("environment") or {}
[print(f"{k}={v}") for k,v in sorted(e.items()) if v is not None]' > "$ENVF" || die "could not read compose env"
EXTRA=(); for kv in "${SERVER[@]}"; do EXTRA+=(-e "$kv"); done

docker run -d --name "$CNAME" -p "$PORT:8080" --env-file "$ENVF" "${EXTRA[@]}" \
  -v "$WT/Project.toml:/app/Project.toml:ro" -v "$WT/src:/app/src:ro" -v "$WT/bench:/app/bench:ro" \
  -v "$WT/cli:/app/cli:ro" -v "$WT/examples:/app/examples:ro" -v "$WT/test:/app/test:ro" \
  -v "$BUILD:/app/pathway_catalog:ro" -v deltasignal_julia-cache:/root/.julia \
  -w /app deltasignal-julia-api julia --project=/app /app/src/api/server.jl --port=8080 --host=0.0.0.0 \
  >/dev/null || die "container failed to start"

for i in $(seq 1 120); do curl -sf "localhost:$PORT/api/health" >/dev/null && break; sleep 5; done
served=$(curl -s "localhost:$PORT/api/health" | python3 -c 'import json,sys; print(json.load(sys.stdin)["catalog"]["build_id"])') \
  || die "API on $PORT never became healthy"
[ "$served" = "$BID" ] || die "API serves $served, expected $BID"
for kv in "${SERVER[@]}"; do
  docker exec "$CNAME" env | grep -qx "$kv" || die "override $kv is NOT in the container environment"
done

mkdir -p "$OUT"
STARTED=$(date -Is); status=0
for gt in curator experimental; do
  env "${BENCH[@]}" DELTASIGNAL_BASE="http://localhost:$PORT" DS_CATALOG_ROOT="$BUILD" \
    "$PY" "$WT/bench/benchmark_vs_mpbiopath.py" --max-edges 40000 "${LIMIT[@]}" --ground-truth "$gt" \
    --report "$OUT/${gt}_report.tsv" --dump-cases "$OUT/${gt}_cases.tsv" > "$OUT/$gt.log" 2>&1 \
    || { echo "run_arm: $gt benchmark failed; see $OUT/$gt.log" >&2; status=1; break; }
done

ARM_LIMIT="${LIMIT[1]:-}" python3 - "$OUT" "$NAME" "$SHA" "$BID" "$STARTED" "$status" "${SERVER[*]:-}" "${BENCH[*]:-}" <<'PY'
import json, os, re, sys, datetime
out, name, sha, bid, started, status, server, bench = sys.argv[1:]
rec = {"name": name, "commit": sha, "build": bid, "started": started,
       "finished": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
       "exit": int(status), "limit": os.environ.get("ARM_LIMIT") or None,
       "server_overrides": server.split() if server else [],
       "bench_overrides": bench.split() if bench else []}
for gt in ("curator", "experimental"):
    log = os.path.join(out, f"{gt}.log")
    if os.path.exists(log):
        text = open(log).read()
        m = re.search(r"^Protocol: .*$", text, re.M); rec[f"{gt}_protocol"] = m.group(0) if m else None
        m = re.search(r"^Pinned: .*$", text, re.M); rec[f"{gt}_pinned"] = m.group(0) if m else None
json.dump(rec, open(os.path.join(out, "ARM.json"), "w"), indent=2)
print(json.dumps(rec, indent=2))
PY
[ "$status" = "0" ] || exit 1
echo "run_arm: done -> $OUT"
