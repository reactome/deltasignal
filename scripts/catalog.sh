#!/usr/bin/env bash
# Durable, provenance-carrying DeltaSignal catalogs.
#
# Why this exists. Every catalog behind specs/018-021 was generated into a
# session scratchpad under /tmp, which on this machine is a tmpfs. A reboot on
# 2026-09-24 wiped all of them, and the specs still cite them by name. Nothing
# recorded which catalog the API was serving either: the dev API was quietly
# serving a two-month-old one while every benchmark ran against newer ones.
#
# Layout (on real disk, outside any repo):
#
#   $DS_CATALOG_HOME/                     default ~/deltasignal-catalogs
#     builds/<YYYYmmdd-HHMM>_<lng-sha>/
#       BUILD.json                        what produced it, and whether it completed
#       pathways.tsv                      the exact list the generator was given
#       generate.log
#       R-HSA-*/                          the networks
#       results/<ds-sha>/                 benchmark output for this catalog + solver
#     current -> builds/<id>              the latest build that VERIFIED complete
#
# `current` is only moved after a build completes with every requested pathway
# present, so a failed or partial build can never become what the API serves.
#
# Usage:
#   scripts/catalog.sh build      generate from LNG HEAD, verify, repoint current
#   scripts/catalog.sh status     what current is, and whether it is stale
#   scripts/catalog.sh list       every build with its manifest summary
#   scripts/catalog.sh use <id>   point current at an existing complete build
#   scripts/catalog.sh bench      benchmark current on both axes, filed under the
#                                 build by solver commit; refuses to run unless
#                                 /api/health reports it is serving that build
set -euo pipefail

HOME_DIR="${DS_CATALOG_HOME:-$HOME/deltasignal-catalogs}"
LNG="${LOGIC_NETWORK_GENERATOR:-$HOME/gitroot/logic-network-generator}"
DS_REPO="$(cd "$(dirname "$0")/.." && pwd)"
LIST="$DS_REPO/bench/catalog_pathways.tsv"
PY="${LNG_PYTHON:-$(cd "$LNG" && poetry env info -p 2>/dev/null)/bin/python}"

die() { echo "catalog: $*" >&2; exit 1; }

case "$HOME_DIR" in
  /tmp/*|/dev/shm/*) die "DS_CATALOG_HOME=$HOME_DIR is on a volatile filesystem; that is the failure this script exists to prevent" ;;
esac

manifest_field() { python3 -c "import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2],''))" "$1" "$2"; }

cmd_build() {
  [ -d "$LNG/.git" ] || die "no generator checkout at $LNG"
  [ -x "$PY" ] || die "no generator python at $PY (set LNG_PYTHON)"
  local sha dirty branch id dir
  sha="$(git -C "$LNG" rev-parse --short HEAD)"
  branch="$(git -C "$LNG" rev-parse --abbrev-ref HEAD)"
  dirty="$(git -C "$LNG" status --porcelain src bin | wc -l | tr -d ' ')"
  [ "$dirty" = "0" ] || die "generator has $dirty uncommitted change(s) in src/ or bin/; a build must be reproducible from a commit"
  id="$(date +%Y%m%d-%H%M)_${sha}"
  dir="$HOME_DIR/builds/$id"
  [ ! -e "$dir" ] || die "$dir already exists"
  mkdir -p "$dir"

  # The exact list the generator is given, comments stripped, kept with the build.
  grep -v '^#' "$LIST" > "$dir/pathways.tsv"
  local requested; requested=$(($(wc -l < "$dir/pathways.tsv") - 1))

  echo "catalog: building $id ($requested pathways) into $dir"
  local t0 rc; t0=$(date +%s)
  set +e
  ( cd "$LNG" && env PYTHONHASHSEED=0 LNG_EMIT_ONE_SIDED=1 \
      "$PY" bin/create-pathways.py --pathway-list "$dir/pathways.tsv" --output-dir "$dir" ) \
      > "$dir/generate.log" 2>&1
  rc=$?
  set -e

  python3 - "$dir" "$sha" "$branch" "$rc" "$requested" "$(( $(date +%s) - t0 ))" "$LIST" <<'PY'
import csv, glob, hashlib, json, os, sys, datetime
d, sha, branch, rc, requested, secs, lst = sys.argv[1:]
built = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "R-HSA-*"))
               if os.path.exists(os.path.join(p, "logic_network.csv")))
wanted = [r["id"] for r in csv.DictReader(open(os.path.join(d, "pathways.tsv")), delimiter="\t")]
nodes = edges = 0
for p in built:
    rows = list(csv.DictReader(open(os.path.join(d, p, "logic_network.csv"))))
    edges += len(rows)
    nodes += len({r["source_id"] for r in rows} | {r["target_id"] for r in rows})
missing = sorted(set(wanted) - set(built))
status = "complete" if (rc == "0" and not missing) else "incomplete"
m = {
    "build_id": os.path.basename(d),
    "created": datetime.datetime.now().isoformat(timespec="seconds"),
    "status": status,
    "generator_commit": sha,
    "generator_branch": branch,
    "generator_exit_code": int(rc),
    "env": {"PYTHONHASHSEED": "0", "LNG_EMIT_ONE_SIDED": "1"},
    "pathway_list": os.path.relpath(lst, os.path.dirname(os.path.dirname(lst))),
    "pathway_list_sha256": hashlib.sha256(open(os.path.join(d, "pathways.tsv"), "rb").read()).hexdigest(),
    "pathways_requested": int(requested),
    "pathways_built": len(built),
    "pathways_missing": missing,
    "nodes": nodes,
    "edges": edges,
    "seconds": int(secs),
}
json.dump(m, open(os.path.join(d, "BUILD.json"), "w"), indent=2)
print(json.dumps({k: m[k] for k in ("status", "pathways_built", "pathways_requested", "nodes", "edges")}))
PY

  if [ "$(manifest_field "$dir/BUILD.json" status)" != "complete" ]; then
    echo "catalog: build $id is INCOMPLETE; current left unchanged. See $dir/generate.log" >&2
    exit 1
  fi
  ln -sfn "builds/$id" "$HOME_DIR/current"
  echo "catalog: current -> builds/$id"
}

cmd_status() {
  local cur="$HOME_DIR/current"
  [ -L "$cur" ] || { echo "no current catalog (run: scripts/catalog.sh build)"; exit 1; }
  local id m head
  id="$(basename "$(readlink "$cur")")"
  m="$cur/BUILD.json"
  head="$(git -C "$LNG" rev-parse --short HEAD 2>/dev/null || echo '?')"
  echo "current:   $id"
  echo "status:    $(manifest_field "$m" status)"
  echo "generator: $(manifest_field "$m" generator_commit)  (generator HEAD now: $head)"
  echo "built:     $(manifest_field "$m" created)"
  echo "size:      $(manifest_field "$m" pathways_built) pathways, $(manifest_field "$m" nodes) nodes, $(manifest_field "$m" edges) edges"
  if [ "$head" != "?" ] && [ "$head" != "$(manifest_field "$m" generator_commit)" ]; then
    local changed
    changed="$(git -C "$LNG" diff --stat "$(manifest_field "$m" generator_commit)" HEAD -- src bin | tail -1)"
    echo "STALE:     generator has moved since this build (${changed:-no src/bin change})"
  fi
  local api
  api="$(docker inspect deltasignal-julia-api-1 --format '{{range .Mounts}}{{if eq .Destination "/app/pathway_catalog"}}{{.Source}}{{end}}{{end}}' 2>/dev/null || true)"
  [ -n "$api" ] && echo "dev API:   serving $api"
}

cmd_list() {
  for m in "$HOME_DIR"/builds/*/BUILD.json; do
    [ -e "$m" ] || { echo "no builds yet"; return; }
    printf '%-26s %-10s %3s pathways  %s\n' "$(manifest_field "$m" build_id)" \
      "$(manifest_field "$m" status)" "$(manifest_field "$m" pathways_built)" \
      "$([ "$(readlink "$HOME_DIR/current" 2>/dev/null)" = "builds/$(manifest_field "$m" build_id)" ] && echo '<- current')"
  done
}

cmd_bench() {
  # Benchmark the CURRENT catalog on both ground-truth axes and file the
  # results under that build, keyed by the solver commit that produced them,
  # so every number is traceable to (catalog build, solver commit).
  local cur="$HOME_DIR/current" api="${DELTASIGNAL_BASE:-http://localhost:8080}"
  [ -L "$cur" ] || die "no current catalog (run: scripts/catalog.sh build)"
  local build ds dirty
  build="$(basename "$(readlink "$cur")")"
  ds="$(git -C "$DS_REPO" rev-parse --short HEAD)"
  dirty="$(git -C "$DS_REPO" status --porcelain src | wc -l | tr -d ' ')"
  [ "$dirty" = "0" ] || die "deltasignal src/ has $dirty uncommitted change(s); results must come from a commit"

  # The guard that makes "what are we running" enforced, not just recorded:
  # refuse to score a catalog the server is not actually serving. Scoring one
  # catalog while the server solves against another makes every observation
  # silently match nothing.
  local served
  served="$(curl -sf -m 5 "$api/api/health" | python3 -c \
    "import json,sys; print(json.load(sys.stdin).get('catalog',{}).get('build_id',''))" 2>/dev/null || true)"
  [ -n "$served" ] || die "cannot read the served catalog from $api/api/health (is the API up, and new enough to report it?)"
  [ "$served" = "$build" ] || die "API at $api serves '$served' but current is '$build'; recreate it: docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api"

  local out="$cur/results/$ds"
  mkdir -p "$out"
  echo "catalog: benchmarking build $build with solver $ds -> $out"
  # The benchmark needs py2neo, which lives in the generator's venv, not the
  # system python -- the first run of this command died on exactly that.
  [ -x "$PY" ] || die "no generator python at $PY (set LNG_PYTHON)"
  for gt in curator experimental; do
    env DELTASIGNAL_BASE="$api" DS_CATALOG_ROOT="$(readlink -f "$cur")" \
      "$PY" "$DS_REPO/bench/benchmark_vs_mpbiopath.py" --max-edges 40000 \
        --ground-truth "$gt" --report "$out/${gt}_report.tsv" --dump-cases "$out/${gt}_cases.tsv" \
        > "$out/${gt}.log" 2>&1 || die "$gt benchmark failed; see $out/${gt}.log"
  done
  "$PY" "$DS_REPO/bench/analysis/holdout_report.py" --cases "$out/curator_cases.tsv" \
      > "$out/holdout.txt" 2>&1 || true
  python3 - "$out" "$build" "$ds" <<'PY'
import json, re, sys, datetime
out, build, ds = sys.argv[1:]
def grab(log):
    t = open(f"{out}/{log}").read()
    acc = re.search(r"end-to-end accuracy:\s+(\d+)/(\d+)", t)
    mf1 = re.search(r"macro-F1 \(DOWN/NORM/UP\):\s+([\d.]+)", t)
    return {"correct": int(acc.group(1)), "cases": int(acc.group(2)),
            "macro_f1": float(mf1.group(1))} if acc and mf1 else None
r = {"catalog_build": build, "solver_commit": ds,
     "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
     "curator": grab("curator.log"), "experimental": grab("experimental.log")}
h = open(f"{out}/holdout.txt").read()
m = re.search(r"HELD-OUT \(report this\)\s+\d+\s+(\d+)\s+([\d.]+)\s+([\d.]+)", h)
if m: r["curator_held_out"] = {"cases": int(m.group(1)), "accuracy": float(m.group(2)), "macro_f1": float(m.group(3))}
json.dump(r, open(f"{out}/RESULTS.json", "w"), indent=2)
print(json.dumps(r, indent=2))
PY
}

cmd_use() {
  local id="${1:?usage: catalog.sh use <build-id>}" dir="$HOME_DIR/builds/$1"
  [ -f "$dir/BUILD.json" ] || die "no build $id"
  [ "$(manifest_field "$dir/BUILD.json" status)" = "complete" ] || die "build $id is not complete"
  ln -sfn "builds/$id" "$HOME_DIR/current"
  echo "catalog: current -> builds/$id"
}

case "${1:-status}" in
  build)  cmd_build ;;
  bench)  cmd_bench ;;
  status) cmd_status ;;
  list)   cmd_list ;;
  use)    shift; cmd_use "$@" ;;
  *) sed -n '2,32p' "$0"; exit 2 ;;
esac
