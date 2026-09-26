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
# present and non-empty, and the move is atomic (rename), so a partial build can
# never become what the API serves.
#
# The API container resolves `current` at every container START -- a restart
# after a crash picks up a newly moved `current` with no recreate. `bench`
# therefore resolves the build once and re-checks, after every axis, that the
# server is still serving it.
#
# Usage:
#   scripts/catalog.sh build      generate from LNG HEAD, verify, repoint current
#   scripts/catalog.sh status     what current is, what the RUNNING API serves, staleness
#   scripts/catalog.sh list       every build with its manifest summary
#   scripts/catalog.sh use <id>   point current at an existing complete build
#   scripts/catalog.sh bench      benchmark current on both axes, filed under the
#                                 build by solver commit; refuses unless the
#                                 running API serves that build and the code
#                                 that produces the numbers is committed
set -euo pipefail

HOME_DIR="$(realpath -m "${DS_CATALOG_HOME:-$HOME/deltasignal-catalogs}")"
LNG="${LOGIC_NETWORK_GENERATOR:-$HOME/gitroot/logic-network-generator}"
DS_REPO="$(cd "$(dirname "$0")/.." && pwd)"
LIST="$DS_REPO/bench/catalog_pathways.tsv"
API="${DELTASIGNAL_BASE:-http://localhost:8080}"
ID_RE='^[0-9]{8}-[0-9]{4}_[0-9a-f]{7,40}$'

die() { echo "catalog: $*" >&2; exit 1; }

# Refuse a volatile filesystem by what it IS, not by what its path looks like:
# a literal "/tmp/*" match missed //tmp, relative paths, /run/user, and any
# symlink into a tmpfs.
guard_durable() {
  local probe="$HOME_DIR" fs
  while [ ! -e "$probe" ] && [ "$probe" != "/" ]; do probe="$(dirname "$probe")"; done
  fs="$(stat -f -c %T "$probe" 2>/dev/null || echo unknown)"
  case "$fs" in
    tmpfs|ramfs) die "DS_CATALOG_HOME=$HOME_DIR is on a $fs filesystem; a reboot would wipe it, which is the failure this script exists to prevent" ;;
  esac
}
guard_durable

# The generator's python is resolved only by the commands that need it, so a
# missing generator checkout cannot break `status`, `list` or `use` (rollback).
lng_python() {
  local py="${LNG_PYTHON:-}"
  if [ -z "$py" ]; then
    [ -d "$LNG" ] || die "no generator checkout at $LNG"
    py="$(cd "$LNG" && poetry env info -p 2>/dev/null)/bin/python"
  fi
  [ -x "$py" ] || die "no generator python at $py (set LNG_PYTHON)"
  echo "$py"
}

manifest_field() {
  python3 -c "import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2],''))" "$1" "$2" 2>/dev/null || true
}

# Atomic repoint: build the link beside the target, then rename over it. A bare
# `ln -sfn` unlinks then creates, leaving a window with no `current`.
point_current_at() {
  ln -sfn "builds/$1" "$HOME_DIR/current.tmp.$$"
  mv -Tf "$HOME_DIR/current.tmp.$$" "$HOME_DIR/current"
  echo "catalog: current -> builds/$1"
}

served_build() {
  curl -sf -m 5 "$API/api/health" | python3 -c \
    "import json,sys; print(json.load(sys.stdin).get('catalog',{}).get('build_id',''))" 2>/dev/null || true
}

cmd_build() {
  # Options (experiments only): --variant NAME builds into <id>_NAME and never
  # moves `current`; --env LNG_X=v passes a generator flag (recorded in
  # BUILD.json). A plain build moves `current` as before.
  local variant="" extra_env=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --variant) [ -n "${2:-}" ] || die "--variant needs a name"; variant="$2"; shift 2 ;;
      --env) [[ "$2" == LNG_*=* ]] || die "--env takes LNG_NAME=value"; extra_env+=("$2"); shift 2 ;;
      *) die "unknown build option $1" ;;
    esac
  done
  [ ${#extra_env[@]} -eq 0 ] || [ -n "$variant" ] || die "--env needs --variant: a flagged build must not become current"
  [ -z "$variant" ] || [[ "$variant" =~ ^[A-Za-z0-9._-]+$ ]] || die "--variant must be [A-Za-z0-9._-]+"
  # Generator flags reach a build only through --env, so they are recorded.
  # LNG_PYTHON only selects the interpreter (lng_python), not a generator flag.
  local leaked; leaked=$( (env | grep -oE '^LNG_[A-Z0-9_]+' | grep -vx LNG_PYTHON | sort -u | tr '\n' ' ') || true)
  [ -z "$leaked" ] || die "the calling shell exports $leaked-- unset them, or pass them with --env"
  [ -d "$LNG/.git" ] || die "no generator checkout at $LNG"
  local py sha dirty branch id dir
  py="$(lng_python)"
  sha="$(git -C "$LNG" rev-parse --short HEAD)"
  branch="$(git -C "$LNG" rev-parse --abbrev-ref HEAD)"
  dirty="$(git -C "$LNG" status --porcelain src bin pyproject.toml poetry.lock | wc -l | tr -d ' ')"
  [ "$dirty" = "0" ] || die "generator has $dirty uncommitted change(s) in src/, bin/ or its dependency files; a build must be reproducible from a commit"
  id="$(date +%Y%m%d-%H%M)_${sha}${variant:+_$variant}"
  dir="$HOME_DIR/builds/$id"
  mkdir -p "$HOME_DIR/builds"
  # No -p on the leaf: mkdir is the atomic existence check, so two builds in the
  # same minute cannot both claim the directory.
  mkdir "$dir" 2>/dev/null || die "$dir already exists (another build this minute?)"

  grep -v '^#' "$LIST" > "$dir/pathways.tsv"
  local requested; requested=$(($(wc -l < "$dir/pathways.tsv") - 1))

  echo "catalog: building $id ($requested pathways) into $dir"
  local t0 rc; t0=$(date +%s)
  set +e
  ( cd "$LNG" && env PYTHONHASHSEED=0 LNG_EMIT_ONE_SIDED=1 "${extra_env[@]}" \
      "$py" bin/create-pathways.py --pathway-list "$dir/pathways.tsv" --output-dir "$dir" ) \
      > "$dir/generate.log" 2>&1
  rc=$?
  set -e

  EXTRA_ENV="$(printf '%s\n' "${extra_env[@]}")" VARIANT="$variant" \
  "$py" - "$dir" "$sha" "$branch" "$rc" "$requested" "$(( $(date +%s) - t0 ))" "$LIST" "$LNG" <<'PY' || die "could not write the manifest for $id; current left unchanged"
import csv, glob, hashlib, json, os, sys, datetime, subprocess
d, sha, branch, rc, requested, secs, lst, lng = sys.argv[1:]

def network_ok(path):
    """A network counts as built only with the expected header and >= 1 edge.
    Existence alone passed a 0-byte or truncated file as complete."""
    try:
        with open(path, newline="") as fh:
            r = csv.DictReader(fh)
            if not {"source_id", "target_id"} <= set(r.fieldnames or []):
                return False
            return next(r, None) is not None
    except (OSError, csv.Error):
        return False

wanted = [r["id"] for r in csv.DictReader(open(os.path.join(d, "pathways.tsv")), delimiter="\t")]
built, bad = [], []
for p in sorted(os.path.basename(x) for x in glob.glob(os.path.join(d, "R-HSA-*"))):
    (built if network_ok(os.path.join(d, p, "logic_network.csv")) else bad).append(p)

# Summed per pathway, the way every earlier measurement in this repo counted.
nodes = edges = 0
for p in built:
    rows = list(csv.DictReader(open(os.path.join(d, p, "logic_network.csv"), newline="")))
    edges += len(rows)
    nodes += len({r["source_id"] for r in rows} | {r["target_id"] for r in rows})

# The same generator commit against a different Reactome release builds a
# different catalog, so the database identity is part of the provenance.
reactome = {}
try:
    from py2neo import Graph
    info = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"),
                 auth=(os.environ.get("NEO4J_USER", "neo4j"),
                       os.environ.get("NEO4J_PASSWORD", ""))).run(
        "MATCH (d:DBInfo) RETURN properties(d) AS p LIMIT 1").evaluate() or {}
    reactome = {k: info.get(k) for k in ("releaseNumber", "releaseDate", "checksum", "neo4j")}
except Exception as e:
    reactome = {"error": f"could not read DBInfo: {e}"}
try:
    reactome["image"] = subprocess.run(
        ["docker", "inspect", "reactome-neo4j", "--format", "{{.Config.Image}}"],
        capture_output=True, text=True, timeout=10).stdout.strip() or None
except Exception:
    pass

def sha256(path):
    try:
        return hashlib.sha256(open(path, "rb").read()).hexdigest()
    except OSError:
        return None

missing = sorted(set(wanted) - set(built))
status = "complete" if (rc == "0" and not missing and not bad) else "incomplete"
m = {
    "build_id": os.path.basename(d),
    "created": datetime.datetime.now().isoformat(timespec="seconds"),
    "status": status,
    "generator_commit": sha,
    "generator_branch": branch,
    "generator_exit_code": int(rc),
    "generator_poetry_lock_sha256": sha256(os.path.join(lng, "poetry.lock")),
    "reactome": reactome,
    "env": {"PYTHONHASHSEED": "0", "LNG_EMIT_ONE_SIDED": "1",
            **dict(kv.split("=", 1) for kv in os.environ.get("EXTRA_ENV", "").splitlines() if "=" in kv)},
    "variant": os.environ.get("VARIANT") or None,
    "pathway_list": os.path.relpath(lst, os.path.dirname(os.path.dirname(lst))),
    "pathway_list_sha256": sha256(os.path.join(d, "pathways.tsv")),
    "pathways_requested": int(requested),
    "pathways_built": len(built),
    "pathways_missing": missing,
    "pathways_empty_or_truncated": bad,
    "nodes": nodes,
    "nodes_note": "summed per pathway, as every earlier measurement counted",
    "edges": edges,
    "seconds": int(secs),
}
# Written atomically, so a crash cannot leave a half-written manifest.
tmp = os.path.join(d, "BUILD.json.tmp")
with open(tmp, "w") as fh:
    json.dump(m, fh, indent=2)
os.replace(tmp, os.path.join(d, "BUILD.json"))
print(json.dumps({k: m[k] for k in ("status", "pathways_built", "pathways_requested", "nodes", "edges")}))
PY

  if [ "$(manifest_field "$dir/BUILD.json" status)" != "complete" ]; then
    echo "catalog: build $id is INCOMPLETE; current left unchanged. See $dir/BUILD.json and generate.log" >&2
    exit 1
  fi
  if [ -n "$variant" ]; then
    echo "catalog: variant build $id is complete; current NOT moved (run arms with scripts/run_arm.sh --catalog $id)"
    return
  fi
  point_current_at "$id"
  echo "catalog: recreate the API to serve it: docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api"
}

cmd_status() {
  local cur="$HOME_DIR/current"
  if [ ! -L "$cur" ]; then
    echo "no current catalog (run: scripts/catalog.sh build)"
    [ -e "$cur" ] && echo "WARNING: $cur exists but is not a symlink -- likely an empty directory docker created for a missing mount; remove it"
    exit 1
  fi
  local id m head served
  id="$(basename "$(readlink "$cur")")"
  m="$cur/BUILD.json"
  head="$(git -C "$LNG" rev-parse --short HEAD 2>/dev/null || echo '?')"
  echo "current:   $id"
  echo "status:    $(manifest_field "$m" status)"
  echo "generator: $(manifest_field "$m" generator_commit)  (generator HEAD now: $head)"
  echo "reactome:  $(python3 -c "import json,sys; r=json.load(open(sys.argv[1])).get('reactome',{}); print(f\"release {r.get('releaseNumber','?')} ({r.get('releaseDate','?')})\")" "$m" 2>/dev/null || echo '?')"
  echo "built:     $(manifest_field "$m" created)"
  echo "size:      $(manifest_field "$m" pathways_built) pathways, $(manifest_field "$m" nodes) nodes, $(manifest_field "$m" edges) edges"
  if [ "$head" != "?" ] && [ -n "$(manifest_field "$m" generator_commit)" ] && [ "$head" != "$(manifest_field "$m" generator_commit)" ]; then
    local changed
    changed="$(git -C "$LNG" diff --stat "$(manifest_field "$m" generator_commit)" HEAD -- src bin 2>/dev/null | tail -1)"
    echo "STALE:     generator has moved since this build (${changed:-no src/bin change})"
  fi
  # What the RUNNING service reports, not what docker was asked to mount: the
  # mount source is the literal `current` path and says nothing about the build.
  served="$(served_build)"
  if [ -z "$served" ]; then
    echo "dev API:   not reachable at $API"
  elif [ "$served" = "$id" ]; then
    echo "dev API:   serving $served (matches current)"
  else
    echo "dev API:   serving $served -- NOT current; recreate it: docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api"
  fi
}

cmd_list() {
  local cur_target any=0
  cur_target="$(readlink "$HOME_DIR/current" 2>/dev/null || true)"
  for dir in "$HOME_DIR"/builds/*/; do
    [ -d "$dir" ] || continue
    any=1
    local id; id="$(basename "$dir")"
    if [ -f "$dir/BUILD.json" ]; then
      printf '%-26s %-10s %3s pathways  %s\n' "$id" "$(manifest_field "$dir/BUILD.json" status)" \
        "$(manifest_field "$dir/BUILD.json" pathways_built)" \
        "$([ "$cur_target" = "builds/$id" ] && echo '<- current')"
    else
      printf '%-26s %-10s               %s\n' "$id" "NO MANIFEST" "(interrupted or failed build)"
    fi
  done
  [ "$any" = 1 ] || echo "no builds yet"
}

cmd_use() {
  local id="${1:?usage: catalog.sh use <build-id>}"
  [[ "$id" =~ $ID_RE ]] || die "'$id' is not a build id (expected <YYYYmmdd-HHMM>_<sha>)"
  local dir="$HOME_DIR/builds/$id"
  [ -f "$dir/BUILD.json" ] || die "no build $id"
  [ "$(manifest_field "$dir/BUILD.json" status)" = "complete" ] || die "build $id is not complete"
  point_current_at "$id"
}

cmd_bench() {
  local cur="$HOME_DIR/current"
  [ -L "$cur" ] || die "no current catalog (run: scripts/catalog.sh build)"
  local py target build ds dirty
  py="$(lng_python)"   # the benchmark needs py2neo, which lives in the generator venv
  # Resolve ONCE. Reading through the symlink per axis let a mid-run repoint
  # send the second axis to a different build.
  target="$(readlink -f "$cur")"
  build="$(basename "$target")"
  ds="$(git -C "$DS_REPO" rev-parse --short HEAD)"
  # Everything that produces the numbers must be committed, not just src/.
  dirty="$(git -C "$DS_REPO" status --porcelain src bench | wc -l | tr -d ' ')"
  [ "$dirty" = "0" ] || die "deltasignal src/ or bench/ has $dirty uncommitted change(s); results must come from a commit"

  check_served() {
    local s; s="$(served_build)"
    [ -n "$s" ] || die "cannot read the served catalog from $API/api/health (is the API up?)"
    [ "$s" = "$build" ] || die "API at $API serves '$s' but this run is scoring '$build'${1:+ ($1)}; recreate it: docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api"
  }
  check_served ""

  local out="$target/results/$ds"
  mkdir -p "$out"
  echo "catalog: benchmarking build $build with solver $ds -> $out"
  for gt in curator experimental; do
    env DELTASIGNAL_BASE="$API" DS_CATALOG_ROOT="$target" \
      "$py" "$DS_REPO/bench/benchmark_vs_mpbiopath.py" --max-edges 40000 \
        --ground-truth "$gt" --report "$out/${gt}_report.tsv" --dump-cases "$out/${gt}_cases.tsv" \
        > "$out/${gt}.log" 2>&1 || die "$gt benchmark failed; see $out/${gt}.log"
    # The container re-resolves `current` on any restart, so a crash mid-run
    # after a repoint would silently switch the served build.
    check_served "it changed during the $gt axis -- the container restarted?"
  done
  "$py" "$DS_REPO/bench/analysis/holdout_report.py" --cases "$out/curator_cases.tsv" \
      > "$out/holdout.txt" 2>&1 || die "holdout report failed; see $out/holdout.txt"

  python3 - "$out" "$build" "$ds" <<'PY'
import json, re, sys, datetime, os
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
# The perturbation protocol is part of the result (specs/023: pins had silently
# become every complex containing the gene for two months).
cur = open(f"{out}/curator.log").read()
for key, pat in (("protocol", r"^Protocol: (.*)$"), ("pinned", r"^Pinned: (.*)$")):
    m = re.search(pat, cur, re.M)
    r[key] = m.group(1) if m else None
h = open(f"{out}/holdout.txt").read()
m = re.search(r"HELD-OUT \(report this\)\s+\d+\s+(\d+)\s+([\d.]+)\s+([\d.]+)", h)
r["curator_held_out"] = ({"cases": int(m.group(1)), "accuracy": float(m.group(2)),
                          "macro_f1": float(m.group(3))} if m else None)
# A result with a missing metric is not a result. Recording nulls and exiting 0
# would look like a successful run.
missing = [k for k in ("curator", "experimental", "curator_held_out", "protocol") if r[k] is None]
if missing:
    sys.exit(f"catalog: could not parse {', '.join(missing)} from the benchmark output; "
             f"the output format may have changed. Nothing recorded as RESULTS.json.")
tmp = f"{out}/RESULTS.json.tmp"
json.dump(r, open(tmp, "w"), indent=2)
os.replace(tmp, f"{out}/RESULTS.json")
print(json.dumps(r, indent=2))
PY
}

case "${1:-status}" in
  build)  shift; cmd_build "$@" ;;
  bench)  cmd_bench ;;
  status) cmd_status ;;
  list)   cmd_list ;;
  use)    shift; cmd_use "$@" ;;
  *) sed -n '2,40p' "$0"; exit 2 ;;
esac
