#!/usr/bin/env bash
# specs/149 FR-004: package a serving catalog = the frozen build plus the lead
# composites under their non-listed names. Run, then point the API at OUT.
#   specs/149-leads-serving/package_leads_catalog.sh OUT_DIR
set -euo pipefail
FROZEN="$HOME/deltasignal-catalogs/builds/20261026-1349_8b6527e"
COMPOSITES="$HOME/deltasignal-catalogs/builds/20261107-arrows97-11c"
OUT="${1:?usage: package_leads_catalog.sh OUT_DIR}"
[ -e "$OUT" ] && { echo "$OUT exists" >&2; exit 1; }
mkdir -p "$OUT"
for d in "$FROZEN"/R-HSA-* "$FROZEN"/BUILD.json; do cp -r "$d" "$OUT/"; done
for d in "$COMPOSITES"/R-HSA-*; do cp -r "$d" "$OUT/ctx-arrows_$(basename "$d")"; done
python3 - "$OUT" <<'PY'
import json, sys
p = sys.argv[1] + "/BUILD.json"
b = json.load(open(p))
b["leads_composites"] = "20261107-arrows97-11c (specs/149; ctx-arrows_ prefix, unlisted)"
json.dump(b, open(p, "w"), indent=1)
PY
echo "packaged -> $OUT"
