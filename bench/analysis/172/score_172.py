#!/usr/bin/env python3
"""specs/172: score the v4 literature arms under DS_PIN_CURATED_ROOT (run name lit172).

  PARITY    v4-single-without (lit172) against v3-without (lit170): every class identical (FR-002), checked first.
  GATE      v4-strict-single-with vs v4-single-without (both lit172): primary net > 0 and above every finished shuffle
            (v4-strict-single-shufNN under lit172, run only if the net is > 0).
  SINGLE    v4-single-with vs v4-single-without (both lit172), described.
  TIES      v3-single-with (lit170) vs v4-single-with (lit172): what the attachment tie-breaking moved, described.
Per pair: the specs/168 numbers (specs/166 three numbers, fixed / broken / net, change calls).

Usage: score_172.py BUILDS_ROOT SHA_170 SHA_172 OUT_JSON
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
for d in ("161", "166", "168", "170"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import score_166 as s  # noqa: E402
import score_168 as s168  # noqa: E402
from score_170 import parity  # noqa: E402

V4 = "20261123-lit-v4"
CONTROL = (f"{V4}-single-without", "new")
PAIRS = {  # name: ((base arm, run), (compared arm, run)); run "old" = lit170, "new" = lit172
    "gate": (CONTROL, (f"{V4}-strict-single-with", "new")),
    "single": (CONTROL, (f"{V4}-single-with", "new")),
    "ties": (("20261120-lit-v3-single-with", "old"), (f"{V4}-single-with", "new")),
}
SHUF = f"{V4}-strict-single-shuf"


def main():
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    import dataset_roles
    dataset_roles.require_role("korkut2015", "fit")
    builds, sha_old, sha_new, out_json = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
    runs = {"old": (sha_old, "lit170"), "new": (sha_new, "lit172")}
    marks, tuning = marker_names(), set(TUNING_PATHWAYS)
    cache = {}

    def arm(name, which, labs=("curator", "experimental")):
        key = (name, which)
        if key not in cache:
            sha, run = runs[which]
            d = builds / name / "results" / sha / run
            if not (d / "korkut_cases.tsv").exists():
                return None
            cache[key] = s.load_arm(d, tuning, marks, labs=labs)
        return cache[key]

    res = {"pairs": {}, "missing": []}
    old_c, new_c = arm("20261120-lit-v3-without", "old"), arm(*CONTROL)
    if old_c is None or new_c is None:
        raise SystemExit("a control run is missing; parity (FR-002) cannot be checked")
    res["parity"] = parity(old_c, new_c)
    if not res["parity"]["identical"]:
        out_json.write_text(json.dumps(res, indent=1, default=str))
        raise SystemExit("FR-002: the control does not reproduce lit170; trace before scoring")
    for name, (a, b) in PAIRS.items():
        base, other = arm(*a), arm(*b)
        if base is None or other is None:
            res["missing"].append(name)
            continue
        res["pairs"][name] = s168.compare(base, other)
    shuf = []
    for d in sorted(builds.glob(SHUF + "*")):
        o = arm(d.name, "new", labs=("experimental",))
        if o is not None and "gate" in res["pairs"]:
            shuf.append(s168.compare(new_c, o)["primary_net"])
    if "gate" in res["pairs"]:
        res["gate"] = s168.gate(res["pairs"]["gate"]["primary_net"], shuf)
    out_json.write_text(json.dumps(res, indent=1, default=str))
    for name, r in res["pairs"].items():
        print(f"{name:>6}: primary net {r['primary_net']:+d}  " + "  ".join(
            f"{ax} {r[ax]['all']['net']:+d} (fixed {r[ax]['all']['fixed']}, broken {r[ax]['all']['broken']}, "
            f"gain {r[ax]['coverage_gained']}, calls {r[ax]['change_calls_base']['n']}->{r[ax]['change_calls_arm']['n']})"
            for ax in sorted(r) if isinstance(r[ax], dict) and "all" in r[ax]))
    print("PARITY identical:", res["parity"]["identical"], " GATE:", res.get("gate"), " missing:", res["missing"])


if __name__ == "__main__":
    main()
