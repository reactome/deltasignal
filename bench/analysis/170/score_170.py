#!/usr/bin/env python3
"""specs/170: score the specs/168 v3 arms re-run under DS_PIN_CURATED_ROOT (run name lit170).

  PARITY    v3-without under lit170 against lit168: every class identical (FR-002), checked before anything else.
  GATE      v3-single-with vs v3-without (both lit170): primary net > 0 and above every finished shuffle
            (v3-single-shufNN under lit170, run only if the net is > 0).
  PROTOCOL  v3-single-with lit168 vs lit170: how much of the specs/168 loss was the pins.
  V3        v3-with vs v3-without (both lit170), described.
Per pair: the specs/168 numbers (specs/166 three numbers, fixed / broken / net, change calls).

Usage: score_170.py BUILDS_ROOT SHA_168 SHA_170 OUT_JSON
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))
sys.path.insert(0, str(HERE.parents[1] / "166"))
sys.path.insert(0, str(HERE.parents[1] / "168"))
import score_166 as s  # noqa: E402
import score_168 as s168  # noqa: E402

V3 = "20261120-lit-v3"
PAIRS = {  # name: ((base arm, run), (compared arm, run)); run "new" = lit170, "old" = lit168
    "gate": ((f"{V3}-without", "new"), (f"{V3}-single-with", "new")),
    "protocol": ((f"{V3}-single-with", "old"), (f"{V3}-single-with", "new")),
    "v3": ((f"{V3}-without", "new"), (f"{V3}-with", "new")),
}
SHUF = f"{V3}-single-shuf"


def parity(old, new):
    """Cases whose class differs between two runs of one arm, per axis; also cases scored in only one."""
    out = {}
    for ax in sorted(set(old) | set(new)):
        a, b = old.get(ax, {}), new.get(ax, {})
        out[ax] = {"cases": len(set(a) & set(b)), "only_old": len(set(a) - set(b)), "only_new": len(set(b) - set(a)),
                   "class_differs": sorted(str(k) for k in set(a) & set(b) if a[k][0] != b[k][0])}
    out["identical"] = all(v["only_old"] == v["only_new"] == 0 and not v["class_differs"]
                           for k, v in out.items() if isinstance(v, dict))
    return out


def main():
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    import dataset_roles
    dataset_roles.require_role("korkut2015", "fit")
    builds, sha_old, sha_new, out_json = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
    runs = {"old": (sha_old, "lit168"), "new": (sha_new, "lit170")}
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
    old_c, new_c = arm(f"{V3}-without", "old"), arm(f"{V3}-without", "new")
    if old_c is None or new_c is None:
        raise SystemExit("control missing under one of the runs; parity (FR-002) cannot be checked")
    res["parity"] = parity(old_c, new_c)
    if not res["parity"]["identical"]:
        out_json.write_text(json.dumps(res, indent=1, default=str))
        raise SystemExit("FR-002: the control does not reproduce lit168; trace before scoring")
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
        print(f"{name:>9}: primary net {r['primary_net']:+d}  " + "  ".join(
            f"{ax} {r[ax]['all']['net']:+d} (fixed {r[ax]['all']['fixed']}, broken {r[ax]['all']['broken']}, "
            f"gain {r[ax]['coverage_gained']}, calls {r[ax]['change_calls_base']['n']}->{r[ax]['change_calls_arm']['n']})"
            for ax in sorted(r) if isinstance(r[ax], dict) and "all" in r[ax]))
    print("PARITY identical:", res["parity"]["identical"], " GATE:", res.get("gate"), " missing:", res["missing"])


if __name__ == "__main__":
    main()
