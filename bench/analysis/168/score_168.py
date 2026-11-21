#!/usr/bin/env python3
"""specs/168: score the single-copy and pilot literature arms (the specs/166 numbers, reused).

Comparisons (each over the cases scored in both arms of the pair, per axis):
  GATE   v3-single vs v3-without: primary net (korkut, own markers excluded, + experimental) > 0, and
         above every finished shuffle of v3-single (shuffles run only if the net is > 0).
  COPY   v3-single vs v3: the effect of attaching each literature edge once rather than to every copy.
  DESC   v3, v3-strict (in-sample selection: chosen from specs/166 slices) vs v3-without;
         pilot, pilot-strict, pilot-single vs pilot-without (three pathways; papers searched with
         korkut drug-target/readout gene pairs, so korkut reach is targeted, labels unused).
Per pair and axis: the specs/166 three numbers (coverage gained, newly covered, previously covered)
and fixed / broken / net, plus change calls and how many are right.

Usage: score_168.py BUILDS_ROOT SHA OUT_JSON
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))
sys.path.insert(0, str(HERE.parents[1] / "166"))
import score_166 as s  # noqa: E402

RUN = "lit168"
V3, PILOT = "20261120-lit-v3", "20261120-lit-pilot"
PAIRS = {  # name: (base arm, compared arm)
    "gate": (f"{V3}-without", f"{V3}-single-with"),
    "copy": (f"{V3}-with", f"{V3}-single-with"),
    "v3": (f"{V3}-without", f"{V3}-with"),
    "v3_strict": (f"{V3}-without", f"{V3}-strict-with"),
    "pilot": (f"{PILOT}-without", f"{PILOT}-with"),
    "pilot_strict": (f"{PILOT}-without", f"{PILOT}-strict-with"),
    "pilot_single": (f"{PILOT}-without", f"{PILOT}-single-with"),
}
SHUF = f"{V3}-single-shuf"


def change_calls(arm, keys):
    calls = [k for k in keys if arm[k][0] != "not"]
    return {"n": len(calls), "right": sum(1 for k in calls if arm[k][0] == arm[k][1])}


def compare(base, arm):
    out = {}
    for ax in sorted(set(base) & set(arm)):
        keys = sorted(set(base[ax]) & set(arm[ax]))
        r = s.three_numbers(base[ax], arm[ax], keys)
        r["change_calls_base"] = change_calls(base[ax], keys)
        r["change_calls_arm"] = change_calls(arm[ax], keys)
        out[ax] = r
    out["primary_net"] = sum(out[ax]["all"]["net"] for ax in s.PRIMARY if ax in out)
    return out


def gate(net_with, shuffle_nets):
    if net_with <= 0:
        return {"net_with": net_with, "pass": False, "shuffles": "not run: net <= 0 (registered)"}
    return s.gate(net_with, shuffle_nets)


def main():
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    import dataset_roles
    dataset_roles.require_role("korkut2015", "fit")
    builds, sha, out_json = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    marks, tuning = marker_names(), set(TUNING_PATHWAYS)
    cache = {}

    def arm(name, labs=("curator", "experimental")):
        if name not in cache:
            run = builds / name / "results" / sha / RUN
            if not (run / "korkut_cases.tsv").exists():
                return None
            cache[name] = s.load_arm(run, tuning, marks, labs=labs)
        return cache[name]

    res = {"pairs": {}, "missing": []}
    for name, (a, b) in PAIRS.items():
        base, other = arm(a), arm(b)
        if base is None or other is None:
            res["missing"].append(name)
            continue
        res["pairs"][name] = compare(base, other)
    shuf = []
    for d in sorted(builds.glob(SHUF + "*")):
        o = arm(d.name, labs=("experimental",))
        if o is not None and "gate" in res["pairs"]:
            shuf.append(compare(arm(PAIRS["gate"][0]), o)["primary_net"])
    if "gate" in res["pairs"]:
        res["gate"] = gate(res["pairs"]["gate"]["primary_net"], shuf)
    out_json.write_text(json.dumps(res, indent=1, default=str))
    for name, r in res["pairs"].items():
        print(f"{name:>13}: primary net {r['primary_net']:+d}  " + "  ".join(
            f"{ax} {r[ax]['all']['net']:+d} (fixed {r[ax]['all']['fixed']}, broken {r[ax]['all']['broken']}, "
            f"gain {r[ax]['coverage_gained']}, calls {r[ax]['change_calls_base']['n']}->{r[ax]['change_calls_arm']['n']})"
            for ax in sorted(r) if isinstance(r[ax], dict) and "all" in r[ax]))
    print("GATE:", res.get("gate"), "missing:", res["missing"])


if __name__ == "__main__":
    main()
