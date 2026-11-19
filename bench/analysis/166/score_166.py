#!/usr/bin/env python3
"""specs/166: score the literature-edge arms against `without`.

The registered gate is on the primary pool (korkut with own markers excluded,
plus MP-BioPath experimental), over the cases scored in every primary arm
(`with`, `without`, all shuffles):
  net(X) = #(correct in X, wrong in without) - #(wrong in X, correct in without)
  PASS iff net(with) > 0 and net(with) > net(shufN) for every seed
  (rank p = (1 + #{seeds with net >= net(with)}) / (1 + seeds) <= 0.05).
Also reported: the three numbers per axis (coverage gained, newly covered,
previously covered) for every arm, and the `only-<slice>` ablations (descriptive).

Usage: score_166.py BUILDS_ROOT SHA UNREACHED_DIR OUT_JSON
(results are read from BUILDS_ROOT/20261118-lit-h12b-<arm>/results/SHA/lit166/)
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))

PREFIX = "20261118-lit-h12b-"
RUN = "lit166"
CODE = {"0": "down", "1": "not", "2": "up"}
PRIMARY = ("korkut", "experimental")
ALPHA = 0.05


def load_mpb(path, tuning, lab):
    """{axis: {(pathway, gene, direction, key_output): (call, truth, fold)}} for valid cases."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r.get("valid") != "1":
            continue
        if lab == "curator":
            axis = "curator_tuning" if r["pathway"] in tuning else "curator_heldout"
        else:
            axis = "experimental"
        k = (r["pathway"], r["gene"], r["direction"], r["key_output"])
        out.setdefault(axis, {})[k] = (CODE[r["predicted"]], CODE[r["expected"]], float(r["pred_ui"]))
    return out


def load_korkut(path, marks):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "scored":
            continue
        if any(r["readout"] in marks.get(d, set()) for d in r["drug_set"].split("+")):
            continue
        out[(r["pathway_id"], r["drug_set"], r["readout"])] = (r["predicted"], r["measured"], float(r["pred_ui"]))
    return out


def load_arm(run_dir, tuning, marks):
    axes = {}
    for lab in ("curator", "experimental"):
        axes.update(load_mpb(run_dir / f"{lab}_cases.tsv", tuning, lab))
    axes["korkut"] = load_korkut(run_dir / "korkut_cases.tsv", marks)
    return axes


def net(base, arm, keys):
    fixed = sum(1 for k in keys if arm[k][0] == arm[k][1] and base[k][0] != base[k][1])
    broken = sum(1 for k in keys if arm[k][0] != arm[k][1] and base[k][0] == base[k][1])
    return {"fixed": fixed, "broken": broken, "net": fixed - broken}


def three_numbers(base, arm, keys):
    new = [k for k in keys if base[k][2] == 1.0 and arm[k][2] != 1.0]
    lost = [k for k in keys if base[k][2] != 1.0 and arm[k][2] == 1.0]
    prev = [k for k in keys if base[k][2] != 1.0]
    chg = [k for k in new if arm[k][0] != "not"]
    acc = lambda ks, a: (sum(a[k][0] == a[k][1] for k in ks) / len(ks)) if ks else None  # noqa: E731
    return {
        "cases": len(keys),
        "unreached_without": sum(1 for k in keys if base[k][2] == 1.0),
        "coverage_gained": len(new), "coverage_lost": len(lost),
        "new": {"accuracy": acc(new, arm), "accuracy_without": acc(new, base),
                "change_calls": len(chg),
                "change_precision": (sum(arm[k][0] == arm[k][1] for k in chg) / len(chg)) if chg else None,
                "calls": dict(Counter(arm[k][0] for k in new)),
                "truth": dict(Counter(arm[k][1] for k in new))},
        "previous": {"n": len(prev), **net(base, arm, prev),
                     "accuracy_without": acc(prev, base), "accuracy_with": acc(prev, arm)},
        "all": {**net(base, arm, keys), "accuracy_without": acc(keys, base), "accuracy_with": acc(keys, arm)},
    }


def gate(nets_real, nets_shuf):
    ge = sum(1 for n in nets_shuf if n >= nets_real)
    p = (1 + ge) / (1 + len(nets_shuf))
    return {"net_with": nets_real, "net_shuffles": sorted(nets_shuf), "seeds_at_or_above": ge,
            "rank_p": p, "net_positive": nets_real > 0,
            "beats_every_seed": ge == 0,
            "pass": nets_real > 0 and ge == 0 and p <= ALPHA}


def common_keys(arms, axis):
    ks = None
    for a in arms.values():
        s = set(a.get(axis, {}))
        ks = s if ks is None else ks & s
    return sorted(ks or ())


def unreached_gains(base, arm, keys, lists):
    """Gains on the lists handed over on 2026-11-18."""
    kk, ce = lists
    out = {}
    korkut_keys = [k for k in keys["korkut"] if k in kk]
    out["korkut_list_scored"] = len(korkut_keys)
    out["korkut_list_gained"] = sum(1 for k in korkut_keys
                                    if base["korkut"][k][2] == 1.0 and arm["korkut"][k][2] != 1.0)
    n = g = 0
    for axis in ("curator_heldout", "curator_tuning", "experimental"):
        for k in keys[axis]:
            if (k[0], k[1], k[3]) in ce:
                n += 1
                g += base[axis][k][2] == 1.0 and arm[axis][k][2] != 1.0
    out["curator_experimental_list_scored"] = n
    out["curator_experimental_list_gained"] = g
    return out


def load_lists(d):
    kk = {(r["pathway_stid"], r["drug_set"], r["readout"])
          for r in csv.DictReader(open(d / "unreached_korkut.tsv"), delimiter="\t")}
    ce = {(r["pathway"], r["perturbed_gene"], r["readout_dbid"])
          for r in csv.DictReader(open(d / "unreached_curator_experimental.tsv"), delimiter="\t")}
    return kk, ce


def score(arms, lists=None):
    """arms: {name: {axis: {key: (call, truth, fold)}}}; names with, without, shufNN, only-*."""
    primary = {n: a for n, a in arms.items() if n in ("with", "without") or n.startswith("shuf")}
    axes = sorted({ax for a in arms.values() for ax in a})
    keys = {ax: common_keys(primary, ax) for ax in axes}
    base = arms["without"]
    res = {"denominators": {n: {ax: len(a.get(ax, {})) for ax in axes} for n, a in arms.items()},
           "common": {ax: len(keys[ax]) for ax in axes}, "arms": {}}
    for n, a in arms.items():
        if n == "without":
            continue
        ks = keys if n in primary else {ax: common_keys({"w": base, "x": a}, ax) for ax in axes}
        r = {ax: three_numbers(base[ax], a[ax], ks[ax]) for ax in axes}
        r["primary_net"] = sum(r[ax]["all"]["net"] for ax in PRIMARY if ax in r)
        if lists is not None:
            r["unreached_lists"] = unreached_gains(base, a, ks, lists)
        res["arms"][n] = r
    shuf = [res["arms"][n]["primary_net"] for n in sorted(res["arms"]) if n.startswith("shuf")]
    res["gate"] = gate(res["arms"]["with"]["primary_net"], shuf)
    return res


def main():
    import dataset_roles
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    dataset_roles.require_role("korkut2015", "fit")
    builds, sha, lists_dir, out_json = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
    marks = marker_names()
    arms = {}
    for d in sorted(builds.glob(PREFIX + "*")):
        run = d / "results" / sha / RUN
        if not (run / "korkut_cases.tsv").exists():
            raise SystemExit(f"{d.name}: no completed run at {run}")
        arms[d.name[len(PREFIX):]] = load_arm(run, set(TUNING_PATHWAYS), marks)
    res = score(arms, load_lists(lists_dir))
    out_json.write_text(json.dumps(res, indent=1, default=str))
    g = res["gate"]
    print(f"common cases: {res['common']}")
    for n in ["with"] + sorted(x for x in res["arms"] if x != "with"):
        r = res["arms"][n]
        print(f"{n:>18}: primary net {r['primary_net']:+d}  " + "  ".join(
            f"{ax} {r[ax]['all']['net']:+d} (gain {r[ax]['coverage_gained']})" for ax in sorted(r)
            if isinstance(r[ax], dict) and "all" in r[ax]))
    print(f"GATE: net(with) {g['net_with']:+d}; shuffles {g['net_shuffles']}; rank p {g['rank_p']:.3f} -> "
          f"{'PASS' if g['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
