#!/usr/bin/env python3
"""specs/175: the Lun total readouts under the abundance rule, scored on half A against the specs/127 predictions.

  control: every phospho-readout row of NEW equals OLD (rows, predicted, reach_predicted, pred_ui within 1e-6 relative)
  gate:    moved total-readout cases, fixed vs broken, sign test p < 0.05; recall of measured changes falls <= 0.01

Usage: score_175.py A_FILE OLD_PREDICTIONS NEW_PREDICTIONS OUT_JSON
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
import lun_a_173 as A  # noqa: E402
from common_173 import split_of  # noqa: E402

REL_TOL = 1e-6
MAX_RECALL_DROP = 0.01
ALPHA = 0.05


def rows(path):
    return {(r["pathway_id"], r["gene"], r["readout"]): r for r in csv.DictReader(open(path), delimiter="\t")}


def control(old, new, totals):
    """Differences on phospho readouts (should be none), as a list of keys with a reason."""
    po = {k: v for k, v in old.items() if k[2] not in totals}
    pn = {k: v for k, v in new.items() if k[2] not in totals}
    bad = [(k, "missing in new") for k in po.keys() - pn.keys()] + [(k, "extra in new") for k in pn.keys() - po.keys()]
    for k in po.keys() & pn.keys():
        a, b = po[k], pn[k]
        for f in ("status", "predicted", "reach_predicted"):
            if a[f] != b[f]:
                bad.append((k, f))
        if a["pred_ui"] and b["pred_ui"]:
            x, y = float(a["pred_ui"]), float(b["pred_ui"])
            if abs(x - y) > REL_TOL * max(abs(x), abs(y), 1e-12):
                bad.append((k, "pred_ui"))
        elif a["pred_ui"] != b["pred_ui"]:
            bad.append((k, "pred_ui"))
    return bad


def combined(pred_rows):
    import lun_score as L
    frozen, reach = defaultdict(list), defaultdict(list)
    for r in pred_rows.values():
        if r["self"] == "1" or r["status"] != "scored" or split_of(r["gene"]) != "A":
            continue
        frozen[(r["gene"], r["readout"])].append(r["predicted"])
        reach[(r["gene"], r["readout"])].append(r["reach_predicted"])
    return ({k: L.combine_networks(v) for k, v in frozen.items()}, {k: L.combine_networks(v) for k, v in reach.items()})


def sign_test(fixed, broken):
    """One-sided P(X >= fixed) under Binomial(fixed + broken, 1/2)."""
    n = fixed + broken
    return sum(comb(n, i) for i in range(fixed, n + 1)) / 2 ** n if n else 1.0


def changed(c):
    return c in ("up", "down")


def recall(cases, col):
    real = [c for c in cases if changed(c["measured"])]
    return sum(c[col] == c["measured"] for c in real) / len(real) if real else None


def main():
    import dataset_roles
    import lun_score as L
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, old_f, new_f, out = sys.argv[1:5]
    totals = {r["readout"] for r in A.read_tsv(A.MAPS / "readouts.tsv") if r["kind"] == "total"}
    old, new = rows(old_f), rows(new_f)
    bad = control(old, new, totals)
    res = {"control_differences": len(bad), "control_examples": [f"{k} {why}" for k, why in bad[:20]]}
    genes = A.marker_genes()
    sl = A.slopes(A.load_a(a_file, genes))
    calls = A.measured_calls(sl, genes, A.thresholds(sl))
    fo, ro = combined(old)
    fn, rn = combined(new)
    keys = sorted((fo.keys() | fn.keys()) & calls.keys())
    cases = [{"gene": k[0], "readout": k[1], "measured": calls[k][0],
              "old": fo.get(k, "not"), "new": fn.get(k, "not"),      # a conflict stays a class, as in specs/173
              "reach_old": ro.get(k, "not"), "reach_new": rn.get(k, "not")} for k in keys]
    moved = [c for c in cases if c["old"] != c["new"]]
    fixed = sum(c["new"] == c["measured"] and c["old"] != c["measured"] for c in moved)
    broken = sum(c["old"] == c["measured"] and c["new"] != c["measured"] for c in moved)
    p = sign_test(fixed, broken)
    r_old, r_new = recall(cases, "old"), recall(cases, "new")
    tot = [c for c in cases if c["readout"] in totals]
    ba = {lab: {col: L.balanced_accuracy([(c["measured"], c[col]) for c in sub])
                for col in ("old", "new", "reach_old", "reach_new")}
          for lab, sub in (("all", cases), ("total_readouts", tot))}
    res.update({
        "pairs": len(cases), "pairs_old_only": len(fo.keys() - fn.keys()), "pairs_new_only": len(fn.keys() - fo.keys()),
        "moved": len(moved), "moved_outside_totals": sum(c["readout"] not in totals for c in moved),
        "moved_transitions": dict(Counter(f"{c['old']}->{c['new']} (measured {c['measured']})" for c in moved)),
        "fixed": fixed, "broken": broken, "sign_test_p": p,
        "recall_old": r_old, "recall_new": r_new,
        "gate": {"paired": fixed > broken and p < ALPHA,
                 "sensitivity_guard": r_new is not None and r_old - r_new <= MAX_RECALL_DROP},
        "balanced_accuracy": ba,
        "total_readout_calls": {"old": dict(Counter(c["old"] for c in tot)), "new": dict(Counter(c["new"] for c in tot)),
                                "measured": dict(Counter(c["measured"] for c in tot))},
    })
    res["gate"]["passes"] = all(res["gate"].values()) and not bad
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
