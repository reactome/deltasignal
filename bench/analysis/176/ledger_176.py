#!/usr/bin/env python3
"""specs/176: one non-overlapping ledger of the half-A misses and wrong calls (each case counted once, first cause
that applies, in the order listed). Usage: ledger_176.py CASES_176 AXIS_JSON TRACES_ALL A_FILE OUT_JSON"""
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402


def main():
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    cases_f, axis_f, traces, a_file, out = sys.argv[1:6]
    cases = list(csv.DictReader(open(cases_f), delimiter="\t"))
    axis = json.load(open(axis_f))
    tk, res = set(axis["tyrosine_kinases"]), axis["readout_residues"]
    m0 = B.measured(a_file, A.marker_genes())
    n_changed = Counter()
    for (g, _), (e, d) in m0.items():
        n_changed[g] += B.changed(B.call_of(e, d))
    best = {}
    for r in csv.DictReader(open(traces), delimiter="\t"):
        if r["hops"] != "":
            k = (r["gene"], r["readout"])
            if k not in best or int(r["hops"]) < best[k][0]:
                best[k] = (int(r["hops"]), B.PIP3 in r["path"])
    miss, wrong = Counter(), Counter()
    for c in cases:
        k = (c["gene"], c["readout"])
        ch, called = B.changed(c["measured"]), B.changed(c["frozen"])
        others = n_changed[c["gene"]] - ch
        if ch and c["frozen"] != c["measured"]:
            if abs(float(c["E"])) / float(c["delta"]) < 1.5:
                miss["1 marginal effect (under 1.5 delta; about 40% do not repeat)"] += 1
            elif c["gene"] in tk and "Y" in res.get(c["readout"], "") and c["measured"] == "up":
                miss["2 tyrosine kinase raising a tyrosine site"] += 1
            elif called:
                miss["3 route exists, model calls the wrong direction"] += 1
            elif c["reach"] == "up":
                miss["4 route exists, model reads no change"] += 1
            elif others >= 10:
                miss["5 no route; protein changes 10+ other readouts"] += 1
            else:
                miss["6 no route; other"] += 1
        if called and c["measured"] == "not":
            h, pip3 = best.get(k, (None, False))
            if others <= 3:
                wrong["1 protein changes at most 3 other readouts (not active here)"] += 1
            elif pip3:
                wrong["2 route through PIP3 (the switch)"] += 1
            elif h is not None and h > B.LONG_ROUTE:
                wrong["3 route longer than 12 steps"] += 1
            else:
                wrong["4 other"] += 1
    res_out = {"changes": sum(B.changed(c["measured"]) for c in cases),
               "caught": sum(B.changed(c["measured"]) and c["frozen"] == c["measured"] for c in cases),
               "missed": dict(sorted(miss.items())), "false_calls": dict(sorted(wrong.items())),
               "false_calls_total": sum(wrong.values())}
    Path(out).write_text(json.dumps(res_out, indent=1))
    print(json.dumps(res_out, indent=1))


if __name__ == "__main__":
    main()
