#!/usr/bin/env python3
"""specs/182: Lun model calls by the over-expressed protein's activation class. Both halves are FIT.

Usage: act_182.py PREDICTIONS BUILD OUT_JSON A_FILE B_FILE [CLASS_DIR]
"""
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
for d in ("173", "176", "177", "179", "181"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import budget_176 as B  # noqa: E402
import engage_179 as E  # noqa: E402
import lun_a_173 as A  # noqa: E402
import routes_181 as R  # noqa: E402
import score_177 as S  # noqa: E402

SEED, DRAWS, MAX_STEPS, MAX_RECALL_DROP, ALPHA = 182, 10_000, 10, 0.01, 0.05
CLASSES = "ABCN"


def read_classes(path=HERE.parent / "activation_classes.tsv"):
    return {r["gene"]: r["class"] for r in A.read_tsv(path)}


def stats(cases):
    calls = [c for c in cases if B.changed(c["model"])]
    return {"pairs": len(cases), "proteins": len({c["gene"] for c in cases}),
            "share_changed": sum(B.changed(c["measured"]) for c in cases) / len(cases) if cases else None,
            "calls": len(calls), "precision": sum(c["model"] == c["measured"] for c in calls) / len(calls) if calls else None}


def precision_diff(cases, a="A", b="C", n=DRAWS, seed=SEED):
    """Model call precision in class a minus class b, proteins resampled."""
    by = defaultdict(list)
    for c in cases:
        if c["cls"] in (a, b):
            by[c["gene"]].append(c)
    genes = sorted(by)

    def diff(sample):
        pa, pb = (stats([c for c in sample if c["cls"] == k])["precision"] for k in (a, b))
        return None if pa is None or pb is None else pa - pb
    rng = random.Random(seed)
    draws = [d for d in (diff([c for g in (rng.choice(genes) for _ in genes) for c in by[g]]) for _ in range(n))
             if d is not None]
    if not draws or diff(cases) is None:
        return None
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {"point": diff(cases), "ci": [float(lo), float(hi)]}


def table(cases):
    return {k: stats([c for c in cases if c["cls"] == k]) for k in CLASSES}


def withhold(routed, everything, cls="C"):
    moved = [c for c in routed if c["cls"] == cls and B.changed(c["model"])]
    fixed, broken = sum(c["measured"] == "not" for c in moved), sum(c["measured"] == c["model"] for c in moved)
    real = [c for c in everything if B.changed(c["measured"])]
    r0 = sum(c["model"] == c["measured"] for c in real) / len(real)
    p = S.sign_test(fixed, broken)
    return {"withheld": len(moved), "fixed": fixed, "broken": broken, "sign_test_p": p, "recall": r0,
            "recall_after": r0 - broken / len(real),
            "passes": bool(fixed > broken and p < ALPHA and broken / len(real) <= MAX_RECALL_DROP)}


def analyse(cases, lengths, cls, tk, eng):
    routed = [dict(c, cls=cls.get(c["gene"], "N"), steps=lengths.get((c["gene"], c["readout"])))
              for c in cases if c["reach"] == "up"]
    short = [c for c in routed if c["steps"] is not None and c["steps"] <= MAX_STEPS]
    no_tk = [c for c in routed if c["gene"] not in tk]
    res = {"routed": table(routed), "routed_up_to_10_steps": table(short), "routed_without_tyrosine_kinases": table(no_tk),
           "A_minus_C": precision_diff(routed), "A_minus_C_up_to_10_steps": precision_diff(short),
           "A_minus_C_without_tyrosine_kinases": precision_diff(no_tk),
           "A_minus_B": precision_diff(routed, "A", "B"),
           "proteins_by_class": dict(Counter(cls.get(g, "N") for g in {c["gene"] for c in cases})),
           "engagement_by_class": {k: {"engaged": sum(v for g, v in eng.items() if cls.get(g, "N") == k),
                                       "not_engaged": sum(not v for g, v in eng.items() if cls.get(g, "N") == k)} for k in CLASSES},
           "withhold_C": withhold(routed, cases), "withhold_B_and_C": None}
    return res


def main():
    import dataset_roles
    pred_file, build, out, a_file, b_file = sys.argv[1:6]
    tk, tyr = S.read_classes(sys.argv[6]) if len(sys.argv) > 6 else S.read_classes()
    genes, cls, direct = A.marker_genes(), read_classes(), E.direct_table()
    halves, m_all = {}, {}
    for half, path in (("A", a_file), ("B", b_file)):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        m = S.measured(path, genes, half)
        halves[half] = (m, S.build_cases(m, S.predictions(pred_file, half), tk, tyr))
        m_all.update(m)
    allc = halves["A"][1] + halves["B"][1]
    lengths = R.route_lengths(pred_file, build, {(c["gene"], c["readout"]): c for c in allc if c["reach"] == "up"})
    res = {f"half_{h}": analyse(cs, lengths, cls, tk, E.engagement(m, direct)) for h, (m, cs) in halves.items()}
    res["pooled"] = analyse(allc, lengths, cls, tk, E.engagement(m_all, direct))
    d = res["pooled"]["A_minus_C"]
    res["supported"] = bool(d and d["ci"][0] > 0 and all((res[f"half_{h}"]["A_minus_C"] or {"point": -1})["point"] > 0 for h in "AB"))
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
