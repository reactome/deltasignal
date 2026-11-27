#!/usr/bin/env python3
"""specs/179: first-hop engagement on Lun 2019. Tests E1-E3 on one half or both pooled.

Usage: engage_179.py PREDICTIONS OUT_JSON A_FILE [B_FILE]
"""
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
sys.path.insert(0, str(HERE.parents[1] / "176"))
sys.path.insert(0, str(HERE.parents[1] / "177"))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
import score_177 as S  # noqa: E402

SEED, DRAWS = 179, 10_000
MAX_RECALL_DROP, ALPHA = 0.01, 0.05


def direct_table(path=HERE.parent / "direct_substrates.tsv"):
    return {(r["gene"], r["readout"]) for r in A.read_tsv(path)}


def engagement(m, direct):
    """{gene: True/False} for proteins with a measured direct pair: any of them reads UP."""
    by = defaultdict(list)
    for k in direct:
        if k in m:
            by[k[0]].append(B.call_of(*m[k]) == "up")
    return {g: any(v) for g, v in by.items()}


def placebo(m, reach_of, direct):
    """{gene: True/False}: the protein raises at least one measured readout it has no route to (direct pairs left out)."""
    by = defaultdict(list)
    for k, (e, d) in m.items():
        if k in direct or reach_of.get(k) == "up" or k not in reach_of:
            continue
        by[k[0]].append(B.call_of(e, d) == "up")
    return {g: any(v) for g, v in by.items()}


def precision(cases, col):
    calls = [c for c in cases if B.changed(c[col])]
    return sum(c[col] == c["measured"] for c in calls) / len(calls) if calls else None


def group_stats(cases, col):
    s = B.scores(cases, col)
    return {k: s[k] for k in ("pairs", "changes", "calls", "recall", "precision_exact", "rate_route", "rate_no_route",
                              "enrichment")}


def boot_diff(cases, key, col, n=DRAWS, seed=SEED):
    """Precision of `col` change calls where c[key] is True minus where it is False; proteins resampled."""
    by = defaultdict(list)
    for c in cases:
        by[c["gene"]].append(c)
    genes = sorted(by)

    def diff(sample):
        a, b = precision([c for c in sample if c[key]], col), precision([c for c in sample if not c[key]], col)
        return None if a is None or b is None else a - b
    rng = random.Random(seed)
    draws = [d for d in (diff([c for g in (rng.choice(genes) for _ in genes) for c in by[g]]) for _ in range(n))
             if d is not None]
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {"point": diff(cases), "ci": [float(lo), float(hi)]}


def gate(cases, col):
    """E3: change calls withheld for not-engaged proteins."""
    import lun_score as L
    for c in cases:
        c["gated"] = c[col] if c["engaged"] else "not"
        c["reach_gated"] = c["reach"] if c["engaged"] else "not"
    moved = [c for c in cases if c["gated"] != c[col]]
    fixed = sum(c["measured"] == "not" for c in moved)
    broken = sum(c["measured"] == c[col] for c in moved)
    p = S.sign_test(fixed, broken)
    real = [c for c in cases if B.changed(c["measured"])]
    r0 = sum(c[col] == c["measured"] for c in real) / len(real) if real else None
    r1 = sum(c["gated"] == c["measured"] for c in real) / len(real) if real else None
    ba = {k: L.balanced_accuracy([(c["measured"], c[k]) for c in cases]) for k in (col, "gated", "reach", "reach_gated")}
    pt, lo, hi = L.bootstrap_difference(cases, "gated", col, seed=SEED)
    return {"withheld": len(moved), "fixed": fixed, "broken": broken, "sign_test_p": p, "recall": r0, "recall_gated": r1,
            "balanced_accuracy": ba, "gated-model": {"point": pt, "ci": [lo, hi]},
            "passes": bool(fixed > broken and p < ALPHA and r0 is not None and r0 - r1 <= MAX_RECALL_DROP)}


def run(label, m, pred, direct, tk, tyr):
    cases = S.build_cases(m, pred, tk, tyr)
    reach_of = {(c["gene"], c["readout"]): c["reach"] for c in cases}
    eng, plc = engagement(m, direct), placebo(m, reach_of, direct)
    test = [dict(c, engaged=eng[c["gene"]], placebo=plc.get(c["gene"], False))
            for c in cases if c["gene"] in eng and (c["gene"], c["readout"]) not in direct]
    res = {"proteins_with_measured_direct_pair": len(eng), "engaged": sum(eng.values()),
           "not_engaged": sum(not v for v in eng.values()), "test_pairs": len(test),
           "test_proteins": len({c["gene"] for c in test})}
    for col in ("model", "candidate"):
        r = {"engaged": group_stats([c for c in test if c["engaged"]], col),
             "not_engaged": group_stats([c for c in test if not c["engaged"]], col),
             "E1_precision_difference": boot_diff(test, "engaged", col),
             "E2_placebo": {"yes": group_stats([c for c in test if c["placebo"]], col),
                            "no": group_stats([c for c in test if not c["placebo"]], col),
                            "precision_difference": boot_diff(test, "placebo", col)},
             "E2_within_placebo": {
                 f"placebo_{'yes' if pv else 'no'}": {
                     "engaged": group_stats([c for c in test if c["placebo"] == pv and c["engaged"]], col),
                     "not_engaged": group_stats([c for c in test if c["placebo"] == pv and not c["engaged"]], col)}
                 for pv in (True, False)},
             "E3_gate": gate(test, col)}
        res[col] = r
    no_tk = [c for c in test if c["gene"] not in tk]
    res["without_tyrosine_kinases"] = {
        "test_pairs": len(no_tk), "proteins": len({c["gene"] for c in no_tk}),
        "engaged": group_stats([c for c in no_tk if c["engaged"]], "model"),
        "not_engaged": group_stats([c for c in no_tk if not c["engaged"]], "model"),
        "E1_precision_difference": boot_diff(no_tk, "engaged", "model") if no_tk else None,
        "E3_gate": gate(no_tk, "model") if no_tk else None}
    res["label"] = label
    return res


def main():
    import dataset_roles
    pred_file, out, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    genes = A.marker_genes()
    tk, tyr = S.read_classes()
    direct = direct_table()
    res, m_all, p_all = {}, {}, {}
    for half, path in zip("AB", files):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        m, pred = S.measured(path, genes, half), S.predictions(pred_file, half)
        res[f"half_{half}"] = run(half, m, pred, direct, tk, tyr)
        m_all.update(m)
        p_all.update(pred)
    if len(files) == 2:
        res["pooled"] = run("pooled", m_all, p_all, direct, tk, tyr)
    Path(out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
