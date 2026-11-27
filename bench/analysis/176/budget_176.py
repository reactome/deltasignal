#!/usr/bin/env python3
"""specs/176: the Lun 2019 half-A failure budget (tests T1, T2, T4, T5, T6; T3 is routes_176.py).

Usage: budget_176.py A_FILE PREDICTIONS TRACES_ALL OUT_JSON
  TRACES_ALL  specs/174's traces_all.tsv (shortest route per reached pair), for T5
Half B is never read: the reader is specs/173's, which refuses any half-B protein.
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
import lun_a_173 as A  # noqa: E402
from common_173 import split_of  # noqa: E402

STRATA = ((1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, float("inf")))
LONG_ROUTE = 12
PIP3 = "57836"            # ChEBI id of PI(3,4,5)P3, as traces label it
AKT_FEEDBACK_UP = ("p-S6", "p-p70S6K", "p-mTOR")


def changed(c):
    return c in ("up", "down")


def call_of(e, d):
    return "up" if e > d else "down" if e < -d else "not"


def measured(a_file, genes, timepoint="0", low=None, high=None):
    """{(gene, readout): (E, delta)} by the specs/173 call; bins may be overridden (T1c)."""
    keep = A.LOW_BINS, A.HIGH_BINS
    if low is not None:
        A.LOW_BINS, A.HIGH_BINS = low, high
    try:
        sl = A.slopes(A.load_a(a_file, genes, timepoint=timepoint))
    finally:
        A.LOW_BINS, A.HIGH_BINS = keep
    delta = A.thresholds(sl)
    return {k: (e, delta[k[1]]) for k, (_, e) in A.measured_calls(sl, genes, delta).items()}


def predictions(pred_file):
    import lun_score as L
    frozen, reach, folds = defaultdict(list), defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if r["self"] == "1" or r["status"] != "scored" or split_of(r["gene"]) != "A":
            continue
        k = (r["gene"], r["readout"])
        frozen[k].append(r["predicted"])
        reach[k].append(r["reach_predicted"])
        folds[k].append(float(r["pred_ui"]))
    return {k: {"frozen": L.as_scored(L.combine_networks(frozen[k])), "reach": L.combine_networks(reach[k]),
                "fold": max(folds[k], key=lambda f: abs(np.log(max(f, 1e-6))))} for k in frozen}


def auroc(scores, labels):
    """Probability a changed pair outranks an unchanged one (ties count half)."""
    s, y = np.asarray(scores, float), np.asarray(labels, bool)
    if y.all() or not y.any():
        return None
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s))
    ranks[order] = np.arange(1, len(s) + 1)
    for v in np.unique(s):                       # average ranks over ties
        m = s == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    n1, n0 = y.sum(), (~y).sum()
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def scores(cases, col="frozen", truth="measured"):
    """Recall (right direction), precision of change calls, three-class balanced accuracy, enrichment."""
    import lun_score as L
    real = [c for c in cases if changed(c[truth])]
    calls = [c for c in cases if changed(c[col])]
    on = [c for c in cases if c["reach"] == "up"]
    off = [c for c in cases if c["reach"] != "up"]
    r_on = sum(changed(c[truth]) for c in on) / len(on) if on else None
    r_off = sum(changed(c[truth]) for c in off) / len(off) if off else None
    return {"pairs": len(cases), "changes": len(real), "calls": len(calls),
            "recall": sum(c[col] == c[truth] for c in real) / len(real) if real else None,
            "precision_exact": sum(c[col] == c[truth] for c in calls) / len(calls) if calls else None,
            "precision_any": sum(changed(c[truth]) for c in calls) / len(calls) if calls else None,
            "balanced_accuracy": L.balanced_accuracy([(c[truth], c[col]) for c in cases]),
            "rate_route": r_on, "rate_no_route": r_off,
            "enrichment": r_on / r_off if r_on is not None and r_off else None}


def t1(a_file, genes, m0, cases):
    out = {}
    # a. the EGF-10-minute well as a second measurement
    m10 = measured(a_file, genes, timepoint="10")
    rep = defaultdict(lambda: [0, 0])
    strat = {s: [0, 0] for s in STRATA}
    for k, (e, d) in m0.items():
        c0 = call_of(e, d)
        if not changed(c0) or k not in m10:
            continue
        same = call_of(*m10[k]) == c0
        rep[k[1]][0] += same
        rep[k[1]][1] += 1
        for lo, hi in STRATA:
            if lo <= abs(e) / d < hi:
                strat[(lo, hi)][0] += same
                strat[(lo, hi)][1] += 1
    tot = [sum(v[0] for v in rep.values()), sum(v[1] for v in rep.values())]
    unchanged = [k for k, (e, d) in m0.items() if not changed(call_of(e, d)) and k in m10]
    out["a_repeat"] = {"basal_changes": tot[1], "repeat_same_direction": tot[0], "share": tot[0] / tot[1],
                       "unchanged_called_changed_at_10": sum(changed(call_of(*m10[k])) for k in unchanged) / len(unchanged),
                       "by_effect_size": {f"{lo}-{hi}": {"n": n, "share": s / n if n else None}
                                          for (lo, hi), (s, n) in strat.items()},
                       "by_readout": {r: {"n": n, "share": round(s / n, 3)} for r, (s, n) in sorted(rep.items())}}
    # b. effect-size strata on the scored pairs, and the scores at stricter calls
    by = {}
    for lo, hi in STRATA:
        sub = [c for c in cases if changed(c["measured"]) and lo <= abs(c["E"]) / c["delta"] < hi]
        on = [c for c in sub if c["reach"] == "up"]
        called = [c for c in on if changed(c["frozen"])]
        by[f"{lo}-{hi}"] = {"changes": len(sub), "with_route": len(on), "share_with_route": len(on) / len(sub) if sub else None,
                            "model_called": len(called), "model_right_direction": sum(c["frozen"] == c["measured"] for c in called)}
    out["b_strata"] = by
    out["b_route_share_of_unchanged"] = (sum(c["reach"] == "up" for c in cases if not changed(c["measured"]))
                                         / sum(not changed(c["measured"]) for c in cases))
    out["b_scores"] = {}
    for mult in (1, 2, 3):
        for c in cases:
            c["truth_m"] = call_of(c["E"], mult * c["delta"])
        out["b_scores"][f"{mult}delta"] = {"model": scores(cases, "frozen", "truth_m"),
                                           "reach_ba": scores(cases, "reach_call", "truth_m")["balanced_accuracy"]}
    # c. bin choice
    alt = measured(a_file, genes, low=(1, 2, 3), high=(6, 7, 8))
    both = [k for k in m0 if k in alt]
    c0 = {k: call_of(*m0[k]) for k in both}
    c1 = {k: call_of(*alt[k]) for k in both}
    ch = [k for k in both if changed(c0[k])]
    out["c_bins"] = {"pairs": len(both), "agree": sum(c0[k] == c1[k] for k in both) / len(both),
                     "changes_kept": sum(c0[k] == c1[k] for k in ch) / len(ch)}
    return out


def t2(m0, cases):
    """Leave-one-out protein rate x readout rate, against route and the model's call."""
    call = {k: changed(call_of(e, d)) for k, (e, d) in m0.items()}
    by_g, by_r = defaultdict(list), defaultdict(list)
    for (g, r), v in call.items():
        by_g[g].append(v)
        by_r[r].append(v)
    two, prot, ro, y, route, model, mag = [], [], [], [], [], [], []
    for c in cases:
        v = call[(c["gene"], c["readout"])]
        pg = (sum(by_g[c["gene"]]) - v) / max(1, len(by_g[c["gene"]]) - 1)
        pr = (sum(by_r[c["readout"]]) - v) / max(1, len(by_r[c["readout"]]) - 1)
        c["protein_rate"], c["readout_rate"] = pg, pr
        two.append(pg * pr)
        prot.append(pg)
        ro.append(pr)
        y.append(v)
        route.append(c["reach"] == "up")
        model.append(changed(c["frozen"]))
        mag.append(abs(np.log(max(c["fold"], 1e-6))))
    out = {"auroc": {"protein_x_readout": auroc(two, y), "protein_rate": auroc(prot, y), "readout_rate": auroc(ro, y),
                     "route": auroc(route, y), "model_call": auroc(model, y), "model_fold_magnitude": auroc(mag, y)}}
    # does a route add anything once the two factors are known? change rate by two-factor quartile and route
    q = np.quantile(two, [0.25, 0.5, 0.75])
    table = defaultdict(lambda: [0, 0])
    for s, v, rt in zip(two, y, route):
        key = (int(np.searchsorted(q, s, side="right")), bool(rt))
        table[key][0] += v
        table[key][1] += 1
    out["by_quartile_and_route"] = {f"q{qi + 1}_{'route' if rt else 'none'}": {"n": n, "rate": round(s / n, 3)}
                                    for (qi, rt), (s, n) in sorted(table.items())}
    return out


def t4(m0):
    calls = defaultdict(dict)
    for (g, r), (e, d) in m0.items():
        calls[g][r] = call_of(e, d)
    tab = Counter()
    for g, cs in calls.items():
        if "p-AKT" not in cs:
            continue
        tab[(cs["p-AKT"] == "down", any(cs.get(r) == "up" for r in AKT_FEEDBACK_UP))] += 1
    from scipy.stats import fisher_exact
    a, b, c, d = tab[(True, True)], tab[(True, False)], tab[(False, True)], tab[(False, False)]
    odds, p = fisher_exact([[a, b], [c, d]])
    up = Counter()
    for g, cs in calls.items():
        up[(cs.get("p-AKT"), tuple(cs.get(r) for r in AKT_FEEDBACK_UP))] += 1
    return {"akt_down_feedback_up": a, "akt_down_no_feedback": b, "akt_not_down_feedback_up": c,
            "akt_not_down_no_feedback": d, "odds_ratio": odds, "fisher_p": p,
            "akt_calls": dict(Counter(cs.get("p-AKT") for cs in calls.values())),
            "akt_down_other_readouts_changed_median": float(np.median(
                [sum(changed(v) for r, v in cs.items() if r != "p-AKT") for cs in calls.values() if cs.get("p-AKT") == "down"])),
            "akt_not_down_other_readouts_changed_median": float(np.median(
                [sum(changed(v) for r, v in cs.items() if r != "p-AKT") for cs in calls.values() if cs.get("p-AKT") != "down"]))}


def t5(cases, traces):
    best = {}
    for r in csv.DictReader(open(traces), delimiter="\t"):
        if r["hops"] == "":
            continue
        k = (r["gene"], r["readout"])
        if k not in best or int(r["hops"]) < best[k][0]:
            best[k] = (int(r["hops"]), PIP3 in r["path"])
    out = {"reached": sum(c["reach"] == "up" for c in cases),
           "reached_traced": sum(c["reach"] == "up" and (c["gene"], c["readout"]) in best for c in cases)}
    rules = {"a_through_pip3": lambda h, p, f: p, "b_longer_than_12": lambda h, p, f: h > LONG_ROUTE,
             "c_either": lambda h, p, f: p or h > LONG_ROUTE,
             "d_fold_railed_posthoc": lambda h, p, f: f >= 40 or f <= 0.025}
    base = scores(cases)
    out["baseline"] = base
    for name, rule in rules.items():
        fixed = broken = held = 0
        for c in cases:
            c["alt"] = c["frozen"]
            k = (c["gene"], c["readout"])
            h, p = best.get(k, (0, False))
            if changed(c["frozen"]) and (k in best or name.startswith("d")) and rule(h, p, c["fold"]):
                c["alt"] = "not"
                held += 1
                fixed += c["measured"] == "not"
                broken += c["measured"] == c["frozen"]
        s = scores(cases, "alt")
        out[name] = {"calls_withheld": held, "fixed": fixed, "broken": broken, "recall": s["recall"],
                     "precision_exact": s["precision_exact"], "balanced_accuracy": s["balanced_accuracy"]}
    return out


def t6(m0, cases):
    call = {k: changed(call_of(e, d)) for k, (e, d) in m0.items()}
    n_changed = Counter()
    for (g, r), v in call.items():
        n_changed[g] += v
    bands = (("0", 0, 0), ("1-3", 1, 3), ("4-9", 4, 9), ("10+", 10, 99))
    out = {}
    for name, lo, hi in bands:
        sub = [c for c in cases if lo <= n_changed[c["gene"]] - call[(c["gene"], c["readout"])] <= hi]
        if sub:
            s = scores(sub)
            out[name] = {"proteins": len({c["gene"] for c in sub}),
                         **{k: s[k] for k in ("pairs", "changes", "calls", "recall", "precision_exact", "rate_route",
                                              "rate_no_route", "enrichment")}}
    dist = Counter(min(v, 10) for v in (n_changed[g] for g in {g for g, _ in call}))
    out["proteins_by_readouts_changed"] = {str(k): dist[k] for k in sorted(dist)}
    return out


def main():
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, pred_file, traces, out = sys.argv[1:5]
    genes = A.marker_genes()
    m0 = measured(a_file, genes)
    pred = predictions(pred_file)
    cases = [{"gene": k[0], "readout": k[1], "E": m0[k][0], "delta": m0[k][1], "measured": call_of(*m0[k]),
              "reach_call": "up" if pred[k]["reach"] == "up" else "not", **pred[k]}
             for k in sorted(pred) if k in m0]
    res = {"pairs": len(cases), "measured_pairs_all": len(m0), "baseline": scores(cases),
           "T1": t1(a_file, genes, m0, cases), "T2": t2(m0, cases), "T4": t4(m0), "T5": t5(cases, traces),
           "T6": t6(m0, cases)}
    Path(out).write_text(json.dumps(res, indent=1, default=str))
    with open(Path(out).with_name("cases_176.tsv"), "w", newline="") as f:
        cols = ["gene", "readout", "measured", "E", "delta", "frozen", "reach", "fold", "protein_rate", "readout_rate"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(cases)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
