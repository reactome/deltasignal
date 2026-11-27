#!/usr/bin/env python3
"""specs/177: the candidate Lun protocol (frozen model under the abundance rule, plus the tyrosine-kinase rule),
scored on one half of Lun 2019 with the specs/173 measured call.

  rule: the protein is tyrosine-specific (tyrosine_kinases.tsv) and the readout has a tyrosine at a registered site
        (readout_sites.tsv) -> UP, where the model itself calls no change

Usage: score_177.py HALF HALF_FILE PREDICTIONS OUT_JSON      HALF = A (development) or B (the sealed test, once)
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
sys.path.insert(0, str(HERE.parents[1] / "176"))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
from common_173 import split_of  # noqa: E402

SEED = 177
ROLE = {"A": ("lun2019_a", "fit"), "B": ("lun2019_b", "score")}
ARMS = ("candidate", "model", "reach", "rule_only", "never")


def read_classes(here=HERE.parent):
    tk = {r["gene"] for r in A.read_tsv(here / "tyrosine_kinases.tsv") if r["tyrosine_specific"] == "1"}
    tyr = {r["readout"] for r in A.read_tsv(here / "readout_sites.tsv") if r["has_tyrosine"] == "1"}
    return tk, tyr


def load_half(path, genes, half, timepoint="0"):
    """As specs/173's reader, for either half: controls and that half's proteins only; anything else is refused."""
    out = defaultdict(lambda: defaultdict(dict))
    for r in csv.DictReader(open(path)):
        m = r["marker"]
        if not A.is_control(m) and (m not in genes or split_of(genes[m]) != half):
            raise SystemExit(f"a sample of the other half or unmapped ({m}) reached the half-{half} reader: refusing")
        if r["timepoint"] != timepoint or r["target"] not in A.TARGETS:
            continue
        out[(r["experiment"], m, r["row_col"])][A.TARGETS[r["target"]]][int(r["bin"])] = (
            float(r["counts"]), float(r["median"]))
    return out


def measured(path, genes, half, timepoint="0"):
    """{(gene, readout): (E, delta)} by the specs/173 call."""
    sl = A.slopes(load_half(path, genes, half, timepoint))
    delta = A.thresholds(sl)
    return {k: (e, delta[k[1]]) for k, (_, e) in A.measured_calls(sl, genes, delta).items()}


def predictions(pred_file, half):
    import lun_score as L
    frozen, reach = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if r["self"] == "1" or r["status"] != "scored" or split_of(r["gene"]) != half:
            continue
        frozen[(r["gene"], r["readout"])].append(r["predicted"])
        reach[(r["gene"], r["readout"])].append(r["reach_predicted"])
    return {k: (L.combine_networks(frozen[k]), L.combine_networks(reach[k])) for k in frozen}


def apply_rule(model_call, gene, readout, tk, tyr):
    """The candidate's call: the model's own change call stands; otherwise UP for a tyrosine kinase at a tyrosine site."""
    if model_call in ("up", "down"):
        return model_call
    return "up" if gene in tk and readout in tyr else "not"


def sign_test(right, wrong):
    n = right + wrong
    return sum(comb(n, i) for i in range(right, n + 1)) / 2 ** n if n else 1.0


def build_cases(m, pred, tk, tyr, mult=1.0):
    cases = []
    for k in sorted(pred):
        if k not in m:
            continue
        e, d = m[k]
        model = pred[k][0] if pred[k][0] in ("up", "down") else "not"
        cases.append({"gene": k[0], "readout": k[1], "measured": B.call_of(e, mult * d), "E": e, "delta": d,
                      "model": model, "frozen": model, "reach": "up" if pred[k][1] == "up" else "not",
                      "candidate": apply_rule(model, k[0], k[1], tk, tyr),
                      "rule_only": "up" if k[0] in tk and k[1] in tyr else "not", "never": "not"})
    return cases


def arm_scores(cases):
    import lun_score as L
    out = {}
    for arm in ARMS:
        real = [c for c in cases if B.changed(c["measured"])]
        calls = [c for c in cases if B.changed(c[arm])]
        out[arm] = {"balanced_accuracy": L.balanced_accuracy([(c["measured"], c[arm]) for c in cases]),
                    "calls": len(calls), "recall": sum(c[arm] == c["measured"] for c in real) / len(real) if real else None,
                    "precision_exact": sum(c[arm] == c["measured"] for c in calls) / len(calls) if calls else None}
    return out


def main():
    import dataset_roles
    import lun_score as L
    half, path, pred_file, out = sys.argv[1:5]
    dataset_roles.require_role(*ROLE[half])
    genes = A.marker_genes()
    tk, tyr = read_classes()
    m = measured(path, genes, half)
    res = {"half": half, "measured_pairs": len(m), "tyrosine_specific_proteins_measured": len({g for g, _ in m} & tk)}

    # positive control, as specs/127 step 7 and specs/173: scoring stops if it fails
    calls = {k: (B.call_of(e, d), e) for k, (e, d) in m.items()}
    selfc = A.self_pairs(calls, A.read_tsv(A.MAPS / "readouts.tsv"))
    passes, ups, nmeas = L.positive_control([c for c, _ in selfc.values()])
    res["positive_control"] = {"passes": passes, "up": ups, "measured": nmeas}
    if not passes:
        Path(out).write_text(json.dumps(res, indent=1))
        raise SystemExit("positive control FAILS: scoring stops")

    pred = predictions(pred_file, half)
    cases = build_cases(m, pred, tk, tyr)
    res["pairs"] = len(cases)
    res["proteins"] = len({c["gene"] for c in cases})
    res["measured_classes"] = dict(Counter(c["measured"] for c in cases))
    res["arms"] = arm_scores(cases)
    res["differences"] = {}
    for a, b in (("candidate", "reach"), ("candidate", "model"), ("model", "reach"), ("model", "never"),
                 ("rule_only", "reach")):
        pt, lo, hi = L.bootstrap_difference(cases, a, b, seed=SEED)
        res["differences"][f"{a}-{b}"] = {"point": pt, "ci": [lo, hi], "verdict": L.verdict(lo, hi)}
    added = [c for c in cases if c["candidate"] != c["model"]]
    right, wrong = sum(c["measured"] == "up" for c in added), sum(c["measured"] != "up" for c in added)
    res["added_calls"] = {"n": len(added), "right": right, "measured_not": sum(c["measured"] == "not" for c in added),
                          "measured_down": sum(c["measured"] == "down" for c in added), "sign_test_p": sign_test(right, wrong),
                          "proteins": len({c["gene"] for c in added})}
    res["gate"] = {"primary_candidate_beats_reach": res["differences"]["candidate-reach"]["ci"][0] > 0,
                   "secondary_candidate_beats_model": res["differences"]["candidate-model"]["ci"][0] > 0,
                   "secondary_added_calls": right > wrong and res["added_calls"]["sign_test_p"] < 0.05}

    # specs/173 R2 for the model, as registered there
    res["R2_model"] = {"enrichment": A.enrichment(cases), "enrichment_ci": A.protein_bootstrap(cases, A.enrichment),
                       "lift": A.lift(cases), "lift_ci": A.protein_bootstrap(cases, A.lift)}

    # stricter measured calls (specs/176 T1b)
    res["by_truth_threshold"] = {}
    for mult in (2, 3):
        cs = build_cases(m, pred, tk, tyr, mult=mult)
        s = arm_scores(cs)
        pt, lo, hi = L.bootstrap_difference(cs, "candidate", "reach", seed=SEED)
        res["by_truth_threshold"][f"{mult}delta"] = {
            "changes": sum(B.changed(c["measured"]) for c in cs),
            **{arm: round(s[arm]["balanced_accuracy"], 4) for arm in ARMS}, "candidate-reach": {"point": pt, "ci": [lo, hi]}}

    # the class table on every measured pair of this half (specs/176's axis)
    tab = defaultdict(Counter)
    for (g, r), (e, d) in m.items():
        tab[("tyrosine kinase" if g in tk else "other protein", "tyrosine site" if r in tyr else "S/T site")][B.call_of(e, d)] += 1
    res["class_table"] = {" x ".join(k): {"pairs": sum(c.values()), "up": round(c["up"] / sum(c.values()), 3),
                                           "down": round(c["down"] / sum(c.values()), 3)} for k, c in sorted(tab.items())}
    per = defaultdict(lambda: [0, 0])
    for (g, r), (e, d) in m.items():
        if g in tk and r in tyr:
            per[g][0] += B.call_of(e, d) == "up"
            per[g][1] += 1
    res["tyrosine_kinases_raising_half_or_more"] = f"{sum(u >= n / 2 for u, n in per.values())}/{len(per)}"

    # specs/176 T2 and T1a on this half
    for c in cases:
        c["fold"] = 1.0
    t2 = B.t2(m, cases)
    res["auroc"] = {**t2["auroc"], "candidate_call": B.auroc([B.changed(c["candidate"]) for c in cases],
                                                            [B.changed(c["measured"]) for c in cases])}
    res["auroc"].pop("model_fold_magnitude", None)
    m10 = measured(path, genes, half, timepoint="10")
    rep = [B.call_of(*m10[k]) == B.call_of(e, d) for k, (e, d) in m.items() if B.changed(B.call_of(e, d)) and k in m10]
    res["repeat_in_egf_well"] = {"basal_changes": len(rep), "share": sum(rep) / len(rep) if rep else None}

    Path(out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
