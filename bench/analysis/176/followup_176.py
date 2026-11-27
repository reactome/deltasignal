#!/usr/bin/env python3
"""specs/176 follow-ups (post hoc, labelled as such in research.md), half A only.

  F1  repeat rate of basal change calls in the EGF well, split by whether EGF itself moves the readout in control wells
  F2  curated enzyme class from the networks: does the protein catalyse its own phosphorylation (self-activating),
      is it a tyrosine kinase, is it a catalyst at all; change rates and model scores by class
  F4  does a route add to protein rate x readout rate? Logistic regression, protein-grouped cross-validation

Usage: followup_176.py A_FILE CASES_176 BUILD OUT_JSON
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402

EGF_SHIFT_SD = 2.0


def control_levels(a_file, genes, timepoint):
    """{readout: [count-weighted mean of bin medians, per FLAG-GFP well]}."""
    out = defaultdict(list)
    for s, per in A.load_a(a_file, genes, timepoint=timepoint).items():
        if not s[1].startswith("FLAG-GFP"):
            continue
        for r, bins in per.items():
            n = sum(c for c, _ in bins.values())
            if n:
                out[r].append(sum(c * m for c, m in bins.values()) / n)
    return out


def f1(a_file, genes):
    l0, l10 = control_levels(a_file, genes, "0"), control_levels(a_file, genes, "10")
    egf = {r: abs(np.median(l10[r]) - np.median(l0[r])) > EGF_SHIFT_SD * np.std(l0[r]) for r in l0 if r in l10}
    m0, m10 = B.measured(a_file, genes), B.measured(a_file, genes, timepoint="10")
    tab = {True: [0, 0], False: [0, 0]}
    strata = {(e, s): [0, 0] for e in (True, False) for s in B.STRATA}
    for k, (e, d) in m0.items():
        c0 = B.call_of(e, d)
        if not B.changed(c0) or k not in m10 or k[1] not in egf:
            continue
        same = B.call_of(*m10[k]) == c0
        tab[egf[k[1]]][0] += same
        tab[egf[k[1]]][1] += 1
        for lo, hi in B.STRATA:
            if lo <= abs(e) / d < hi:
                strata[(egf[k[1]], (lo, hi))][0] += same
                strata[(egf[k[1]], (lo, hi))][1] += 1
    return {"egf_moves_readout": sorted(r for r, v in egf.items() if v),
            "repeat_egf_readouts": {"n": tab[True][1], "share": tab[True][0] / max(1, tab[True][1])},
            "repeat_other_readouts": {"n": tab[False][1], "share": tab[False][0] / max(1, tab[False][1])},
            "other_readouts_by_effect_size": {f"{lo}-{hi}": {"n": n, "share": round(s / n, 3) if n else None}
                                              for (e, (lo, hi)), (s, n) in strata.items() if not e}}


def phospho(mods):
    return {m.strip() for m in mods.split(";") if "phospho" in m}


def enzyme_classes(build, accessions):
    """{accession: {"catalyst", "self_phos", "tyr_kinase", "st_kinase"}} from every network of the build.
    self_phos: the protein is in a catalyst and in an input of a reaction whose output carries it with a phospho
    residue no input copy has. tyr/st_kinase: catalyst of a reaction whose output gains such a residue on any protein."""
    want = set(accessions)
    cls = defaultdict(set)
    for pdir in sorted(p for p in Path(build).iterdir() if (p / "logic_network.csv").exists()):
        ref, mods = {}, {}
        for r in csv.DictReader(open(pdir / "leaf_refs.csv")):
            if r["reference_stid"].startswith("uniprot:"):
                ref[r["leaf_stid"]] = r["reference_stid"].split(":")[1].split("-")[0]
                mods[r["leaf_stid"]] = phospho(r["modifications"])
        leaves = {r["uuid"]: [x for x in r["member_leaves"].split("|") if x in ref]
                  for r in csv.DictReader(open(pdir / "nodes.csv"))}
        ins, cats, outs = defaultdict(set), defaultdict(set), defaultdict(set)
        for e in csv.DictReader(open(pdir / "logic_network.csv")):
            if e["edge_type"] == "input":
                ins[e["target_id"]].add(e["source_id"])
            elif e["edge_type"] == "catalyst":
                cats[e["target_id"]].add(e["source_id"])
            elif e["edge_type"] == "output":
                outs[e["source_id"]].add(e["target_id"])
        for rx, cs in cats.items():
            cat_acc = {ref[x] for u in cs for x in leaves.get(u, ())} & want
            if not cat_acc:
                continue
            in_mods = defaultdict(set)
            for u in ins[rx] | cs:                       # a catalyst that is also the substrate may be listed only once
                for x in leaves.get(u, ()):
                    in_mods[ref[x]] |= mods[x]
            gained = defaultdict(set)
            for u in outs[rx]:
                for x in leaves.get(u, ()):
                    new = mods[x] - in_mods.get(ref[x], set())
                    if new and ref[x] in in_mods:
                        gained[ref[x]] |= new
            for a in cat_acc:
                cls[a].add("catalyst")
                if gained.get(a):
                    cls[a].add("self_phos")
                allnew = set().union(*gained.values()) if gained else set()
                if any("tyrosine" in m for m in allnew):
                    cls[a].add("tyr_kinase")
                if any(re.search("serine|threonine", m) for m in allnew):
                    cls[a].add("st_kinase")
    return cls


def f2(a_file, genes, cases, build):
    pois = A.read_tsv(A.MAPS / "pois.tsv")
    acc = {r["gene"]: r["accession"] for r in pois}
    m0 = B.measured(a_file, genes)
    rate = defaultdict(list)
    for (g, _), (e, d) in m0.items():
        rate[g].append(B.changed(B.call_of(e, d)))
    cls = enzyme_classes(build, {acc[g] for g in rate if acc.get(g)})
    out = {"proteins_measured": len(rate)}
    groups = {"self_phos": lambda g: "self_phos" in cls.get(acc.get(g), ()),
              "tyr_kinase": lambda g: "tyr_kinase" in cls.get(acc.get(g), ()),
              "catalyst_not_self_phos": lambda g: "catalyst" in cls.get(acc.get(g), ()) and "self_phos" not in cls.get(acc.get(g), ()),
              "not_a_catalyst": lambda g: "catalyst" not in cls.get(acc.get(g), ())}
    for name, f in groups.items():
        gs = [g for g in rate if f(g)]
        rs = [np.mean(rate[g]) for g in gs]
        sub = [c for c in cases if f(c["gene"])]
        s = B.scores(sub) if sub else {}
        out[name] = {"proteins": len(gs), "median_share_of_readouts_changed": float(np.median(rs)) if rs else None,
                     "mean_share": float(np.mean(rs)) if rs else None,
                     "proteins_changing_10_plus": sum(sum(rate[g]) >= 10 for g in gs),
                     "scored": {k: s.get(k) for k in ("pairs", "changes", "calls", "recall", "precision_exact",
                                                      "rate_route", "rate_no_route", "enrichment", "balanced_accuracy")}}
    out["self_phos_genes"] = sorted(g for g in rate if groups["self_phos"](g))
    # what withholding calls for proteins that do not self-activate would do (sizing only)
    for rule, f in (("calls_only_self_phos", groups["self_phos"]),
                    ("calls_only_self_phos_or_tyr", lambda g: groups["self_phos"](g) or groups["tyr_kinase"](g))):
        fixed = broken = 0
        for c in cases:
            c["alt"] = c["frozen"]
            if B.changed(c["frozen"]) and not f(c["gene"]):
                c["alt"] = "not"
                fixed += c["measured"] == "not"
                broken += c["measured"] == c["frozen"]
        s = B.scores(cases, "alt")
        out[rule] = {"fixed": fixed, "broken": broken, **{k: s[k] for k in ("calls", "recall", "precision_exact",
                                                                             "balanced_accuracy")}}
    out["auroc_protein_level"] = {
        "self_phos_predicts_10_plus": B.auroc([groups["self_phos"](g) for g in rate], [sum(rate[g]) >= 10 for g in rate]),
        "tyr_kinase_predicts_10_plus": B.auroc([groups["tyr_kinase"](g) for g in rate], [sum(rate[g]) >= 10 for g in rate])}
    y = [B.changed(c["measured"]) for c in cases]
    out["auroc_pair_level"] = {"self_phos": B.auroc([groups["self_phos"](c["gene"]) for c in cases], y),
                               "self_phos_and_route": B.auroc([groups["self_phos"](c["gene"]) + (c["reach"] == "up") for c in cases], y)}
    return out


def f4(cases):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    y = np.array([B.changed(c["measured"]) for c in cases])
    g = np.array([c["gene"] for c in cases])
    base = np.array([[float(c["protein_rate"]), float(c["readout_rate"])] for c in cases])
    route = np.array([[c["reach"] == "up"] for c in cases], float)
    model = np.array([[B.changed(c["frozen"])] for c in cases], float)
    out = {}
    for name, x in (("protein+readout", base), ("protein+readout+route", np.hstack([base, route])),
                    ("protein+readout+model_call", np.hstack([base, model])), ("route", route)):
        pred = np.zeros(len(y))
        for tr, te in GroupKFold(n_splits=5).split(x, y, g):
            pred[te] = LogisticRegression(max_iter=1000).fit(x[tr], y[tr]).predict_proba(x[te])[:, 1]
        out[name] = B.auroc(pred, y)
    return out


def main():
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, cases_f, build, out = sys.argv[1:5]
    genes = A.marker_genes()
    cases = list(csv.DictReader(open(cases_f), delimiter="\t"))
    res = {"F1": f1(a_file, genes), "F2": f2(a_file, genes, cases, build), "F4": f4(cases)}
    Path(out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
