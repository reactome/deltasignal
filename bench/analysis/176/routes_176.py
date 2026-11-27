#!/usr/bin/env python3
"""specs/176 T3 and two post hoc structure checks, half A only.

  T3   measured changes with no route (N1): does any form of the protein reach the readout in a network that scored
       the pair, or in any network of the build?
  S1   route beyond protein rate x readout rate, with the measured call at 1, 2 and 3 delta (post hoc)
  S2   how many independent directions the response matrix has (singular values of E / delta, proteins x readouts)

Usage: routes_176.py A_FILE PREDICTIONS BUILD OUT_JSON
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
sys.path.insert(0, str(HERE.parents[1] / "174"))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
from common_173 import split_of  # noqa: E402
from trace_174 import READOUT_ACC  # noqa: E402


class Net:
    """A network's forward edges and its nodes by UniProt accession (isoform suffix dropped) and phosphosite."""

    def __init__(self, pdir):
        self.fwd = defaultdict(list)
        self.indeg = Counter()
        for e in csv.DictReader(open(pdir / "logic_network.csv")):
            self.fwd[e["source_id"]].append(e["target_id"])
            if e["edge_type"] != "depletion":
                self.indeg[e["target_id"]] += 1
        ref, mods = {}, {}
        for r in csv.DictReader(open(pdir / "leaf_refs.csv")):
            if r["reference_stid"].startswith("uniprot:"):
                ref[r["leaf_stid"]] = r["reference_stid"].split(":")[1].split("-")[0]
                mods[r["leaf_stid"]] = r["modifications"]
        self.by_acc, self.by_site, self.unmodified = defaultdict(set), defaultdict(set), defaultdict(set)
        for r in csv.DictReader(open(pdir / "nodes.csv")):
            if r["node_kind"] == "reaction":
                continue
            for leaf in r["member_leaves"].split("|"):
                if leaf not in ref:
                    continue
                self.by_acc[ref[leaf]].add(r["uuid"])
                if not mods[leaf]:
                    self.unmodified[ref[leaf]].add(r["uuid"])
                for pos in re.findall(r"phospho[^;]* at (\d+)", mods[leaf]):
                    self.by_site[(ref[leaf], int(pos))].add(r["uuid"])

    def reach(self, sources):
        seen, todo = set(sources), list(sources)
        while todo:
            for v in self.fwd.get(todo.pop(), ()):
                if v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    def readout_nodes(self, row):
        if row["kind"] == "total":
            return set(self.unmodified.get(READOUT_ACC.get(row["sites"]), ()))
        out = set()
        for s in row["sites"].split(","):
            g, pos = s.split(":")
            out |= self.by_site.get((READOUT_ACC.get(g), int(pos)), set())
        return out


def t3(cases, pred_file, build):
    n1 = {(c["gene"], c["readout"]) for c in cases if c["reach"] != "up" and B.changed(c["measured"])}
    scored = defaultdict(set)
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if (r["gene"], r["readout"]) in n1 and r["status"] == "scored":
            scored[(r["gene"], r["readout"])].add(r["pathway_id"])
    acc = {r["gene"]: r["accession"] for r in A.read_tsv(A.MAPS / "pois.tsv")}
    ro = {r["readout"]: r for r in A.read_tsv(A.MAPS / "readouts.tsv")}
    genes = {g for g, _ in n1}
    state = {k: {"matched": False, "any_form_scored": False, "any_form_anywhere": False, "root_anywhere": False}
             for k in n1}
    for pdir in sorted(p for p in Path(build).iterdir() if (p / "logic_network.csv").exists()):
        net = Net(pdir)
        targets = {r: net.readout_nodes(row) for r, row in ro.items() if row["kind"] in ("phospho", "total")}
        for g in genes:
            forms = net.by_acc.get(acc[g], set())
            mine = [k for k in n1 if k[0] == g and targets.get(k[1])]
            for k in mine:
                state[k]["matched"] = True
            if not forms or not mine:
                continue
            seen = net.reach(forms)
            roots = [u for u in forms if net.indeg[u] == 0]
            seen_root = net.reach(roots) if roots else set()
            for k in mine:
                t = targets[k[1]] - forms
                if t & seen:
                    state[k]["any_form_anywhere"] = True
                    if pdir.name in scored[k]:
                        state[k]["any_form_scored"] = True
                if t & seen_root:
                    state[k]["root_anywhere"] = True

    def label(s):
        if not s["matched"]:
            return "readout site not matched in any network"
        if s["any_form_scored"]:
            return "a non-pinned form reaches it in a scored network"
        if s["any_form_anywhere"]:
            return "a form reaches it only in another network"
        return "no form reaches it in any network"
    out = Counter(label(s) for s in state.values())
    by_readout = defaultdict(Counter)
    for k, s in state.items():
        by_readout[k[1]][label(s)] += 1
    return {"n1": len(n1), "classes": dict(out), "root_reaches_in_some_network": sum(s["root_anywhere"] for s in state.values()),
            "by_readout": {r: dict(v) for r, v in sorted(by_readout.items())},
            "pairs_with_a_form_route": sorted(f"{g}|{r}" for (g, r), s in state.items() if s["any_form_anywhere"])}


def s1(m0, cases):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    out = {}
    g = np.array([c["gene"] for c in cases])
    for mult in (1, 2, 3):
        call = {k: B.changed(B.call_of(e, mult * d)) for k, (e, d) in m0.items()}
        by_g, by_r = defaultdict(list), defaultdict(list)
        for (gg, r), v in call.items():
            by_g[gg].append(v)
            by_r[r].append(v)
        y = np.array([call[(c["gene"], c["readout"])] for c in cases])
        base = np.array([[(sum(by_g[c["gene"]]) - v) / (len(by_g[c["gene"]]) - 1),
                          (sum(by_r[c["readout"]]) - v) / (len(by_r[c["readout"]]) - 1)] for c, v in zip(cases, y)])
        route = np.array([[c["reach"] == "up"] for c in cases], float)
        res = {"changes": int(y.sum()), "route_alone": B.auroc(route[:, 0], y)}
        for name, x in (("protein+readout", base), ("protein+readout+route", np.hstack([base, route]))):
            pred = np.zeros(len(y))
            for tr, te in GroupKFold(n_splits=5).split(x, y, g):
                pred[te] = LogisticRegression(max_iter=1000).fit(x[tr], y[tr]).predict_proba(x[te])[:, 1]
            res[name] = B.auroc(pred, y)
        out[f"{mult}delta"] = res
    return out


def s2(m0):
    genes = sorted({g for g, _ in m0})
    ros = sorted({r for _, r in m0})
    x = np.array([[m0[(g, r)][0] / m0[(g, r)][1] if (g, r) in m0 else 0.0 for r in ros] for g in genes])
    sv = np.linalg.svd(x, compute_uv=False)
    var = sv ** 2 / (sv ** 2).sum()
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    first = sorted(zip(vt[0] * np.sign(vt[0].sum()), ros), reverse=True)
    return {"proteins": len(genes), "readouts": len(ros), "variance_explained": [round(float(v), 3) for v in var[:6]],
            "first_direction_loadings": [(r, round(float(w), 2)) for w, r in first]}


def main():
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, pred_file, build, out = sys.argv[1:5]
    genes = A.marker_genes()
    m0 = B.measured(a_file, genes)
    pred = B.predictions(pred_file)
    cases = [{"gene": k[0], "readout": k[1], "measured": B.call_of(*m0[k]), **pred[k]} for k in sorted(pred) if k in m0]
    assert all(split_of(c["gene"]) == "A" for c in cases)
    res = {"S1": s1(m0, cases), "S2": s2(m0), "T3": t3(cases, pred_file, build)}
    Path(out).write_text(json.dumps(res, indent=1, default=str))
    res["T3"].pop("pairs_with_a_form_route")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
