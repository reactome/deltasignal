#!/usr/bin/env python3
"""specs/176 post hoc: is the dominant response axis "an over-expressed tyrosine kinase raises tyrosine sites"?

Protein class from the networks (followup_176.enzyme_classes: curated catalyst of a reaction that adds a
phosphotyrosine). Readout class from the canonical sequence: the residue at each registered site.

Usage: axis_176.py A_FILE PREDICTIONS BUILD FASTA OUT_JSON
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
sys.path.insert(0, str(HERE.parents[1] / "174"))
import budget_176 as B  # noqa: E402
import followup_176 as F  # noqa: E402
import lun_a_173 as A  # noqa: E402
from trace_174 import READOUT_ACC  # noqa: E402


def sequences(fasta, wanted):
    seqs, acc = {}, None
    for ln in open(fasta):
        if ln.startswith(">"):
            acc = ln.split("|")[1] if ln.startswith(">sp|") else None
            if acc in wanted:
                seqs[acc] = []
            else:
                acc = None
        elif acc:
            seqs[acc].append(ln.strip())
    return {a: "".join(s) for a, s in seqs.items()}


def readout_residues(fasta):
    """{readout: residues at its sites, e.g. 'Y' or 'TY'}; a site whose residue is not S/T/Y is reported as '?'."""
    rows = [r for r in A.read_tsv(A.MAPS / "readouts.tsv") if r["kind"] == "phospho"]
    seqs = sequences(fasta, set(READOUT_ACC.values()))
    out = {}
    for r in rows:
        res = set()
        for s in r["sites"].split(","):
            g, pos = s.split(":")
            seq = seqs.get(READOUT_ACC.get(g), "")
            aa = seq[int(pos) - 1] if len(seq) >= int(pos) else "?"
            res.add(aa if aa in "STY" else "?")
        out[r["readout"]] = "".join(sorted(res))
    return out


def main():
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, pred_file, build, fasta, out = sys.argv[1:6]
    genes = A.marker_genes()
    m0 = B.measured(a_file, genes)
    pred = B.predictions(pred_file)
    acc = {r["gene"]: r["accession"] for r in A.read_tsv(A.MAPS / "pois.tsv")}
    cls = F.enzyme_classes(build, {acc[g] for g, _ in m0 if acc.get(g)})
    res_of = readout_residues(fasta)
    tk = {g: "tyr_kinase" in cls.get(acc.get(g), ()) for g, _ in m0}
    tab, in_p = defaultdict(Counter), defaultdict(Counter)
    for (g, r), (e, d) in m0.items():
        if r not in res_of:
            continue
        site = "tyrosine only" if res_of[r] == "Y" else "tyrosine + S/T" if "Y" in res_of[r] else "S/T only"
        key = ("tyrosine kinase" if tk[g] else "other protein", site)
        tab[key][B.call_of(e, d)] += 1
        if (g, r) in pred:
            in_p[key][(B.call_of(e, d), pred[(g, r)]["reach"] == "up", pred[(g, r)]["frozen"])] += 1
    res = {"readout_residues": res_of, "tyrosine_kinases": sorted(g for g, v in tk.items() if v), "table": {}}
    for key, c in sorted(tab.items()):
        n = sum(c.values())
        p = in_p[key]
        npairs = sum(p.values())
        res["table"][" x ".join(key)] = {
            "pairs": n, "up": round(c["up"] / n, 3), "down": round(c["down"] / n, 3),
            "scored_pairs": npairs,
            "scored_with_route": sum(v for (m, rt, f), v in p.items() if rt),
            "scored_up": sum(v for (m, rt, f), v in p.items() if m == "up"),
            "scored_up_called_up": sum(v for (m, rt, f), v in p.items() if m == "up" and f == "up")}
    # sizing: the model's calls plus "tyrosine kinase x tyrosine-containing site -> UP" on the scored pairs
    cases = []
    for k in sorted(pred):
        if k in m0:
            rule = tk[k[0]] and "Y" in res_of.get(k[1], "")
            cases.append({"gene": k[0], "readout": k[1], "measured": B.call_of(*m0[k]), "reach": pred[k]["reach"],
                          "frozen": pred[k]["frozen"], "reach_call": "up" if pred[k]["reach"] == "up" else "not",
                          "plus_rule": "up" if rule and pred[k]["frozen"] == "not" else pred[k]["frozen"],
                          "rule_only": "up" if rule else "not"})
    res["sizing"] = {name: {k: v for k, v in B.scores(cases, name).items()
                            if k in ("calls", "recall", "precision_exact", "balanced_accuracy")}
                     for name in ("frozen", "plus_rule", "rule_only")}
    import lun_score as L
    for x, y in (("plus_rule", "frozen"), ("plus_rule", "reach_call"), ("rule_only", "reach_call")):
        pt, lo, hi = L.bootstrap_difference(cases, x, y, seed=176)
        res["sizing"][f"ba_{x}-{y}"] = {"point": pt, "ci": [lo, hi]}
    per = defaultdict(lambda: [0, 0])
    for (g, r), (e, d) in m0.items():
        if tk[g] and "Y" in res_of.get(r, ""):
            per[g][0] += B.call_of(e, d) == "up"
            per[g][1] += 1
    res["per_tyrosine_kinase_up_of_tyr_readouts"] = {g: f"{u}/{n}" for g, (u, n) in sorted(per.items(), key=lambda kv: -kv[1][0])}
    res["tyrosine_kinases_raising_half_or_more"] = sum(u >= n / 2 for u, n in per.values())
    added = [c for c in cases if c["plus_rule"] != c["frozen"]]
    res["sizing"]["added_calls"] = {"n": len(added), "right": sum(c["measured"] == "up" for c in added),
                                    "measured_not": sum(c["measured"] == "not" for c in added),
                                    "measured_down": sum(c["measured"] == "down" for c in added)}
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
