#!/usr/bin/env python3
"""specs/181 part K: enzyme–substrate edges (OmniPath) as reachability on Lun 2019. Both halves are FIT.

Usage: ks_181.py PREDICTIONS ENZSUB_TSV OUT_JSON A_FILE B_FILE [CLASS_DIR]
"""
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
for d in ("173", "176", "177", "179"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
import score_177 as S  # noqa: E402

SEED, NULL_DRAWS, BOOT = 181, 1000, 10_000
MODS = {"phosphorylation": "up", "dephosphorylation": "down"}
ARMS = ("K1", "K1p", "K2")
MIN_COVERAGE = 0.20


def read_enzsub(path):
    """{enzyme gene: {(substrate gene, position, direction)}} for phosphorylation / dephosphorylation rows."""
    out = defaultdict(set)
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["modification"] in MODS and r["enzyme_genesymbol"] != r["substrate_genesymbol"]:
            out[r["enzyme_genesymbol"]].add((r["substrate_genesymbol"], int(r["residue_offset"]), MODS[r["modification"]]))
    return out


def readout_sites(path=A.MAPS / "readouts.tsv"):
    """{readout: {(gene, position)}} for phospho readouts."""
    out = {}
    for r in A.read_tsv(path):
        if r["kind"] == "phospho":
            out[r["readout"]] = {(g, int(p)) for g, p in (s.split(":") for s in r["sites"].split(","))}
    return out


def one(dirs):
    return next(iter(dirs)) if len(dirs) == 1 else None


def arm_calls(ks, sites, reach_up, lun_genes):
    """{arm: {(gene, readout): direction}}. reach_up: {gene: readouts it reaches by curated routes}."""
    by_gene = defaultdict(set)
    for ro, ss in sites.items():
        for g, _ in ss:
            by_gene[g].add(ro)
    k1, k1p, k2 = defaultdict(set), defaultdict(set), {}
    for e, subs in ks.items():
        if e not in lun_genes:
            continue
        for g, pos, d in subs:
            for ro in by_gene.get(g, ()):
                k1p[(e, ro)].add(d)
                if (g, pos) in sites[ro]:
                    k1[(e, ro)].add(d)
            if d == "up" and g in lun_genes:
                for ro in reach_up.get(g, ()):
                    k2[(e, ro)] = "up"
    return {"K1": {k: one(v) for k, v in k1.items() if one(v)}, "K1p": {k: one(v) for k, v in k1p.items() if one(v)},
            "K2": k2}


def added(cases, calls):
    """Cases the arm calls where curated reachability has no route."""
    return [dict(c, arm=calls[(c["gene"], c["readout"])]) for c in cases
            if c["reach"] != "up" and (c["gene"], c["readout"]) in calls]


def precision(add):
    return sum(c["arm"] == c["measured"] for c in add) / len(add) if add else None


def relabel(ks, lun_genes, rng):
    """Substrate lists permuted among Lun enzymes of the same type (has a phosphorylation row or not)."""
    out = dict(ks)
    for is_kinase in (True, False):
        es = sorted(e for e in ks if e in lun_genes and any(d == "up" for *_, d in ks[e]) == is_kinase)
        sh = es[:]
        rng.shuffle(sh)
        out.update({a: ks[b] for a, b in zip(es, sh)})
    return out


def with_arm(cases, calls, base):
    return [dict(c, plus=c[base] if B.changed(c[base]) else calls.get((c["gene"], c["readout"]), "not")) for c in cases]


def analyse(cases, ks, sites, lun_genes):
    import lun_score as L
    reach_up = defaultdict(set)
    for c in cases:
        if c["reach"] == "up":
            reach_up[c["gene"]].add(c["readout"])
    noroute = [c for c in cases if c["reach"] != "up"]
    real = [c for c in noroute if B.changed(c["measured"])]
    res = {"pairs": len(cases), "no_route_pairs": len(noroute), "no_route_real_changes": len(real),
           "base_rate_changed": len(real) / len(noroute), "base_rate_up": sum(c["measured"] == "up" for c in noroute) / len(noroute),
           "lun_enzymes_in_table": len(set(ks) & lun_genes)}
    calls = arm_calls(ks, sites, reach_up, lun_genes)
    rng = random.Random(SEED)
    nulls = defaultdict(list)
    for _ in range(NULL_DRAWS):
        nc = arm_calls(relabel(ks, lun_genes, rng), sites, reach_up, lun_genes)
        for arm in ARMS:
            p = precision(added(cases, nc[arm]))
            if p is not None:
                nulls[arm].append(p)
    for arm in ARMS:
        add = added(cases, calls[arm])
        covered = sum(B.changed(c["measured"]) for c in add)
        p = precision(add)
        r = {"calls_total": len(calls[arm]), "calls_on_scored_pairs": sum((c["gene"], c["readout"]) in calls[arm] for c in cases),
             "already_curated_reach": sum(c["reach"] == "up" and (c["gene"], c["readout"]) in calls[arm] for c in cases),
             "added_pairs": len(add), "added_proteins": len({c["gene"] for c in add}),
             "covered_real_changes": covered, "coverage": covered / len(real) if real else None,
             "added_changed": covered / len(add) if add else None, "added_direction_right": p,
             "added_by_direction": {d: {"n": sum(c["arm"] == d for c in add), "right": sum(c["arm"] == d == c["measured"] for c in add)}
                                    for d in ("up", "down")},
             "null_p95": float(np.percentile(nulls[arm], 95)) if nulls[arm] else None,
             "null_median": float(np.median(nulls[arm])) if nulls[arm] else None}
        r["above_null"] = bool(p is not None and r["null_p95"] is not None and p > r["null_p95"])
        for base in ("candidate", "model"):
            cs = with_arm(cases, calls[arm], base)
            pt, lo, hi = L.bootstrap_difference(cs, "plus", base, seed=SEED)
            r[f"balanced_accuracy_{base}_plus_arm"] = L.balanced_accuracy([(c["measured"], c["plus"]) for c in cs])
            r[f"gain_over_{base}"] = {"point": pt, "ci": [lo, hi]}
        res[arm] = r
    k1 = res["K1"]
    res["registered_reading"] = {"K1_coverage_at_least_20pct": bool((k1["coverage"] or 0) >= MIN_COVERAGE),
                                 "K1_above_null": k1["above_null"],
                                 "covers_the_gap": bool((k1["coverage"] or 0) >= MIN_COVERAGE and k1["above_null"])}
    return res


def main():
    import dataset_roles
    pred_file, enzsub, out, a_file, b_file = sys.argv[1:6]
    tk, tyr = S.read_classes(sys.argv[6]) if len(sys.argv) > 6 else S.read_classes()
    genes = A.marker_genes()
    lun_genes = set(genes.values())
    ks, sites = read_enzsub(enzsub), readout_sites()
    res, allc = {}, []
    for half, path in (("A", a_file), ("B", b_file)):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        cases = S.build_cases(S.measured(path, genes, half), S.predictions(pred_file, half), tk, tyr)
        res[f"half_{half}"] = analyse(cases, ks, sites, lun_genes)
        allc += cases
    res["pooled"] = analyse(allc, ks, sites, lun_genes)
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
