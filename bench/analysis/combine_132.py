#!/usr/bin/env python3
"""specs/132: score the curated-connection arms against the control, as registered.

Arms are run by `scripts/run_arm.sh` on bundles built by `context_build.py --curated`: the control `c0` (each pathway
alone, same build) and the rules (L1, L2). For each arm this reports, on cases valid (scored) in both arms:

- the ceiling count: of the control's wrong cases, how many have no route in the control but a route in the arm
  (experimental and curator: reachability from the case's pins to its readout in the bundle; korkut and decryptM: the
  bench's reachability column);
- the fallback F: the control's answer wherever the control has a route, the arm's elsewhere;
- the candidate gate: three-class balanced accuracy over experimental (RAF set aside) and korkut pooled, arm minus
  control, with a 95% interval from resampling perturbations (experimental: gene; korkut: drug set) within each
  dataset, 10,000 times, seed 132; and each dataset's own difference. decryptM (by drug) and curator held-out are
  reported, never gated.

Usage: combine_132.py CONTROL_RESULTS ARM_RESULTS CONTROL_BUNDLE ARM_BUNDLE > out.json
"""
import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
from holdout_report import TUNING_PATHWAYS  # noqa: E402

CLASSES = ("up", "not", "down")
CODE = {"0": "down", "1": "not", "2": "up"}
RAF = "RAF_MAP_kinase_cascade"
SEED = 132
BOOTSTRAP = 10000
CATALOG = Path(__file__).resolve().parents[1] / "catalog_pathways.tsv"


def pathway_ids():
    return {r["pathway_name"]: r["id"] for r in csv.DictReader((ln for ln in open(CATALOG) if not ln.startswith("#")),
                                                               delimiter="\t")}


class Reach:
    """Route from a case's pins to its readout in one bundle, adjacency cached per pathway."""

    def __init__(self, bundle):
        self.bundle, self.adj, self.ids = Path(bundle), {}, pathway_ids()

    def __call__(self, pathway, pins, readout):
        pid = self.ids[pathway]
        if pid not in self.adj:
            self.adj[pid] = bm.build_adjacency(self.bundle / pid)
        return bool(set(readout) & bm.reachable_from(self.adj[pid], set(pins)))


def load_cell_cases(path, reach, kind):
    """{case: dict(cluster, measured, predicted, route, pathway)} over valid cases (experimental or curator)."""
    out = {}
    for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
        if r["valid"] != "1":
            continue
        pins, ro = [u for u in r["gene_uuids"].split("|") if u], [u for u in r["output_uuids"].split("|") if u]
        out[(r["pathway"], r["gene"], r["direction"], r["key_output"])] = {
            "cluster": r["gene"] if kind == "experimental" else (r["gene"], r["direction"]),
            "measured": CODE[r["expected"]], "predicted": CODE[r["predicted"]],
            "route": reach(r["pathway"], pins, ro), "pathway": r["pathway"]}
    return out


def load_drug_cases(path, kind):
    """{case: ...} over scored cases with a measured call (korkut or decryptM)."""
    out = {}
    for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
        if r["status"] != "scored" or r["measured"] not in CLASSES:
            continue
        if kind == "korkut":
            key, cluster = (r["pathway_id"], r["drug_set"], r["readout"]), r["drug_set"]
        else:
            key, cluster = (r["pathway_id"], r["drug"], r["uniprot"], r["position"]), r["drug"]
        out[key] = {"cluster": cluster, "measured": r["measured"], "predicted": r["predicted"] or "not",
                    "route": (r["reach_predicted"] or "not") != "not", "pathway": r["pathway_id"]}
    return out


def fallback(control, arm):
    """F: the control's answer where the control routes, the arm's elsewhere (keys valid in both)."""
    return {k: dict(arm[k], predicted=control[k]["predicted"]) if control[k]["route"] else arm[k]
            for k in control.keys() & arm.keys()}


def ceiling(control, arm):
    """(control wrong, of which no route in the control but a route in the arm) over keys in both."""
    wrong = [k for k in control.keys() & arm.keys() if control[k]["predicted"] != control[k]["measured"]]
    return len(wrong), sum(1 for k in wrong if not control[k]["route"] and arm[k]["route"])


def balanced_accuracy(pairs):
    hit, tot = Counter(), Counter()
    for m, p in pairs:
        tot[m] += 1
        hit[m] += m == p
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def gain(datasets, control, arm, keys):
    """Pooled balanced accuracy of `arm` minus `control` over the listed datasets and keys."""
    a = [(arm[d][k]["measured"], arm[d][k]["predicted"]) for d in datasets for k in keys[d]]
    c = [(control[d][k]["measured"], control[d][k]["predicted"]) for d in datasets for k in keys[d]]
    return balanced_accuracy(a) - balanced_accuracy(c)


def bootstrap(datasets, control, arm, n=BOOTSTRAP, seed=SEED):
    """95% interval of the pooled gain, perturbations resampled within each dataset."""
    by = {d: defaultdict(list) for d in datasets}
    for d in datasets:
        for k in control[d].keys() & arm[d].keys():
            by[d][control[d][k]["cluster"]].append(k)
    rng, draws = random.Random(seed), []
    clusters = {d: sorted(by[d], key=str) for d in datasets}
    for _ in range(n):
        keys = {d: [k for c in rng.choices(clusters[d], k=len(clusters[d])) for k in by[d][c]] for d in datasets}
        draws.append(gain(datasets, control, arm, keys))
    draws.sort()
    return draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


def report(control, arm, n=BOOTSTRAP):
    """control, arm: {dataset: cases}. Returns the registered figures for the arm and for arm + F."""
    out = {}
    for name, a in (("rule", arm), ("rule+F", {d: fallback(control[d], arm[d]) for d in arm})):
        res = {}
        for d in control:
            keys = sorted(control[d].keys() & a[d].keys(), key=str)
            wrong, reached = ceiling(control[d], a[d])
            moved = [k for k in keys if a[d][k]["predicted"] != control[d][k]["predicted"]]
            res[d] = {"cases": len(keys), "control_wrong": wrong, "ceiling": reached,
                      "ceiling_share": reached / wrong if wrong else None,
                      "changed": len(moved),
                      "fixed": sum(a[d][k]["predicted"] == a[d][k]["measured"] for k in moved),
                      "broken": sum(control[d][k]["predicted"] == control[d][k]["measured"] for k in moved),
                      "balanced_gain": gain([d], control, a, {d: keys})}
        gated = [d for d in ("experimental", "korkut") if d in control]
        keys = {d: sorted(control[d].keys() & a[d].keys(), key=str) for d in gated}
        g = gain(gated, control, a, keys)
        lo, hi = bootstrap(gated, control, a, n=n)
        res["gate"] = {"pooled_gain": g, "ci95": [lo, hi],
                       "passes": lo > 0 and all(res[d]["balanced_gain"] >= 0 for d in gated)}
        out[name] = res
    return out


def load_all(results, bundle):
    r, reach = Path(results), Reach(bundle)
    data = {}
    exp = load_cell_cases(r / "experimental_cases.tsv", reach, "experimental")
    data["experimental"] = {k: v for k, v in exp.items() if v["pathway"] != RAF}
    data["experimental_with_raf"] = exp
    cur = load_cell_cases(r / "curator_cases.tsv", reach, "curator")
    data["curator_heldout"] = {k: v for k, v in cur.items() if v["pathway"] not in TUNING_PATHWAYS}
    data["korkut"] = load_drug_cases(r / "korkut_cases.tsv", "korkut")
    data["decryptm"] = load_drug_cases(r / "decryptm" / "decryptm_cases.tsv", "decryptm")
    return data


def main():
    c_res, a_res, c_bundle, a_bundle = sys.argv[1:]
    control, arm = load_all(c_res, c_bundle), load_all(a_res, a_bundle)
    json.dump(report(control, arm), sys.stdout, indent=1, default=str)
    print()


if __name__ == "__main__":
    main()
