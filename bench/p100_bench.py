#!/usr/bin/env python3
"""specs/113: reachability, with and without Reactome's interaction records, on LINCS P100.

The protocol is pre-registered in specs/113-interaction-edges/spec.md. In short:

- `predict` (no P100 value is read): each probe is scored by its protein for each compound with targets
  (p100_targets.py):
    A  1 if, in one network of the frozen catalog, a node holding a target reaches a node holding the probe's
       protein along logic_network.csv edges (length >= 0), else 0;
    B  2 if A; else 1 if A succeeds once the targets are widened by their IntAct partners and the probe by its
       protein's IntAct partners; else 0;
    C  IntAct distance from a target to the probe's protein: 0-1 hops 2, 2 hops 1, else 0.
- `score`: per (compound, cell line) at the highest dose, a probe is changed if its replicate median is at least
  3 x 1.4826 x MAD of the cell line's DMSO wells from their median. AUROC per context (ties half), mean over
  contexts, a bootstrap over compounds, and a random-interaction control for B.

Usage:
  p100_bench.py predict CATALOG TARGETS_TSV ANALYTE_GZ INTERACTIONS_CSV OUT_DIR
  p100_bench.py score INST_GZ LEVEL4_GCTX OUT_DIR [CATALOG INTERACTIONS_CSV TARGETS_TSV ANALYTE_GZ]
  (the four extra arguments add the random-interaction control)

INTERACTIONS_CSV is the Neo4j export (interaction dbId, accession, taxId); only human interactors are used.
`score` needs h5py and numpy.
"""
import collections
import csv
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dataset_roles  # noqa: E402
from hpn_dream_bench import auroc  # noqa: E402

DATASET = "lincs_p100"
HUMAN = "9606"
MAD_K = 3 * 1.4826
BOOT, SEED, CONTROLS = 10000, 113, 200


def acc(x):
    return x.strip().strip('"').split("-")[0]


def load_partners(path):
    """accession -> set of other human accessions sharing an interaction."""
    members = collections.defaultdict(set)
    with open(path) as f:
        next(f)
        for line in f:
            i, a, tax = (x.strip().strip('"') for x in line.split(","))
            if tax == HUMAN:
                members[i].add(acc(a))
    partners = collections.defaultdict(set)
    for ms in members.values():
        for a in ms:
            partners[a] |= ms - {a}
    return dict(partners)


def network_reach(catalog):
    """accession -> set of accessions it reaches, the union over networks (each searched on its own)."""
    reach = collections.defaultdict(set)
    for net in sorted(Path(catalog).glob("R-HSA-*")):
        leaf = {}
        for r in csv.DictReader(open(net / "leaf_refs.csv")):
            if r["reference_stid"].lower().startswith("uniprot:"):
                leaf[r["leaf_stid"]] = acc(r["reference_stid"].split(":", 1)[1])
        holds = {}
        for r in csv.DictReader(open(net / "nodes.csv")):
            ps = {leaf[x] for x in r["member_leaves"].split("|") if x in leaf}
            if ps:
                holds[r["uuid"]] = ps
        out = collections.defaultdict(set)
        for r in csv.DictReader(open(net / "logic_network.csv")):
            out[r["source_id"]].add(r["target_id"])
        nodes_of = collections.defaultdict(set)
        for n, ps in holds.items():
            for p in ps:
                nodes_of[p].add(n)
        for p, start in nodes_of.items():
            seen, stack = set(start), list(start)
            while stack:
                for m in out.get(stack.pop(), ()):
                    if m not in seen:
                        seen.add(m)
                        stack.append(m)
            for n in seen:
                reach[p] |= holds.get(n, set())
    return dict(reach)


def reached(sources, reach):
    out = set()
    for s in sources:
        out |= reach.get(s, set())
    return out


def score_a_b(targets, probe, reach, partners):
    return scores_a_b(targets, [probe], reach, partners)[0]


def scores_a_b(targets, prots, reach, partners):
    """[(A, B)] for each protein in `prots`; reach is computed once per compound."""
    direct = reached(targets, reach)
    wide = reached(set(targets) | set().union(*(partners.get(t, set()) for t in targets)), reach)
    out = []
    for p in prots:
        if p in direct:
            out.append((1, 2))
        else:
            out.append((0, int(p in wide or bool(partners.get(p, set()) & wide))))
    return out


def score_c(targets, probe, partners):
    one = set(targets) | set().union(*(partners.get(t, set()) for t in targets))
    if probe in one:
        return 2
    two = one | set().union(*(partners.get(x, set()) for x in one))
    return int(probe in two)


def load_targets(path):
    t = collections.defaultdict(set)
    for r in csv.DictReader(open(path), delimiter="\t"):
        t[r["pert_id"]].add(acc(r["uniprot"]))
    return dict(t)


def probes(analyte_gz):
    with gzip.open(analyte_gz, "rt") as f:
        return {r["pr_analyte_id"]: acc(r["pr_uniprot_id"]) for r in csv.DictReader(f, delimiter="\t")}


def predict(catalog, targets_tsv, analyte_gz, interactions_csv, out_dir):
    dataset_roles.require_role(DATASET, "score")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    reach, partners = network_reach(catalog), load_partners(interactions_csv)
    tg, pr = load_targets(targets_tsv), probes(analyte_gz)
    with open(out / "p100_predictions.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["pert_id", "probe", "uniprot", "A", "B", "C"])
        for pid in sorted(tg):
            for probe, p in sorted(pr.items()):
                a, b = score_a_b(tg[pid], p, reach, partners)
                w.writerow([pid, probe, p, a, b, score_c(tg[pid], p, partners)])
    summary = {"compounds": len(tg), "probes": len(pr), "network_proteins": len(reach),
               "compounds_with_target_in_networks": sum(bool(t & reach.keys()) for t in tg.values())}
    (out / "predict.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


def permuted(partners, rng):
    """Partner sets relabelled by a permutation of the proteins that have interactions (degrees kept)."""
    ids = sorted(partners)
    pi = dict(zip(ids, (ids[i] for i in rng.permutation(len(ids)))))
    return {pi[a]: {pi[b] for b in bs} for a, bs in partners.items()}


def read_gctx(path):
    """(matrix as {(col_id, row_id): value}, row ids, col ids) from a GCTX (HDF5) file."""
    import h5py
    import numpy as np
    with h5py.File(path, "r") as h:
        m = np.array(h["0/DATA/0/matrix"])
        rows = [x.decode() for x in h["0/META/ROW/id"][:]]
        cols = [x.decode() for x in h["0/META/COL/id"][:]]
    if m.shape == (len(cols), len(rows)):
        return m, rows, cols
    if m.shape == (len(rows), len(cols)):
        return m.T, rows, cols
    raise SystemExit(f"matrix {m.shape} does not fit {len(rows)} rows x {len(cols)} cols")


def truth(inst_gz, level4_gctx):
    """{(pert_id, cell): {probe: 0/1}} at each context's highest dose, against the cell line's DMSO wells."""
    import numpy as np
    m, rows, cols = read_gctx(level4_gctx)
    col = {c: i for i, c in enumerate(cols)}
    with gzip.open(inst_gz, "rt") as f:
        all_inst = list(csv.DictReader(f, delimiter="\t"))
    inst = [r for r in all_inst if r["inst_id"] in col]
    if len(inst) != len(cols) or len(inst) != len(all_inst):
        raise SystemExit(f"{len(inst)} instances match the matrix's {len(cols)} columns ({len(all_inst)} listed)")
    dmso = collections.defaultdict(list)
    top = {}
    for r in inst:
        if r["pert_type"] == "ctl_vehicle":
            dmso[r["cell_id"]].append(col[r["inst_id"]])
        elif r["pert_type"] == "trt_cp":
            k, d = (r["pert_id"], r["cell_id"]), float(r["pert_dose"])
            if k not in top or d > top[k][0]:
                top[k] = (d, [])
            if d == top[k][0]:
                top[k][1].append(col[r["inst_id"]])
    out = {}
    for (pid, cell), (_, idx) in sorted(top.items()):
        ref = m[dmso[cell]]
        vec = {}
        for j, probe in enumerate(rows):
            ctl = ref[:, j][~np.isnan(ref[:, j])]
            vals = m[idx, j][~np.isnan(m[idx, j])]
            if len(ctl) < 3 or len(vals) == 0:
                continue
            med = np.median(ctl)
            mad = np.median(np.abs(ctl - med))
            if mad == 0:
                continue
            vec[probe] = int(abs(np.median(vals) - med) >= MAD_K * mad)
        out[(pid, cell)] = vec
    return out


def contexts_auroc(pred, tv, col):
    res = []
    for (pid, cell), vec in tv.items():
        if pid not in pred:
            continue
        keep = [p for p in sorted(vec) if p in pred[pid]]
        a = auroc([pred[pid][p][col] for p in keep], [vec[p] for p in keep])
        if a is not None:
            res.append({"pert_id": pid, "cell": cell, "n": len(keep), "positives": sum(vec[p] for p in keep),
                        col: a})
    return res


def score(inst_gz, level4_gctx, out_dir, catalog=None, interactions_csv=None, targets_tsv=None, analyte_gz=None):
    import numpy as np
    role = dataset_roles.require_role(DATASET, "score")
    out = Path(out_dir)
    pred = collections.defaultdict(dict)
    for r in csv.DictReader(open(out / "p100_predictions.tsv"), delimiter="\t"):
        pred[r["pert_id"]][r["probe"]] = {k: float(r[k]) for k in ("A", "B", "C")}
    tv = truth(inst_gz, level4_gctx)
    names = {n for v in pred.values() for n in v}
    matrix_probes = {n for v in tv.values() for n in v}
    if not matrix_probes <= names or len(names) != 96:
        raise SystemExit(f"probe ids differ: {len(matrix_probes - names)} matrix probes not predicted")
    missing = sorted({pid for pid, _ in tv if pid in pred} ^ set(pred))
    if missing:
        raise SystemExit(f"{len(missing)} predicted compounds have no context")
    per = {}
    for col in ("A", "B", "C"):
        for r in contexts_auroc(pred, tv, col):
            per.setdefault((r["pert_id"], r["cell"]), {k: v for k, v in r.items() if k not in ("A", "B", "C")})[col] = r[col]
    ctx = sorted(per)
    comps = sorted({p for p, _ in ctx})
    by_comp = {c: [per[k] for k in ctx if k[0] == c] for c in comps}
    mean = lambda rs, k: float(np.mean([r[k] for r in rs]))  # noqa: E731
    rng = np.random.default_rng(SEED)
    boot = {k: [] for k in ("A", "B", "C", "B-A")}
    for _ in range(BOOT):
        rs = [r for i in rng.integers(0, len(comps), len(comps)) for r in by_comp[comps[i]]]
        for k in ("A", "B", "C"):
            boot[k].append(mean(rs, k))
        boot["B-A"].append(boot["B"][-1] - boot["A"][-1])
    rows = list(per.values())
    summary = {"dataset": DATASET, "role": role, "contexts": len(rows), "compounds": len(comps),
               "mean_auroc": {k: mean(rows, k) for k in ("A", "B", "C")},
               "B_minus_A": mean(rows, "B") - mean(rows, "A"),
               "interval_95": {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
                               for k, v in boot.items()},
               "boot": BOOT, "seed": SEED}
    if catalog:
        reach, partners = network_reach(catalog), load_partners(interactions_csv)
        tg, pr = load_targets(targets_tsv), probes(analyte_gz)
        crng, draws = np.random.default_rng(SEED), []
        for _ in range(CONTROLS):
            pp = permuted(partners, crng)
            names = sorted(pr)
            cp = {pid: {n: {"B": float(ab[1])} for n, ab in zip(names, scores_a_b(tg[pid], [pr[n] for n in names],
                                                                                 reach, pp))}
                  for pid in pred}
            draws.append(mean(contexts_auroc(cp, {k: tv[k] for k in ctx}, "B"), "B"))
        b = summary["mean_auroc"]["B"]
        summary["random_interaction_control"] = {
            "draws": CONTROLS, "seed": SEED, "mean": float(np.mean(draws)),
            "p95": float(np.percentile(draws, 95)), "B_percentile": float(np.mean([d < b for d in draws]) * 100)}
    summary["per_context"] = rows
    (out / "score.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "per_context"}, indent=1))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "predict" and len(sys.argv) == 7:
        return predict(*sys.argv[2:7])
    if mode == "score" and len(sys.argv) in (5, 9):
        return score(*sys.argv[2:5], *sys.argv[5:9])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
