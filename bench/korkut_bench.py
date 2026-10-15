#!/usr/bin/env python3
"""specs/062: score the solver against korkut2015 (SkMel-133, 12 drugs and
pairs, RPPA readouts). The protocol is registered in specs/062-korkut-test/spec.md.

- A case is (drug set, pathway, readout). A drug set is one row of
  expr_index.txt without doses; its target genes (bench/analysis/062/drugs.tsv)
  are pinned to DS_PERTURB_UI_DOWN at the forms the curator bench pins
  (`pin_occurrences`), one solve per drug set per pathway.
- The readout (bench/analysis/062/readouts.tsv) is the mean over the nodes that
  are, or contain, an EWAS of the gene (total) or of the gene phosphorylated at
  the coordinate (phospho), pinned nodes left out.
- `--predict-only` never opens expr.csv: it writes predictions only, for the
  registered values-blind step. Without it, each case's measured call is made
  from expr.csv by the registered rule.

Usage: DS_CATALOG_ROOT=<build> korkut_bench.py [--predict-only] OUT_DIR
"""
import collections
import csv
import json
import math
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import dataset_roles  # noqa: E402
from decryptm_bench import git_state  # noqa: E402

DATA = Path(os.environ.get("KORKUT_DIR", Path.home() / "deltasignal-catalogs/external/korkut2015"))
MAPS = Path(__file__).resolve().parent / "analysis" / "062"
THRESHOLD = math.log2(1.5)        # registered primary
THRESHOLD_SECONDARY = math.log2(1.15)
NAMES = {bm.DOWN: "down", bm.NORMAL: "not", bm.UP: "up"}


def read_tsv(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def conditions(path=None):
    """[(row index, drug set as a sorted tuple of codes)] in expr_index.txt order."""
    out = []
    for i, line in enumerate(open(path or DATA / "expr_index.txt")):
        label = line.split()[0]
        out.append((i, tuple(sorted(part.split("|")[0] for part in label.split(",")))))
    return out


def call(value, threshold):
    if value >= threshold:
        return "up"
    if value <= -threshold:
        return "down"
    return "not"


def drug_set_call(values, threshold):
    """Several conditions (doses) of one drug set -> one call: agreeing regulated
    calls give that call, disagreeing ones `conflict`, none `not` (specs/062)."""
    called = {c for c in (call(v, threshold) for v in values) if c != "not"}
    return "conflict" if len(called) > 1 else (called.pop() if called else "not")


def neo4j_site_stids(sites):
    """{(gene, coordinate): [EWAS stIds]} for human EWAS phosphorylated there."""
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $sites AS s
        MATCH (e:EntityWithAccessionedSequence)-[:referenceEntity]->(r:ReferenceGeneProduct),
              (e)-[:hasModifiedResidue]->(m)-[:psiMod]->(p)
        WHERE e.speciesName = 'Homo sapiens' AND s[0] IN r.geneName AND m.coordinate = s[1]
          AND p.name[0] IN ['O-phospho-L-serine', 'O-phospho-L-threonine', "O4'-phospho-L-tyrosine"]
        RETURN s[0] AS gene, s[1] AS coord, COLLECT(DISTINCT e.stId) AS stids
        """,
        sites=[[g, c] for g, c in sites],
    ).data()
    return {(r["gene"], r["coord"]): r["stids"] for r in rows}


def readout_stids(readouts, gene_to_stids, site_to_stids):
    """{readout: [stIds]}: the EWAS whose nodes (and containers) are read."""
    out = {}
    for r in readouts:
        if r["kind"] == "total":
            out[r["readout"]] = sorted({s for g in r["sites"].split(",") for s in gene_to_stids.get(g, [])})
        elif r["kind"] == "phospho":
            out[r["readout"]] = sorted({s for gc in r["sites"].split(",")
                                        for s in site_to_stids.get((gc.split(":")[0], int(gc.split(":")[1])), [])})
    return out


def measured_calls(threshold, design):
    """{(drug set, readout): call} from expr.csv. Only called without --predict-only."""
    names = [l.strip() for l in open(DATA / "node_Index.csv")]
    rows = list(csv.reader(open(DATA / "expr.csv")))
    if len(rows) != len(design):
        raise SystemExit(f"expr.csv has {len(rows)} rows, expr_index.txt {len(design)}: stop and amend (specs/062)")
    by = collections.defaultdict(list)
    for (i, ds), row in zip(design, rows):
        for name, v in zip(names, row):
            by[(ds, name)].append(float(v))
    return {k: drug_set_call(v, threshold) for k, v in by.items()}


def main():
    args = sys.argv[1:]
    predict_only = "--predict-only" in args
    args = [a for a in args if a != "--predict-only"]
    out_dir = Path(args[0])
    role = dataset_roles.require_role("korkut2015", "score")   # specs/051
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")

    drugs = {d["code"]: d["genes"].split(",") for d in read_tsv(MAPS / "drugs.tsv")}
    readouts = [r for r in read_tsv(MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    design = conditions()
    drug_sets = sorted({ds for _, ds in design})
    unknown = {c for ds in drug_sets for c in ds} - set(drugs)
    if unknown:
        raise SystemExit(f"drug codes not in drugs.tsv: {sorted(unknown)}")
    target_genes = sorted({g for gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    site_to_stids = neo4j_site_stids(sites)
    r_stids = readout_stids(readouts, gene_to_stids, site_to_stids)

    print(f"Protocol: korkut2015 drug-set bench (specs/062), DS_PIN_SCOPE={bm.PIN_SCOPE} "
          f"DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} DS_PERTURB_UI_DOWN={bm.PERTURB_UI_DOWN:g} "
          f"readout=mean(unpinned) cutoffs {bm.DOWN_CUTOFF}/{bm.UP_CUTOFF}; measured |log2| >= "
          f"{THRESHOLD:.3f}; build {build.name}; code {sha}; dataset korkut2015 registered {role}; "
          f"{'PREDICT ONLY (values unread)' if predict_only else 'scoring'}", flush=True)
    rows, solves, conv = [], 0, 0
    pathways = sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists())
    for pdir in pathways:
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = {name: list(dict.fromkeys(u for s in st for u in stid_to_uuids.get(s, [])))
                    for name, st in r_stids.items()}
        ro_nodes = {k: v for k, v in ro_nodes.items() if v}
        if not ro_nodes:
            continue
        parsed = None
        for ds in drug_sets:
            genes = sorted({g for c in ds for g in drugs[c]})
            g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
            if not any(g2u.values()):
                continue
            pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
            pinned = set().union(*pins.values()) if pins else set()
            if not pinned:
                continue
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                adj = bm.build_adjacency(pdir)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed: {parsed.get('message')}")
            obs = {u: [bm.PERTURB_UI_DOWN, bm.PIN_CONFIDENCE] for u in sorted(pinned)}
            res = bm.solve_via_ds_api(None, obs, network_id=parsed["network_id"])
            if res.get("status") != "success":
                raise SystemExit(f"{pdir.name} {ds}: solve failed: {res.get('message')}")
            bm.refuse_unless_default_model(res, pdir.name)
            solves += 1
            conv += bool(res.get("converged"))
            act = res.get("node_activities", {})
            reach = bm.reachable_from(adj, pinned)
            for name, nodes in sorted(ro_nodes.items()):
                free = [u for u in nodes if u not in pinned]
                if not free:
                    status, ui, pred, rp = "trivial", None, None, None
                else:
                    ui = 100.0 * sum(act.get(u, 0.01) for u in free) / len(free)
                    status, pred = "scored", bm.classify(ui)
                    # specs/021 reachability baseline: DOWN when a readout node is reachable
                    rp = bm.DOWN if any(u in reach for u in free) else bm.NORMAL
                rows.append({"pathway_id": pdir.name, "drug_set": "+".join(ds), "readout": name,
                             "targets_pinned": ";".join(sorted(g for g, us in pins.items() if us)),
                             "pinned": len(pinned), "readout_nodes": len(nodes), "status": status,
                             "converged": int(bool(res.get("converged"))),
                             "pred_ui": "" if ui is None else f"{ui:.6g}",
                             "predicted": NAMES.get(pred, ""),
                             "reach_predicted": NAMES.get(rp, "")})
        print(f"  {pdir.name}: {len(ro_nodes)} readouts", flush=True)
    after = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if after.get("build_id") != build.name:
        raise SystemExit(f"the API switched to {after.get('build_id')!r} during the run; nothing recorded")
    if not predict_only:
        for thr, col in ((THRESHOLD, "measured"), (THRESHOLD_SECONDARY, "measured_secondary")):
            m = measured_calls(thr, design)
            for r in rows:
                r[col] = m[(tuple(r["drug_set"].split("+")), r["readout"])]
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "korkut_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Converged: {conv} of {solves} solves")
    print(f"Cases: {len(rows)} ({sum(r['status'] == 'scored' for r in rows)} scored) in "
          f"{len({r['pathway_id'] for r in rows})} pathways, {len(drug_sets)} drug sets")


if __name__ == "__main__":
    main()
