#!/usr/bin/env python3
"""specs/127: predictions for Lun et al. 2019 (649 kinases and phosphatases over-expressed one at a time in HEK293T,
CyTOF readouts). The scoring is registered in specs/127-lun2019-scoring/spec.md.

This script never opens the dataset: it only writes predictions (the values-blind step). The measured calls are made by
`bench/lun_score.py`, once, at specs/125 stage 5.

- A protein of interest (bench/analysis/127/pois.tsv, Table S1) is resolved by its UniProt accession and pinned at the
  bench's over-expression level at the forms the curator bench pins (`pin_occurrences`), one solve per protein per
  pathway that holds a root form of it.
- A readout (bench/analysis/127/readouts.tsv, Table S2) is the mean over the nodes that are, or contain, an EWAS of the
  gene phosphorylated at the coordinate (phospho) or any EWAS of the gene (total), pinned nodes left out; classified
  with the curator cutoffs. Reachability: UP when a readout node is reachable from the pinned roots.
- A case row is (protein, pathway, readout); `self` marks a protein that is one of the readout's genes.

Usage: DS_CATALOG_ROOT=<build> lun_bench.py OUT_DIR
"""
import csv
import json
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import dataset_roles  # noqa: E402
from decryptm_bench import git_state  # noqa: E402
from korkut_bench import NAMES, neo4j_site_stids, read_tsv, readout_stids  # noqa: E402

MAPS = Path(__file__).resolve().parent / "analysis" / "127"


def pois(path=None):
    """[(gene, accession)] once per gene, in Table S1 order; entries without an accession are left out."""
    out, seen = [], set()
    for r in read_tsv(path or MAPS / "pois.tsv"):
        if r["accession"] and r["gene"] not in seen:
            seen.add(r["gene"])
            out.append((r["gene"], r["accession"]))
    return out


def readout_genes(readout):
    return {s.split(":")[0] for s in readout["sites"].split(",") if s}


def neo4j_accession_to_stids(accessions):
    """{accession: [stIds]}: every physical entity whose reference entity is that UniProt accession."""
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $accs AS acc
        MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE re.identifier = acc
        RETURN acc AS acc, COLLECT(DISTINCT pe.stId) AS stids
        """,
        accs=list(accessions),
    ).data()
    return {r["acc"]: r["stids"] for r in rows}


def main():
    out_dir = Path(sys.argv[1])
    role = dataset_roles.require_role("lun2019", "score")   # specs/051; predictions read no value
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")

    proteins = pois()
    readouts = [r for r in read_tsv(MAPS / "readouts.tsv") if r["kind"] in ("phospho", "total")]
    acc_to_stids = neo4j_accession_to_stids(sorted({a for _, a in proteins}))
    gene_to_stids = {g: acc_to_stids.get(a, []) for g, a in proteins}
    ro_genes = sorted({g for r in readouts for g in readout_genes(r)})
    ro_gene_stids = bm.neo4j_gene_to_stids(ro_genes)
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = readout_stids(readouts, ro_gene_stids, neo4j_site_stids(sites))
    selfs = {r["readout"]: readout_genes(r) for r in readouts}
    up = bm.PERTURB_FOLD_UP / bm.CEILING * 100.0

    print(f"Protocol: Lun 2019 over-expression bench (specs/127), PREDICT ONLY (the dataset is never opened); "
          f"DS_PIN_SCOPE={bm.PIN_SCOPE} DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} "
          f"DS_BASELINE={bm.BASELINE:g} over-expression pin {bm.PERTURB_FOLD_UP:g}x (UI {up:g}); readout=mean(unpinned) "
          f"cutoffs {bm.DOWN_CUTOFF}/{bm.UP_CUTOFF}; build {build.name}; code {sha}; dataset lun2019 registered {role}; "
          f"proteins {len(proteins)} ({sum(1 for g in gene_to_stids.values() if g)} in Reactome); "
          f"readouts {len(readouts)}", flush=True)
    rows, solves, conv = [], 0, 0
    pathways = sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists())
    for pdir in pathways:
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = {name: list(dict.fromkeys(u for s in st for u in stid_to_uuids.get(s, [])))
                    for name, st in r_stids.items()}
        ro_nodes = {k: v for k, v in ro_nodes.items() if v}
        if not ro_nodes:
            continue
        parsed = adj = None
        for gene, _ in proteins:
            g2u = {gene: [u for s in gene_to_stids[gene] for u in stid_to_uuids.get(s, [])]}
            if not g2u[gene]:
                continue
            pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
            roots = sorted(pins.get(gene) or [])
            if not roots:
                continue
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed: {parsed.get('message')}")
                adj = bm.build_adjacency(pdir)
            obs = {u: [up, bm.PIN_CONFIDENCE] for u in roots}
            res = bm.solve_via_ds_api(None, obs, network_id=parsed["network_id"])
            if res.get("status") != "success":
                raise SystemExit(f"{pdir.name} {gene}: solve failed: {res.get('message')}")
            bm.refuse_unless_default_model(res, pdir.name)
            served_base = float(res.get("resting_baseline", 0.01))
            if abs(served_base - bm.BASELINE) > 1e-12 * bm.BASELINE:
                raise SystemExit(f"the server rests nodes at {served_base}, the bench expects {bm.BASELINE}")
            solves += 1
            conv += bool(res.get("converged"))
            act = res.get("node_activities", {})
            reach = bm.reachable_from(adj, set(roots))
            pinned = set(roots)
            for name, nodes in sorted(ro_nodes.items()):
                free = [u for u in nodes if u not in pinned]
                if not free:
                    status, ui, pred, rp = "trivial", None, None, None
                else:
                    ui = sum(act.get(u, bm.BASELINE) for u in free) / len(free) * bm.CEILING   # fold, as pred_ui
                    status, pred = "scored", bm.classify(ui)
                    rp = bm.UP if any(u in reach for u in free) else bm.NORMAL
                rows.append({"pathway_id": pdir.name, "gene": gene, "readout": name,
                             "self": int(gene in selfs[name]), "pinned": len(pinned), "readout_nodes": len(nodes),
                             "status": status, "converged": int(bool(res.get("converged"))),
                             "pred_ui": "" if ui is None else f"{ui:.6g}",
                             "predicted": NAMES.get(pred, ""), "reach_predicted": NAMES.get(rp, "")})
        print(f"  {pdir.name}: {len(ro_nodes)} readouts", flush=True)
    after = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if after.get("build_id") != build.name:
        raise SystemExit(f"the API switched to {after.get('build_id')!r} during the run; nothing recorded")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "lun_predict_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Converged: {conv} of {solves} solves")
    print(f"Cases: {len(rows)} ({sum(r['status'] == 'scored' for r in rows)} scored) in "
          f"{len({r['pathway_id'] for r in rows})} pathways, {len({r['gene'] for r in rows})} proteins")


if __name__ == "__main__":
    main()
