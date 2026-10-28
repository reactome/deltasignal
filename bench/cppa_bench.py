#!/usr/bin/env python3
"""specs/121: the frozen model's predictions for CPPA single drugs (values-blind; no CPPA measurement is opened).

As korkut_bench.py, one drug at a time: the drug's target genes (from cppa_targets_121.py: gene symbols for class and
ligand rows, UniProt accessions otherwise, translated by Reactome 97's accession-to-gene table) are pinned to
DS_PERTURB_UI_DOWN at the forms the curator bench pins, one solve per drug per pathway. A readout (the CPPA map, the
korkut format) is the mean over the nodes that are, or contain, its EWAS, pinned nodes left out. Stimuli are not
solved. The case log has no measured column: scoring is a separate, registered step.

Usage: DS_CATALOG_ROOT=<build> cppa_bench.py DRUG_TARGETS_TSV ACC_GENE_TSV READOUTS_TSV OUT_DIR
"""
import csv
import json
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
from decryptm_bench import git_state  # noqa: E402
from korkut_bench import NAMES, neo4j_site_stids, read_tsv, readout_stids  # noqa: E402


def drug_genes(drug_targets_tsv, acc_gene_tsv):
    """{drug: (role, [genes])} for every drug with a target Reactome knows; stimuli left out."""
    acc_gene = {r["acc"]: r["gene"] for r in csv.DictReader(open(acc_gene_tsv), delimiter="\t")}
    out = {}
    for r in csv.DictReader(open(drug_targets_tsv), delimiter="\t"):
        if r["role"] == "stimulus":
            continue
        gene = r["target"] if r["source"] in ("class", "ligand") else acc_gene.get(r["target"])
        if gene:
            out.setdefault(r["drug"], (r["role"], set()))[1].add(gene)
    return {d: (role, sorted(gs)) for d, (role, gs) in sorted(out.items())}


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    drug_targets, acc_gene, readouts_tsv, out_dir = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")

    drugs = drug_genes(drug_targets, acc_gene)
    readouts = [r for r in read_tsv(readouts_tsv) if r["kind"] != "exclude"]
    target_genes = sorted({g for _, gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = readout_stids(readouts, gene_to_stids, neo4j_site_stids(sites))
    print(f"Protocol: CPPA single-drug predictions (specs/121), DS_PIN_SCOPE={bm.PIN_SCOPE} "
          f"DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} DS_PERTURB_UI_DOWN={bm.PERTURB_UI_DOWN:g} "
          f"readout=mean(unpinned) cutoffs {bm.DOWN_CUTOFF}/{bm.UP_CUTOFF}; build {build.name}; code {sha}; "
          f"{len(drugs)} drugs, {len(readouts)} readouts; PREDICT ONLY (values unread)", flush=True)
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
        for drug, (role, genes) in drugs.items():
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
                raise SystemExit(f"{pdir.name} {drug}: solve failed: {res.get('message')}")
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
                    rp = bm.DOWN if any(u in reach for u in free) else bm.NORMAL
                rows.append({"pathway_id": pdir.name, "drug": drug, "role": role, "readout": name,
                             "targets_pinned": ";".join(sorted(g for g, us in pins.items() if us)),
                             "pinned": len(pinned), "readout_nodes": len(nodes), "status": status,
                             "converged": int(bool(res.get("converged"))),
                             "pred_ui": "" if ui is None else f"{ui:.6g}",
                             "predicted": NAMES.get(pred, ""), "reach_predicted": NAMES.get(rp, "")})
        print(f"  {pdir.name}: {len(ro_nodes)} readouts", flush=True)
    after = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if after.get("build_id") != build.name:
        raise SystemExit(f"the API switched to {after.get('build_id')!r} during the run; nothing recorded")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "cppa_predict_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Converged: {conv} of {solves} solves")
    print(f"Cases: {len(rows)} ({sum(r['status'] == 'scored' for r in rows)} scored) in "
          f"{len({r['pathway_id'] for r in rows})} pathways, {len({r['drug'] for r in rows})} drugs")


if __name__ == "__main__":
    main()
