#!/usr/bin/env python3
"""specs/131 S2: collect the korkut solve problems (network, pins, readout nodes) for fitting the surrogate.

The pins, readout nodes and held nodes are korkut_bench's own (same helpers, same order); nothing is solved here.
Measured calls are read from a korkut case file (korkut is FIT). Writes a pickle outside the repository.

Usage: collect_131.py KORKUT_CASES_TSV OUT.pickle   (DELTASIGNAL_BASE, DS_CATALOG_ROOT as for the bench)
"""
import csv
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import dataset_roles  # noqa: E402
import drug_activity as da  # noqa: E402
import korkut_bench as kb  # noqa: E402


def main():
    cases_tsv, out = sys.argv[1:]
    dataset_roles.require_role("korkut2015", "fit")
    measured = {(r["pathway_id"], r["drug_set"], r["readout"]): r["measured"]
                for r in csv.DictReader(open(cases_tsv), delimiter="\t") if r["status"] == "scored"}
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    drug_sets = sorted({ds for _, ds in kb.conditions()})
    target_genes = sorted({g for gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = kb.readout_stids(readouts, gene_to_stids, kb.neo4j_site_stids(sites))
    units = da.active_units() if da.needs_units() else {}
    actions = da.drug_actions("korkut", sorted(drugs))

    networks, solves = {}, []
    build = bm.CATALOG_ROOT.resolve()
    for pdir in sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists()):
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
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed: {parsed.get('message')}")
            action = da.gene_actions([(actions[c], drugs[c]) for c in ds])
            pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, gene_to_stids, stid_to_uuids, units, action)
            if not roots:
                continue
            if cut:
                raise SystemExit("edge redirects are not modelled by the surrogate")
            rows = []
            for name, nodes in sorted(ro_nodes.items()):
                free = [u for u in nodes if u not in obs]
                key = (pdir.name, "+".join(ds), name)
                if free and key in measured:
                    rows.append((name, free, measured[key]))
            if rows:
                networks.setdefault(pdir.name, {k: parsed[k] for k in
                                                ("nodes", "edges", "cofactor_stids", "drug_stids", "pathogen_stids")})
                solves.append({"pathway": pdir.name, "drug_set": "+".join(ds), "obs": obs, "rows": rows})
        print(f"  {pdir.name}: {sum(s['pathway'] == pdir.name for s in solves)} solves", flush=True)
    with open(out, "wb") as f:
        pickle.dump({"networks": networks, "solves": solves}, f)
    print(f"{len(solves)} solves, {sum(len(s['rows']) for s in solves)} rows, {len(networks)} networks")


if __name__ == "__main__":
    main()
