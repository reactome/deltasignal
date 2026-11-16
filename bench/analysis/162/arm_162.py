#!/usr/bin/env python3
"""specs/162: CPPA solves for both arms: full-knockout pins (baseline, one
solve per (drug, pathway)) and engagement pins (one solve per (drug, cell line,
pathway), pins at e from engagement_162.py). The valuation is cppa_bench.py's
own, applied identically to both arms: a readout is the mean activity over
its unpinned nodes (specs/134 total-readout nodes). Values-blind: no CPPA
measurement is read here. Output stays OUTSIDE the repository (CPPA terms).

CPPA_SHARD=i/n runs every n-th pathway from the i-th.

Usage: DELTASIGNAL_BASE=... DS_CATALOG_ROOT=<served build> arm_162.py \
  REGISTRY_TSV ENGAGEMENT_TSV DRUG_TARGETS_TSV ACC_GENE_TSV READOUTS_TSV OUT_TSV
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "157"))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import drug_activity as da  # noqa: E402
from cppa_bench import drug_genes  # noqa: E402
from ensemble_157 import load_registry, verify_env  # noqa: E402
from korkut_bench import (neo4j_site_stids, read_tsv, readout_nodes, readout_stids,  # noqa: E402
                          total_gene_stids, total_readout_nodes)


def main():
    registry, engagement_tsv, drug_targets, acc_gene, readouts_tsv, out_tsv = sys.argv[1:7]
    verify_env(load_registry(registry), None)
    eng = defaultdict(dict)
    for r in csv.DictReader(open(engagement_tsv), delimiter="\t"):
        eng[r["drug"]][r["cell_line"]] = float(r["engagement"])
    drugs = {d: v for d, v in drug_genes(drug_targets, acc_gene).items() if d in eng}
    actions = da.drug_actions("cppa", sorted(drugs))
    readouts = [r for r in read_tsv(readouts_tsv) if r["kind"] != "exclude"]
    tg = sorted({g for _, gs in drugs.values() for g in gs})
    rg = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(sorted(set(tg) | set(rg)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = readout_stids(readouts, total_gene_stids(readouts, g2s), neo4j_site_stids(sites))
    totals = {r["readout"] for r in readouts if r["kind"] == "total"}
    units = da.active_units() if da.needs_units() else {}
    pathways = sorted(p for p in bm.CATALOG_ROOT.resolve().iterdir()
                      if (p / "stid_to_uuid_mapping.csv").exists())
    shard = os.environ.get("CPPA_SHARD", "")
    if shard:
        i, n = map(int, shard.split("/"))
        pathways = pathways[i::n]
    solves = 0
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["arm", "pathway_id", "drug", "cell_line", "readout", "pred_ui", "converged"])
        for pdir in pathways:
            s2u = bm.load_stid_to_uuids(pdir)
            ro_nodes = readout_nodes(r_stids, s2u)
            if not ro_nodes:
                continue
            ro_nodes = total_readout_nodes(ro_nodes, totals, pdir)
            parsed = None
            for drug, (_, genes) in drugs.items():
                g2u = {g: [u for s in g2s.get(g, []) for u in s2u.get(s, [])] for g in genes}
                if not any(g2u.values()):
                    continue
                if parsed is None:
                    parsed = bm.parse_pathway_via_ds_api(pdir.name)
                    if parsed.get("status") != "success":
                        raise SystemExit(f"{pdir.name}: parse failed")
                pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, g2s, s2u, units,
                                                   actions[drug])
                if not roots:
                    continue
                pinned = set(obs)
                arms = [("baseline", "", obs)] + [
                    ("engagement", line, {u: [e, v[1]] for u, v in obs.items()})
                    for line, e in sorted(eng[drug].items())]
                for arm, line, o in arms:
                    res = bm.solve_via_ds_api(None, o, network_id=parsed["network_id"],
                                              edge_redirects=cut)
                    if res.get("status") != "success":
                        raise SystemExit(f"{pdir.name} {drug} {line}: solve failed")
                    bm.refuse_unless_default_model(res, pdir.name)
                    solves += 1
                    act = res.get("node_activities", {})
                    for name, nodes in sorted(ro_nodes.items()):
                        free = [u for u in nodes if u not in pinned]
                        if not free:
                            continue
                        ui = 100.0 * sum(act.get(u, 0.01) for u in free) / len(free)
                        w.writerow([arm, pdir.name, drug, line, name, f"{ui:.6g}",
                                    int(bool(res.get("converged")))])
            print(f"  {pdir.name}: done ({solves} solves)", flush=True)
    print(f"{solves} solves -> {out_tsv}")


if __name__ == "__main__":
    main()
