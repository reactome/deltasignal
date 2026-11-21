#!/usr/bin/env python3
"""specs/161: the korkut replay with engagement-scaled pins.

Identical to the specs/157 canonical run's korkut leg (same pin resolution,
the same /api/readouts valuation) except each drug's pins are set to its
measured engagement (engagement_161.py) instead of 0. Korkut phospho readouts
are valued in `sites` mode (issue #278, korkut_bench's valuation);
`--phospho-mode forms` reproduces the 2026-11-16 arm. The env is verified by
the specs/157 positive probe (every registry rule canonical) before any
measurement.

Usage: DELTASIGNAL_BASE=... DS_CATALOG_ROOT=<served build> \
  arm_161.py REGISTRY_TSV KORKUT_CASES_TSV OUT_TSV [--phospho-mode sites|forms]
"""
import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "157"))
sys.path.insert(0, str(HERE.parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import drug_activity as da  # noqa: E402
import dataset_roles  # noqa: E402
import korkut_bench as kb  # noqa: E402
from ensemble_157 import korkut_query_mode, load_registry, post, verify_env  # noqa: E402
from engagement_161 import engagement  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("registry_tsv"); ap.add_argument("cases_tsv"); ap.add_argument("out_tsv")
    ap.add_argument("--phospho-mode", default="sites", choices=["sites", "forms"])
    a = ap.parse_args()
    registry_tsv, cases_tsv, out_tsv, phospho_mode = a.registry_tsv, a.cases_tsv, a.out_tsv, a.phospho_mode
    print(f"korkut phospho mode {phospho_mode}", flush=True)
    dataset_roles.require_role("korkut2015", "fit")
    verify_env(load_registry(registry_tsv), None)
    eng = {d: e for d, (e, _) in engagement().items()}

    recorded = set()
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] == "scored":
            recorded.add((r["pathway_id"], r["drug_set"], r["readout"]))
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    gene_drug = {g: d for d, gs in drugs.items() for g in gs}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    tg = sorted({g for gs in drugs.values() for g in gs})
    rg = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(sorted(set(tg) | set(rg)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    s2s = kb.neo4j_site_stids(sites)
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, g2s), s2s)
    kinds = {r["readout"]: r["kind"] for r in readouts}
    units = da.active_units() if da.needs_units() else {}
    actions = da.drug_actions("korkut", sorted(drugs))
    drug_sets = sorted({ds for _, ds in kb.conditions()})
    pathways = sorted(p for p in bm.CATALOG_ROOT.resolve().iterdir()
                      if (p / "stid_to_uuid_mapping.csv").exists())

    solves = 0
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["axis", "pathway", "case_a", "case_b", "readout", "fold", "call", "converged"])
        for pdir in pathways:
            s2u = bm.load_stid_to_uuids(pdir)
            parsed = None
            for ds in drug_sets:
                ds_name = "+".join(ds)
                names = sorted(n for n in r_stids if (pdir.name, ds_name, n) in recorded)
                if not names:
                    continue
                genes = sorted({g for c in ds for g in drugs[c]})
                if parsed is None:
                    parsed = bm.parse_pathway_via_ds_api(pdir.name)
                    if parsed.get("status") != "success":
                        raise SystemExit(f"{pdir.name}: parse failed")
                action = da.gene_actions([(actions[c], drugs[c]) for c in ds])
                pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, g2s, s2u, units, action)
                if not roots:
                    continue
                # THE ONLY CHANGE: each pinned uuid at its drug's engagement
                # (a uuid shared by two drugs' targets takes the smaller e).
                pin_value = {}
                for g, us in pins.items():
                    e = eng[gene_drug[g]]
                    for u in us:
                        pin_value[u] = min(pin_value.get(u, 1.0), e)
                obs = {u: [pin_value[u], v[1]] for u, v in obs.items()}
                payload = {"network_id": parsed["network_id"], "observations": obs}
                if cut:
                    payload["edge_redirects"] = [list(e) for e in sorted(cut)]
                sol = post("solve", payload)
                if sol.get("status") != "success":
                    raise SystemExit(f"{pdir.name} {ds_name}: solve failed")
                solves += 1
                queries = [{"id": n, "stable_ids": sorted(r_stids[n]),
                            "mode": korkut_query_mode(kinds.get(n), phospho_mode)}
                           for n in names]
                ro = post("readouts", {"network_id": parsed["network_id"],
                                       "activities": sol.get("node_activities") or {},
                                       "pinned": sorted(obs), "queries": queries,
                                       "converged": bool(sol.get("converged", False))})
                down, up = ro["protocol"]["down_cutoff"], ro["protocol"]["up_cutoff"]
                for e in ro["results"]:
                    fold = e.get("fold")
                    call = ("" if fold is None else
                            "down" if fold < down else "up" if fold >= up else "not")
                    w.writerow(["korkut", pdir.name, ds_name, "", e["id"], fold, call,
                                int(bool(sol.get("converged")))])
            print(f"  {pdir.name}: done ({solves} solves)", flush=True)
    print(f"{solves} solves -> {out_tsv}")


if __name__ == "__main__":
    main()
