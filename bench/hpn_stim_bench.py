#!/usr/bin/env python3
"""specs/184: predictions for the HPN-DREAM ligand stimulations (Hill et al. 2016), registered in
specs/184-ligand-stimulation/spec.md. The Lun protocol (bench/lun_bench.py) with ligands in place of proteins.

This script never opens the dataset: it only writes predictions. `bench/analysis/184/stim_184.py` scores them.

- A condition pins the root forms of its ligand(s) at the bench's over-expression level; one solve per condition per
  network that holds a root form. The withdrawal contrast (reported only) pins the same roots at 0, and the stimulation
  call is the opposite of that solve's class.
- A readout (specs/112's site map) is the mean over the unpinned nodes that are, or contain, an EWAS phosphorylated at
  one of its sites, classified with the curator cutoffs. Reachability: UP when a readout node is reachable from the
  pinned roots.

Usage: DS_CATALOG_ROOT=<build> hpn_stim_bench.py OUT_DIR
"""
import csv
import json
import os
import re
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import dataset_roles  # noqa: E402
from decryptm_bench import git_state  # noqa: E402
from korkut_bench import NAMES, neo4j_site_stids, read_tsv, readout_stids  # noqa: E402

SITE_MAP = Path(__file__).resolve().parents[1] / "specs" / "112-hpn-dream" / "site_map.tsv"
CONDITIONS = {"EGF": ["EGF"], "FGF1": ["FGF1"], "HGF": ["HGF"], "IGF1": ["IGF1"], "Insulin": ["INS"],
              "NRG1": ["NRG1"], "Serum": ["EGF", "FGF1", "HGF", "IGF1", "INS", "NRG1"]}
OPPOSITE = {bm.UP: bm.DOWN, bm.DOWN: bm.UP, bm.NORMAL: bm.NORMAL}


def readouts(path=SITE_MAP):
    """specs/112's site map in the korkut readout format; readouts without a site (a total protein) are left out."""
    out = []
    for r in read_tsv(path):
        sites = [p for p in r["sites"].split(";") if p]
        if sites:
            out.append({"readout": r["species"], "kind": "phospho",
                        "sites": ",".join(f"{g}:{int(re.sub(r'^[A-Z]', '', s))}" for g, s in (p.split(":") for p in sites))})
    return out


def condition_roots(pins, genes):
    return sorted({u for g in genes for u in pins.get(g) or []})


def main():
    out_dir = Path(sys.argv[1])
    role = dataset_roles.require_role("hpn_dream", "score")   # specs/051; predictions read no value
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")

    ros = readouts()
    genes = sorted({g for gs in CONDITIONS.values() for g in gs})
    gene_to_stids = bm.neo4j_gene_to_stids(genes)
    sites = sorted({(s.split(":")[0], int(s.split(":")[1])) for r in ros for s in r["sites"].split(",")})
    r_stids = readout_stids(ros, {}, neo4j_site_stids(sites))
    up = bm.PERTURB_FOLD_UP / bm.CEILING * 100.0

    print(f"Protocol: HPN-DREAM ligand stimulation bench (specs/184), PREDICT ONLY (the dataset is never opened); "
          f"DS_PIN_SCOPE={bm.PIN_SCOPE} DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} "
          f"DS_BASELINE={bm.BASELINE:g} stimulation pin {bm.PERTURB_FOLD_UP:g}x (UI {up:g}), withdrawal pin "
          f"{bm.PERTURB_UI_DOWN:g}; readout=mean(unpinned) cutoffs {bm.DOWN_CUTOFF}/{bm.UP_CUTOFF}; build {build.name}; "
          f"code {sha}; dataset hpn_dream registered {role}; ligands {len(genes)} "
          f"({sum(1 for g in genes if gene_to_stids.get(g))} in Reactome); readouts {len(ros)} "
          f"({sum(1 for v in r_stids.values() if v)} with a phospho form in Reactome)", flush=True)
    rows, solves, conv = [], 0, 0
    for pdir in sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists()):
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = {name: list(dict.fromkeys(u for s in st for u in stid_to_uuids.get(s, [])))
                    for name, st in r_stids.items()}
        ro_nodes = {k: v for k, v in ro_nodes.items() if v}
        if not ro_nodes:
            continue
        g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
        pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
        parsed = adj = None
        for cond, cgenes in CONDITIONS.items():
            roots = condition_roots(pins, cgenes)
            if not roots:
                continue
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed: {parsed.get('message')}")
                adj = bm.build_adjacency(pdir)
            acts, convs = {}, {}
            for arm, level in (("stim", up), ("withdraw", bm.PERTURB_UI_DOWN)):
                res = bm.solve_via_ds_api(None, {u: [level, bm.PIN_CONFIDENCE] for u in roots},
                                          network_id=parsed["network_id"])
                if res.get("status") != "success":
                    raise SystemExit(f"{pdir.name} {cond} {arm}: solve failed: {res.get('message')}")
                bm.refuse_unless_default_model(res, pdir.name)
                served_base = float(res.get("resting_baseline", 0.01))
                if abs(served_base - bm.BASELINE) > 1e-12 * bm.BASELINE:
                    raise SystemExit(f"the server rests nodes at {served_base}, the bench expects {bm.BASELINE}")
                solves += 1
                conv += bool(res.get("converged"))
                acts[arm], convs[arm] = res.get("node_activities", {}), int(bool(res.get("converged")))
            reach, pinned = bm.reachable_from(adj, set(roots)), set(roots)
            for name, nodes in sorted(ro_nodes.items()):
                free = [u for u in nodes if u not in pinned]
                row = {"pathway_id": pdir.name, "condition": cond, "readout": name,
                       "ligands_pinned": sum(1 for g in cgenes if pins.get(g)), "pinned": len(pinned),
                       "readout_nodes": len(nodes), "status": "scored" if free else "trivial",
                       "converged": convs["stim"], "converged_withdraw": convs["withdraw"],
                       "pred_ui": "", "predicted": "", "reach_predicted": "", "pred_ui_withdraw": "",
                       "predicted_withdraw": ""}
                if free:
                    ui = sum(acts["stim"].get(u, bm.BASELINE) for u in free) / len(free) * bm.CEILING
                    uw = sum(acts["withdraw"].get(u, bm.BASELINE) for u in free) / len(free) * bm.CEILING
                    row.update(pred_ui=f"{ui:.6g}", predicted=NAMES[bm.classify(ui)],
                               reach_predicted=NAMES[bm.UP if any(u in reach for u in free) else bm.NORMAL],
                               pred_ui_withdraw=f"{uw:.6g}", predicted_withdraw=NAMES[OPPOSITE[bm.classify(uw)]])
                rows.append(row)
        print(f"  {pdir.name}: {len(ro_nodes)} readouts, "
              f"{sum(1 for c, g in CONDITIONS.items() if condition_roots(pins, g))} conditions", flush=True)
    after = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if after.get("build_id") != build.name:
        raise SystemExit(f"the API switched to {after.get('build_id')!r} during the run; nothing recorded")
    if not rows:
        raise SystemExit("no network holds a ligand root and a readout; nothing recorded")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "hpn_stim_predict_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Converged: {conv} of {solves} solves")
    print(f"Cases: {len(rows)} ({sum(r['status'] == 'scored' for r in rows)} scored) in "
          f"{len({r['pathway_id'] for r in rows})} pathways, {len({r['condition'] for r in rows})} conditions")


if __name__ == "__main__":
    main()
