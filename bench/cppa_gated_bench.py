#!/usr/bin/env python3
"""specs/121 phase 2: predictions for CPPA (drug, cell line) pairs with expression gating (values-blind).

As cppa_bench.py, plus, for each cell line, every node it cannot form (expression_gate.gated_nodes: a required
protein's gene below log2(TPM+1) = 1 in DepMap 24Q2) pinned at baseline, UI 1; the drug's target pins win where they
overlap. The pairs are read from the CPPA design (drug and cell line only) for one portion of the split; no
measurement is opened. A cell line maps to DepMap by exact name, else its longest-prefix parent (>= 4 characters);
a cell line with no expression is not solved here (its samples keep the frozen prediction). CPPA_SHARD as cppa_bench.

Usage: DS_CATALOG_ROOT=<build> cppa_gated_bench.py PORTION CPPA_UNPACKED SPLIT_TSV MODEL_CSV EXPRESSION_CSV
                                                    DRUG_TARGETS_TSV ACC_GENE_TSV READOUTS_TSV OUT_DIR
"""
import csv
import json
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import cppa_data_121  # noqa: E402
import expression_gate  # noqa: E402
import split_121  # noqa: E402
from cppa_bench import drug_genes  # noqa: E402
from decryptm_bench import git_state  # noqa: E402
from korkut_bench import NAMES, neo4j_site_stids, read_tsv, readout_stids  # noqa: E402

THRESHOLD = 1.0      # log2(TPM + 1): TPM below 1 is absent (specs/121 amendment 2)


def depmap_models(model_csv):
    """{normalised name: ModelID}; only the identifier and name columns are read (the table has patient columns)."""
    out = {}
    with open(model_csv, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        i_id, i_strip, i_name = (header.index(c) for c in ("ModelID", "StrippedCellLineName", "CellLineName"))
        for r in rd:
            for i in (i_strip, i_name):
                out.setdefault(split_121.norm(r[i]), r[i_id])
    return out


def model_for(cell, models):
    if cell in models:
        return models[cell], "exact"
    parents = [p for p in models if len(p) >= 4 and cell.startswith(p)]
    return (models[max(parents, key=len)], "parent") if parents else (None, "none")


def expression_rows(expression_csv, model_ids, genes):
    """{ModelID: {gene: log2(TPM+1)}} for the models and genes asked for, one pass over the table."""
    out = {}
    with open(expression_csv) as f:
        header = f.readline().rstrip("\n").split(",")
        col = {h.split(" (")[0]: i for i, h in enumerate(header) if i and h.split(" (")[0] in genes}
        for line in f:
            mid = line.split(",", 1)[0]
            if mid in model_ids:
                vals = line.rstrip("\n").split(",")
                out[mid] = {g: float(vals[i]) for g, i in col.items() if vals[i] != ""}
    return out


def main():
    if len(sys.argv) != 10:
        raise SystemExit(__doc__)
    portion, unpacked, split_tsv, model_csv, expr_csv, drug_targets, acc_gene, readouts_tsv, out_dir = sys.argv[1:]
    out_dir = Path(out_dir)
    sha, dirty = git_state(Path(__file__).resolve().parents[1])
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")

    drugs = drug_genes(drug_targets, acc_gene)
    roles = cppa_data_121.roles(split_tsv)
    pairs = sorted({(d, c) for d, c, _, _ in cppa_data_121.design(unpacked).values()
                    if d in drugs and cppa_data_121.portion_of(roles["drug"][d], roles["cell_line"][c]) == portion})
    models = depmap_models(model_csv)
    cell_model = {c: model_for(c, models) for c in sorted({c for _, c in pairs})}
    pairs = [(d, c) for d, c in pairs if cell_model[c][0]]

    readouts = [r for r in read_tsv(readouts_tsv) if r["kind"] != "exclude"]
    all_genes = {r["gene"] for r in csv.DictReader(open(acc_gene), delimiter="\t")}
    expr = expression_rows(expr_csv, {m for m, _ in cell_model.values() if m}, all_genes)
    # targets and readouts resolve exactly as in cppa_bench (the frozen arm); absent genes use the one-pass table
    target_genes = sorted({g for _, gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    gene_ewas = expression_gate.human_gene_ewas()
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = readout_stids(readouts, gene_to_stids, neo4j_site_stids(sites))
    absent = {}
    for c, (m, _) in cell_model.items():
        if m in expr:
            gs = {g for g, v in expr[m].items() if v < THRESHOLD}
            absent[c] = set().union(*(set(gene_ewas.get(g, [])) for g in gs)) if gs else set()
    pairs = [(d, c) for d, c in pairs if c in absent]
    print(f"Protocol: CPPA gated predictions (specs/121 amendment 2), portion {portion}, log2(TPM+1) < {THRESHOLD:g} "
          f"absent, DS_PIN_SCOPE={bm.PIN_SCOPE} DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} "
          f"DS_PERTURB_UI_DOWN={bm.PERTURB_UI_DOWN:g}; build {build.name}; code {sha}; {len(pairs)} pairs, "
          f"{len(absent)} cell lines with expression; PREDICT ONLY (values unread)", flush=True)

    rows, solves, conv = [], 0, 0
    pathways = sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists())
    shard = os.environ.get("CPPA_SHARD", "")
    if shard:
        i, n = map(int, shard.split("/"))
        pathways = pathways[i::n]
    for pdir in pathways:
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = {name: list(dict.fromkeys(u for s in st for u in stid_to_uuids.get(s, [])))
                    for name, st in r_stids.items()}
        ro_nodes = {k: v for k, v in ro_nodes.items() if v}
        if not ro_nodes:
            continue
        with open(pdir / "nodes.csv") as f:
            nrows = list(csv.DictReader(f))
        below = expression_gate.entity_leaves({x for r in nrows if r["node_kind"] == "set_variant"
                                               for col in ("source_sets", "chosen_members")
                                               for x in (r.get(col) or "").split("|") if x})
        gated = {}
        parsed = adj = None
        for drug, cell in pairs:
            genes = drugs[drug][1]
            g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
            if not any(g2u.values()):
                continue
            pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
            pinned = set().union(*pins.values()) if pins else set()
            if not pinned:
                continue
            if cell not in gated:
                gated[cell] = expression_gate.gated_nodes(pdir, absent[cell], below)
            held = gated[cell] - pinned
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                adj = bm.build_adjacency(pdir)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed: {parsed.get('message')}")
            obs = {u: [bm.PERTURB_UI_DOWN, bm.PIN_CONFIDENCE] for u in sorted(pinned)}
            obs.update({u: [1.0, bm.PIN_CONFIDENCE] for u in sorted(held)})
            res = bm.solve_via_ds_api(None, obs, network_id=parsed["network_id"])
            if res.get("status") != "success":
                raise SystemExit(f"{pdir.name} {drug} {cell}: solve failed: {res.get('message')}")
            bm.refuse_unless_default_model(res, pdir.name)
            solves += 1
            conv += bool(res.get("converged"))
            act = res.get("node_activities", {})
            reach = bm.reachable_from(adj, pinned)
            for name, nodes in sorted(ro_nodes.items()):
                free = [u for u in nodes if u not in pinned and u not in held]
                if not free:
                    status, ui, pred, rp = "trivial", None, None, None
                else:
                    ui = 100.0 * sum(act.get(u, 0.01) for u in free) / len(free)
                    status, pred = "scored", bm.classify(ui)
                    rp = bm.DOWN if any(u in reach for u in free) else bm.NORMAL
                rows.append({"pathway_id": pdir.name, "drug": drug, "cell_line": cell,
                             "model_match": cell_model[cell][1], "readout": name, "pinned": len(pinned),
                             "held": len(held), "readout_nodes": len(nodes), "status": status,
                             "converged": int(bool(res.get("converged"))),
                             "pred_ui": "" if ui is None else f"{ui:.6g}",
                             "predicted": NAMES.get(pred, ""), "reach_predicted": NAMES.get(rp, "")})
        print(f"  {pdir.name}: {len(ro_nodes)} readouts", flush=True)
    after = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if after.get("build_id") != build.name:
        raise SystemExit(f"the API switched to {after.get('build_id')!r} during the run; nothing recorded")
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"cppa_gated_{portion}" + (f".{shard.replace('/', 'of')}" if shard else "") + ".tsv"
    with open(out_dir / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Converged: {conv} of {solves} solves")
    print(f"Cases: {len(rows)} ({sum(r['status'] == 'scored' for r in rows)} scored), "
          f"{len({(r['drug'], r['cell_line']) for r in rows})} pairs")


if __name__ == "__main__":
    main()
