#!/usr/bin/env python3
"""specs/156 FR-004: the registered arm — does solver-level conditioning help?

Leave-one-readout-out on korkut (FIT): for each (pathway, drug set, readout)
the canonical bench scored, pin the drug targets hard EXACTLY as the recorded
protocol does (korkut_bench's own perturb machinery), soft-observe the drug
set's OTHER measured readouts at their measured folds, and predict the
held-out readout. Weights {0.25, 0.5, 0.75}, registered, no fitting. The gate
(spec.md FR-004): balanced accuracy on held-out readouts beats the
unconditioned record with a 95% bootstrap interval over DRUG SETS excluding
zero, at some registered weight.

Protocol details registered here, before the run:
- A conditioning readout's soft value is its measured fold, 2^mean(log2) over
  the drug set's conditions (the same expr.csv rows the measured call reads).
- NO LEAKAGE: a conditioning node that is also one of the held-out readout's
  nodes is dropped from the soft set.
- Two conditioning readouts sharing a node: the node takes 2^mean of their
  log2 folds (one value; the API refuses contradictions by design).
- Pinned nodes are never soft-observed (the API reports them overridden; the
  arm simply does not send them).
- The baseline is the recorded canonical scoring's `predicted` on exactly the
  keys the arm scores; measured truth is the recorded primary-threshold call.

Usage:
  DELTASIGNAL_BASE=http://localhost:PORT arm_156.py CANONICAL_CASES_TSV OUT_DIR
      [--pathways A,B,...]   # shard for parallel servers
The server must run DS_SOFT_OBS=blend; the script VERIFIES it positively (a
soft-carrying probe must succeed) before any measurement, and refuses a
server whose catalog is not the serving record.
"""
import collections
import csv
import json
import math
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import drug_activity as da  # noqa: E402
import dataset_roles  # noqa: E402
import korkut_bench as kb  # noqa: E402

WEIGHTS = (0.25, 0.5, 0.75)


def post(ep, body):
    rq = urllib.request.Request(f"{bm.DS_BASE}/api/{ep}", json.dumps(body).encode(),
                                {"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(rq, timeout=600).read())


def main():
    args = [a for a in sys.argv[1:]]
    shard = None
    if "--pathways" in args:
        i = args.index("--pathways")
        shard = set(args[i + 1].split(","))
        del args[i:i + 2]
    cases_tsv, out_dir = Path(args[0]), Path(args[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset_roles.require_role("korkut2015", "fit")

    served = (json.loads(urllib.request.urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read())
              .get("catalog") or {})
    build = served.get("build_id")
    print(f"server {bm.DS_BASE} serving {build}", flush=True)

    # The recorded canonical scoring: baseline calls and the measured truth.
    recorded = {}
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored" or not r["pred_ui"]:
            continue
        recorded[(r["pathway_id"], r["drug_set"], r["readout"])] = r

    # Measured folds per (drug set tuple, readout): 2^mean(log2) over conditions.
    names = [l.strip() for l in open(kb.DATA / "node_Index.csv")]
    design = kb.conditions()
    rows = list(csv.reader(open(kb.DATA / "expr.csv")))
    if len(rows) != len(design) or any(len(row) != len(names) for row in rows):
        raise SystemExit("expr.csv misaligned (specs/062 guard)")
    vals = collections.defaultdict(list)
    for (_, ds), row in zip(design, rows):
        for name, v in zip(names, row):
            vals[(ds, name)].append(float(v))
    meas_fold = {k: 2.0 ** (sum(v) / len(v)) for k, v in vals.items()}

    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    drug_sets = sorted({ds for _, ds in design})
    target_genes = sorted({g for gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    site_to_stids = kb.neo4j_site_stids(sites)
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, gene_to_stids), site_to_stids)
    totals = {r["readout"] for r in readouts if r["kind"] == "total"}
    units = da.active_units() if da.needs_units() else {}
    actions = da.drug_actions("korkut", sorted(drugs))

    # VERIFY the flag positively on a real network before measuring anything.
    build_root = bm.CATALOG_ROOT.resolve()
    pathways = sorted(p for p in build_root.iterdir() if (p / "stid_to_uuid_mapping.csv").exists())
    probe = bm.parse_pathway_via_ds_api(pathways[0].name)
    pu = probe["nodes"][0]["uuid"]; pv = probe["nodes"][1]["uuid"]
    pr = post("solve", {"network_id": probe["network_id"],
                        "observations": {pu: [50.0, 1.0]},
                        "soft_observations": {pv: [2.0, 0.5]}})
    if pr.get("soft_obs_rule") != "blend" or pr.get("soft_observed") != 1:
        raise SystemExit("the server did not apply a soft observation; is DS_SOFT_OBS=blend set?")
    print("soft-observation probe: blend applied", flush=True)

    out_rows = []
    solves = 0
    for pdir in pathways:
        if shard is not None and pdir.name not in shard:
            continue
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = kb.readout_nodes(r_stids, stid_to_uuids, bm.context_own_uuids(pdir))
        if not ro_nodes:
            continue
        ro_nodes = kb.total_readout_nodes(ro_nodes, totals, pdir)
        parsed = None
        for ds in drug_sets:
            genes = sorted({g for c in ds for g in drugs[c]})
            g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
            if not any(g2u.values()):
                continue
            ds_name = "+".join(ds)
            scored_here = [name for name in sorted(ro_nodes)
                           if (pdir.name, ds_name, name) in recorded]
            if len(scored_here) < 2:
                continue  # nothing to condition on, or nothing to hold out
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed")
            action = da.gene_actions([(actions[c], drugs[c]) for c in ds])
            pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, gene_to_stids,
                                               stid_to_uuids, units, action)
            if not roots:
                continue
            pinned = set(obs)
            for held in scored_here:
                held_nodes = set(ro_nodes[held])
                # conditioning values per node: 2^mean(log2) over the OTHER
                # measured readouts naming it; never a pinned or held-out node
                node_logs = collections.defaultdict(list)
                for other in scored_here:
                    if other == held:
                        continue
                    f = meas_fold.get((ds, other))
                    if f is None or f <= 0:
                        continue
                    for u in ro_nodes[other]:
                        if u in pinned or u in held_nodes:
                            continue
                        node_logs[u].append(math.log2(f))
                if not node_logs:
                    continue
                soft_vals = {u: min(100.0, max(0.0, 2.0 ** (sum(ls) / len(ls))))
                             for u, ls in node_logs.items()}
                free = [u for u in ro_nodes[held] if u not in pinned]
                if not free:
                    continue
                for w in WEIGHTS:
                    payload = {"network_id": parsed["network_id"], "observations": obs,
                               "soft_observations": {u: [v, w] for u, v in soft_vals.items()}}
                    if cut:   # specs/122, exactly as solve_via_ds_api sends them
                        payload["edge_redirects"] = [list(e) for e in sorted(cut)]
                    res = post("solve", payload)
                    if res.get("status") != "success":
                        raise SystemExit(f"{pdir.name} {ds_name} {held} w={w}: solve failed")
                    solves += 1
                    act = res.get("node_activities", {})
                    ui = 100.0 * sum(act.get(u, 0.01) for u in free) / len(free)
                    rec = recorded[(pdir.name, ds_name, held)]
                    out_rows.append({"pathway_id": pdir.name, "drug_set": ds_name,
                                     "readout": held, "weight": w,
                                     "conditioned_nodes": len(soft_vals),
                                     "soft_observed": res.get("soft_observed"),
                                     "converged": int(bool(res.get("converged"))),
                                     "pred_ui": f"{ui:.6g}",
                                     "predicted": kb.NAMES.get(bm.classify(ui), ""),
                                     "baseline_predicted": rec["predicted"],
                                     "measured": rec["measured"]})
        print(f"  {pdir.name}: done ({solves} solves)", flush=True)

    suffix = "" if shard is None else "_" + sorted(shard)[0][:20]
    out = out_dir / f"arm_cases{suffix}.tsv"
    with open(out, "w", newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        wtr.writeheader()
        wtr.writerows(out_rows)
    print(f"{len(out_rows)} held-out rows ({solves} solves) -> {out}")


if __name__ == "__main__":
    main()
