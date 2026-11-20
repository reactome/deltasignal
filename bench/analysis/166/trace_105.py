#!/usr/bin/env python3
"""specs/166 side question: korkut pairs predicted at exactly 1 although a directed path runs from a
node containing a drug target to a node containing the readout. Classify where the perturbation stops.

Env: DS_CATALOG_ROOT=<without build>, DELTASIGNAL_BASE=<server serving it>.
Usage: trace_105.py REPO KORKUT_CASES UNREACHED_KORKUT OUT_TSV
"""
import csv
import os
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "bench"))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import korkut_bench as kb  # noqa: E402
import drug_activity as da  # noqa: E402

BASE = 0.01
EPS = 1e-12


def rev_reach(radj, targets):
    seen, q = set(targets), deque(targets)
    while q:
        v = q.popleft()
        for u in radj.get(v, ()):
            if u not in seen:
                seen.add(u)
                q.append(u)
    return seen


def main():
    cases_path, unreached_path, out_path = sys.argv[2:5]
    build = bm.CATALOG_ROOT.resolve()
    unreached = {(r["pathway_stid"], r["drug_set"], r["readout"])
                 for r in csv.DictReader(open(unreached_path), delimiter="\t")}
    cases = {(r["pathway_id"], r["drug_set"], r["readout"]): r
             for r in csv.DictReader(open(cases_path), delimiter="\t")}
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    target_genes = sorted({g for gs in drugs.values() for g in gs})
    readout_genes = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(sorted(set(target_genes) | set(readout_genes)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, gene_to_stids), kb.neo4j_site_stids(sites))
    totals = {r["readout"] for r in readouts if r["kind"] == "total"}
    out = []
    for pdir in sorted(p for p in build.iterdir() if (p / "logic_network.csv").exists()):
        mine = sorted(k for k in unreached if k[0] == pdir.name)
        if not mine:
            continue
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = kb.total_readout_nodes(kb.readout_nodes(r_stids, stid_to_uuids, None), totals, pdir)
        rows = list(csv.DictReader(open(pdir / "logic_network.csv")))
        adj, radj = defaultdict(list), defaultdict(list)
        etype = defaultdict(set)
        for r in rows:
            adj[r["source_id"]].append(r["target_id"])
            radj[r["target_id"]].append(r["source_id"])
            etype[(r["source_id"], r["target_id"])].add(r["edge_type"])
        kinds = {r["uuid"]: r["node_kind"] for r in csv.DictReader(open(pdir / "nodes.csv"))}
        held_tables = {}
        for name in ("cofactors", "drugs", "pathogens"):
            f = pdir / f"{name}.csv"
            if f.exists():
                rr = list(csv.DictReader(open(f)))
                col = "uuid" if rr and "uuid" in rr[0] else (list(rr[0])[0] if rr else None)
                held_tables[name] = {r[col] for r in rr} if col else set()
        parsed = None
        for k in mine:
            _, ds, ro = k
            c = cases.get(k)
            genes = sorted({g for code in ds.split("+") for g in drugs[code]})
            t_all = {u for g in genes for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])}
            R = set(ro_nodes.get(ro, []))
            rec = {"pathway": pdir.name, "drug_set": ds, "readout": ro,
                   "status": c["status"] if c else "absent", "pred_ui": c["pred_ui"] if c else ""}
            if not t_all or not R:
                rec["class"] = "no_target_or_readout_node"
                out.append(rec)
                continue
            ro_genes = sorted({x.split(":")[0] for r in readouts if r["readout"] == ro for x in r["sites"].split(",")})
            ro_kind = next(r["kind"] for r in readouts if r["readout"] == ro)
            R_any = {u for g in ro_genes for s_ in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s_, [])}
            reach_t = bm.reachable_from(adj, list(t_all))
            rec["readout_kind"] = ro_kind
            rec["any_form_path"] = int(bool(reach_t & R_any))
            if not (reach_t & R):
                rec["class"] = ("no_structural_path" if not rec["any_form_path"]
                                else f"path_reaches_other_forms_not_the_scored_{ro_kind}_form")
                out.append(rec)
                continue
            rec["structural_path"] = 1
            if not c or c["status"] != "scored" or float(c["pred_ui"]) != 1.0:
                rec["class"] = "not_exactly_1"
                out.append(rec)
                continue
            no_solve = os.environ.get("TRACE_NO_SOLVE") == "1"
            if parsed is None and not no_solve:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
            action = da.gene_actions([(da.drug_actions("korkut", [code])[code], drugs[code])
                                      for code in ds.split("+")])
            pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, gene_to_stids, stid_to_uuids, {}, action)
            free_R = R - set(obs)
            if not roots:
                rec["class"] = "no_root_pinned"
            elif not (bm.reachable_from(adj, list(roots)) & free_R):
                rec["class"] = "pinned_roots_do_not_reach_readout"
                rec["non_root_target_forms_reaching"] = len(
                    {u for u in t_all - roots if bm.reachable_from(adj, [u]) & free_R})
            elif no_solve:
                rec["class"] = "roots_reach_readout_needs_solve"
            else:
                res = bm.solve_via_ds_api(None, obs, network_id=parsed["network_id"], edge_redirects=cut)
                act = res.get("node_activities", {})
                moved = {u for u, a in act.items() if abs(a - BASE) > EPS}
                to_R = rev_reach(radj, list(free_R))
                # frontier: a moved node with an edge into an unmoved node that still leads to the readout
                front = Counter()
                examples = []
                for u in moved & to_R:
                    for v in adj.get(u, ()):
                        if v in to_R and v not in moved:
                            for et in etype[(u, v)]:
                                held = [n for n, s in held_tables.items() if v in s]
                                lab = f"{et}->{kinds.get(v, '?')}" + (f"[held:{','.join(held)}]" if held else "")
                                front[lab] += 1
                                if len(examples) < 3:
                                    examples.append(f"{u}->{v}:{lab}")
                rec["class"] = "stops_on_route" if front else "moved_set_never_reaches"
                rec["frontier"] = ";".join(f"{a}={n}" for a, n in front.most_common(5))
                rec["examples"] = " ".join(examples)
                rec["moved_nodes"] = len(moved)
            out.append(rec)
        print(f"  {pdir.name}: {len(mine)} pairs", flush=True)
    fields = ["pathway", "drug_set", "readout", "readout_kind", "any_form_path", "status", "pred_ui", "structural_path", "class",
              "non_root_target_forms_reaching", "moved_nodes", "frontier", "examples"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(Counter(r["class"] for r in out))


if __name__ == "__main__":
    main()
