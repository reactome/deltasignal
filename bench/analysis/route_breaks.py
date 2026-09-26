#!/usr/bin/env python3
"""For each failed case where Reactome has a route and our network does not,
find the first step of Reactome's route our network does not reach.

Input is case_oracle.py's --out table. For every missed case with our reach
`no_path` and Reactome route `strict`, Reactome's shortest route (case_oracle's
graph) is walked against the nodes our network reaches from the pinned uuids.
A set counts as reached if any member is. The first unreached step is typed:

  component -> complex   an assembly step (node exists but edge missing, or
                         the complex is absent from our network)
  member -> set          a member standing in for a set
  entity -> reaction     an input, catalyst or regulator
  reaction -> output
and two outcomes that are not a break: the route starts from a form of the
gene we did not pin, or the whole route is reached (readout mapping).

Usage: route_breaks.py --oracle ARM/curator_case_oracle.tsv --catalog BUILD
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from case_oracle import MADE, USE, faithful_comp_pairs, fetch_one_level, oracle_graph  # noqa: E402
from curator_oracle import fetch_reaction_rows, shortest_path  # noqa: E402

RX = ("Reaction", "BlackBoxEvent", "Polymerisation", "Depolymerisation", "FailedReaction")


def first_break(base_path, reached_ids, expand, present_ids, klass):
    """(label, a, b) for the first step whose target is not reached."""
    i = next((k for k, s in enumerate(base_path) if not (expand(s) & reached_ids)), None)
    if i is None:
        return "whole route reached", "", ""
    if i == 0:
        return "route starts from another form of the gene", "", ""
    a, b = base_path[i - 1], base_path[i]
    ca, cb = klass(a), klass(b)
    step = ("entity -> reaction" if cb in RX else "reaction -> output" if ca in RX
            else "member -> set" if "Set" in cb else "component -> complex")
    where = "node exists, edge missing" if expand(b) & present_ids else "node absent"
    return f"{step} | {where}", a, b


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--oracle", type=Path, required=True)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--steps", type=Path, help="write one row per case: pathway, gene, break, from, to")
    ap.add_argument("--faithful", action="store_true",
                    help="walk the route using only representable hops (faithful_comp_pairs) where one exists")
    a = ap.parse_args()
    import benchmark_vs_mpbiopath as BM
    from py2neo import Graph
    g = Graph(BM.NEO4J_URL, auth=(BM.NEO4J_USER, BM.NEO4J_PASSWORD))
    cls: dict = {}

    def klass(s):
        if s not in cls:
            r = g.run("MATCH (n {stId:$s}) RETURN n.schemaClass AS c", s=s).data()
            cls[s] = r[0]["c"] if r else "?"
        return cls[s]

    names = {}
    for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv"):
        if line.startswith("#") or line.startswith("id\t"):
            continue
        sid, name = line.rstrip("\n").split("\t")
        names[name] = sid
    rows = [r for r in csv.DictReader(open(a.oracle, newline=""), delimiter="\t")
            if r["error"] == "missed" and r["reach"] == "no_path" and r["reactome_route"] == "strict"]
    graphs: dict = {}
    out = collections.Counter()
    per = collections.Counter()
    steps = []
    for r in rows:
        pw, pid = r["pathway"], names[r["pathway"]]
        if pw not in graphs:
            comp, sets = fetch_one_level(pid)
            rxn_rows = fetch_reaction_rows(pid)
            adj = oracle_graph(rxn_rows, comp, sets)
            # Reactome's route using only hops the generator can represent.
            fadj = oracle_graph(rxn_rows, faithful_comp_pairs(rxn_rows, comp), sets)
            members = collections.defaultdict(set)
            for x, m in sets:
                members[x].add(m)

            def expand(s, depth=0, members=members):
                if s in members and depth < 5:
                    return {s}.union(*(expand(m, depth + 1) for m in members[s]))
                return {s}

            nodes = list(csv.DictReader(open(a.catalog / pid / "nodes.csv", newline="")))
            ident = {x["uuid"]: x["diagram_entity_id"].split("::")[0] for x in nodes}
            fwd = collections.defaultdict(list)
            for e in csv.DictReader(open(a.catalog / pid / "logic_network.csv", newline="")):
                fwd[e["source_id"]].append(e["target_id"])
            graphs[pw] = ({k: {t for t, _ in v} for k, v in adj.items()}, ident, fwd,
                          set(ident.values()), expand,
                          {k: {t for t, _ in v} for k, v in fadj.items()})
        uadj, ident, fwd, present, expand, fuadj = graphs[pw]
        gene = BM.GENE_NAME_CORRECTIONS.get((pw, r["gene"]), r["gene"])
        starts = {v for x in BM.neo4j_gene_to_stids([gene]).get(gene, []) for v in (x, x + USE, x + MADE)}
        goal = BM.neo4j_dbid_to_stid({r["key_output"]}).get(r["key_output"])
        goals = {goal, goal + USE, goal + MADE} if goal else set()
        path = shortest_path(uadj, starts, goals) if goal else None
        if not path:
            out["(no route found on re-walk)"] += 1
            continue
        fpath = shortest_path(fuadj, starts, goals)
        if a.faithful and fpath:
            path = fpath      # break on the route the generator could have built
        seen = set(r["gene_uuids"].split("|"))
        stack = list(seen)
        while stack:
            n = stack.pop()
            for t in fwd[n]:
                if t not in seen:
                    seen.add(t)
                    stack.append(t)
        label, frm, to = first_break([p.split("::")[0] for p in path], {ident[u] for u in seen},
                                  expand, present, klass)
        if a.faithful:
            label = f"{label} | {'faithful route' if fpath else 'composition-only'}"
        out[label] += 1
        per[(pw, label)] += 1
        steps.append((pw, r["gene"], r["direction"], r["key_output"], label, frm, to))
    print(f"{len(rows)} cases where Reactome has a route and ours does not; first unreached step:")
    for k, v in out.most_common():
        print(f"  {v:4d}  {k}")
    print("top (pathway, break):", per.most_common(6))
    if a.steps:
        with open(a.steps, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pathway", "gene", "direction", "key_output", "break", "from", "to"])
            w.writerows(steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
