#!/usr/bin/env python3
"""Inventory of how catalysts and regulators are wired into a catalog build
(specs/034). Every catalyst/regulator EDGE into a reaction node is classified on
four axes:

  role    catalyst | positive | negative | requirement | pos_gene_expr | neg_gene_expr
          (from the Regulation subclass in Neo4j; the generator only knows
          catalyst / PositiveRegulation / NegativeRegulation)
  shape   of the curated entity: protein | complex | complex_with_set | set |
          set_of_complexes | small_molecule | drug-derived (drugs.csv when present)
  copy    of the SOURCE node that is wired in:
            root          no incoming edge (what the benchmark pins)
            upstream      produced, and the reaction does not reach it
            downstream    produced, and the reaction reaches it (a loop closes)
            pool          an LNG_SET_POOL node
  combine the edge's and_or flag

With --cases, each curator/experimental case is attributed to the regulator
edges its perturbed gene's pinned nodes can reach through a regulator edge, and
accuracy is tallied per (role, shape, copy) cell.

Usage: regulator_inventory.py --catalog BUILD_DIR [--cases ARM_DIR]
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from network_files import open_network  # noqa: E402  (specs/044)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROLE_Q = """
UNWIND $rids AS rid MATCH (r:ReactionLikeEvent {stId: rid})
CALL {
  WITH r MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(e)
  RETURN 'catalyst' AS role, e.stId AS e
  UNION WITH r MATCH (r)-[:regulatedBy]->(rg)-[:regulator]->(e)
  RETURN CASE WHEN rg:Requirement THEN 'requirement'
              WHEN rg:PositiveGeneExpressionRegulation THEN 'pos_gene_expr'
              WHEN rg:NegativeGeneExpressionRegulation THEN 'neg_gene_expr'
              WHEN rg:PositiveRegulation THEN 'positive'
              ELSE 'negative' END AS role, e.stId AS e
}
RETURN rid AS r, role, e
"""

SHAPE_Q = """
UNWIND $ids AS s MATCH (e:PhysicalEntity {stId: s})
RETURN s,
  CASE WHEN e:SimpleEntity THEN 'small_molecule'
       WHEN e:EntitySet AND size([(e)-[:hasMember|hasCandidate]->(m:Complex) | m]) > 0 THEN 'set_of_complexes'
       WHEN e:EntitySet THEN 'set'
       WHEN e:Complex AND size([(e)-[:hasComponent*1..6]->(x:EntitySet) | x]) > 0 THEN 'complex_with_set'
       WHEN e:Complex THEN 'complex'
       ELSE 'protein' END AS shape,
  [(e)-[:hasComponent|hasMember|hasCandidate*0..8]->(l:PhysicalEntity) | l.stId] AS leaves
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--cases", type=Path, help="arm dir with {curator,experimental}_cases.tsv")
    ap.add_argument("--cell", help="role/shape/copy: list its edges and the FAILED cases routed through them")
    a = ap.parse_args()
    import benchmark_vs_mpbiopath as BM
    from py2neo import Graph
    g = Graph(BM.NEO4J_URL, auth=(BM.NEO4J_USER, BM.NEO4J_PASSWORD))

    names = {}
    for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv"):
        if line.startswith(("#", "id\t")):
            continue
        sid, name = line.rstrip("\n").split("\t")
        names[sid] = name

    cells = collections.Counter()          # (role, shape, copy, and_or) -> edges
    want = tuple(a.cell.split("/")) if a.cell else None
    cell_edges = collections.Counter()     # (pathway, regulator stid, reaction stid) -> failed cases
    case_cells = collections.defaultdict(lambda: [0, 0])   # (role, shape, copy) -> [n, correct]
    cases = {}
    if a.cases:
        for ax in ("curator", "experimental"):
            p = a.cases / f"{ax}_cases.tsv"
            if p.exists():
                for r in csv.DictReader(open(p), delimiter="\t"):
                    if r["valid"] == "1" and not r.get("exclusion"):
                        cases.setdefault(r["pathway"], []).append((ax, r))

    for d in sorted(a.catalog.glob("R-HSA-*")):
        pid = d.name
        lab = {x["uuid"]: x["stable_id"] for x in csv.DictReader(open(d / "stid_to_uuid_mapping.csv"))}
        kinds = {x["uuid"]: x["node_kind"] for x in csv.DictReader(open(d / "nodes.csv"))}
        drugs = set()
        if (d / "drugs.csv").exists():
            drugs = {x["stable_id"] for x in csv.DictReader(open(d / "drugs.csv"))}
        E = list(csv.DictReader(open_network(d)))
        fwd = collections.defaultdict(set)
        indeg = collections.Counter()
        for e in E:
            fwd[e["source_id"]].add(e["target_id"])
            indeg[e["target_id"]] += 1
        reach_cache: dict = {}

        def reach(u):
            if u not in reach_cache:
                seen, st = {u}, [u]
                while st:
                    for t in fwd[st.pop()]:
                        if t not in seen:
                            seen.add(t)
                            st.append(t)
                reach_cache[u] = seen
            return reach_cache[u]

        reg_edges = [e for e in E if e["edge_type"] in ("catalyst", "regulator")]
        rx_stids = sorted({lab.get(e["target_id"], "") for e in reg_edges} - {""})
        roles = collections.defaultdict(list)       # (reaction stid) -> [(role, entity stid)]
        for r in g.run(ROLE_Q, rids=rx_stids).data():
            roles[r["r"]].append((r["role"], r["e"]))
        ent = sorted({e for v in roles.values() for _, e in v})
        shape, leaves = {}, {}
        for r in g.run(SHAPE_Q, ids=ent).data():
            shape[r["s"]], leaves[r["s"]] = r["shape"], set(r["leaves"])
        edge_cell = {}
        for e in reg_edges:
            src, tgt = e["source_id"], e["target_id"]
            rstid, sstid = lab.get(tgt, ""), lab.get(src, "")
            want_neg = e["pos_neg"] == "neg"
            # the curated annotation this edge comes from: the entity equal to,
            # or containing, the source node's entity, with a matching sign
            cands = [(ro, en) for ro, en in roles.get(rstid, [])
                     if (ro in ("negative", "neg_gene_expr")) == want_neg
                     and (en == sstid or sstid in leaves.get(en, set()))]
            if e["edge_type"] == "catalyst":
                cands = [c for c in cands if c[0] == "catalyst"] or cands
            else:
                cands = [c for c in cands if c[0] != "catalyst"] or cands
            role, en = cands[0] if cands else ("unmatched", sstid)
            shp = "drug" if (en in drugs or sstid in drugs) else shape.get(en, "?")
            if kinds.get(src) == "set_pool":
                copy = "pool"
            elif indeg[src] == 0:
                copy = "root"
            elif src in reach(tgt):
                copy = "downstream"
            else:
                copy = "upstream"
            cells[(role, shp, copy, e["and_or"])] += 1
            edge_cell[(src, tgt)] = (role, shp, copy)
        # attribute cases: a case uses a regulator edge if its pins reach the source
        for ax, r in cases.get(names.get(pid, ""), []):
            pins = [u for u in r["gene_uuids"].split("|") if u]
            if not pins:
                continue
            got = set()
            seen, st = set(pins), list(pins)
            while st:
                for t in fwd[st.pop()]:
                    if t not in seen:
                        seen.add(t)
                        st.append(t)
            outs = set(u for u in r["output_uuids"].split("|") if u)
            ok = r["predicted"] == r["expected"]
            for (s, t), c in edge_cell.items():
                if s in seen and outs & reach(t):
                    got.add(c)
                    if want and c == want and not ok:
                        cell_edges[(ax, names.get(pid, pid)[:30], r["gene"], r["direction"],
                                    lab.get(s, s), lab.get(t, t))] += 1
            for c in got:
                key = (ax,) + c
                case_cells[key][0] += 1
                case_cells[key][1] += ok

    print("role / shape / copy / and_or -> edges")
    for k, v in sorted(cells.items(), key=lambda kv: -kv[1]):
        print(f"  {v:7d}  {' / '.join(k)}")
    if want:
        print(f"\nFAILED cases through {a.cell} (axis, pathway, gene, dir, regulator node, reaction):")
        for k, v in cell_edges.most_common(25):
            print(f"  {v:4d}  {k}")
        return 0
    if case_cells:
        print("\ncases whose route uses a regulator edge of the cell (a case may use many cells):")
        for k, (n, ok) in sorted(case_cells.items(), key=lambda kv: -kv[1][0]):
            if n >= 20:
                print(f"  {k[0]:12} {' / '.join(k[1:]):45} n={n:6d} acc={ok / n:.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
