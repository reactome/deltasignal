#!/usr/bin/env python3
"""For each failed case, ask Reactome itself whether the perturbed gene reaches
the readout, and so whether the logic network, the solver, or neither is at
fault.

Reactome's graph is built from the pathway in Neo4j: reactions (input, output,
catalyst, signed regulator), complex assembly one level at a time (component
-> complex, chaining through nested complexes), and sets with SPLIT use and
production: a member may stand in for a set where the set is used
(member -> set_use), and a set a reaction produces hands on to its members
(set_made -> member), but no path runs member -> set -> sibling member.
curator_oracle.py's reading treats a set as its members in both directions,
which lets a signal hop between siblings (PMAIP1 -> "Bcl-XL interacting BH3-
only proteins" -> tBID, traced) and manufactured routes. Small molecules
cannot carry a signal. LENIENT adds complex -> component, the direction our
`dissociation` edges take.

Each case gets `reactome_route`:
  strict      Reactome connects gene -> readout without dissociation
              (sign: matches / opposite_only / both_parities -- the last means
              a walk through a cycle with an inhibition gives either sign, so the
              sign is not evidence)
  lenient     only through a complex releasing a component
  none        Reactome does not connect them within this pathway at all

and, for strict routes, `reactome_sign`: whether some route's sign parity
matches the expected direction (+1 = the readout should follow the gene).

Combined with our own reach (failure_structure.py's `reach`):
  none                   -> the curator reasoned beyond this pathway; no
                            within-pathway model can get it
  strict, ours no_path   -> the logic network severs a route Reactome has
  strict, ours routed    -> the propagator gets it wrong on a real route

Usage:
  case_oracle.py --annotated ARM/curator_failure_structure.tsv --catalog BUILD \
      [--split held-out|all] [--out annotated.tsv]
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from curator_oracle import fetch_reaction_rows  # noqa: E402
from holdout_report import TUNING_PATHWAYS  # noqa: E402


USE, MADE = "::use", "::made"


def oracle_graph(rows, comp_pairs, set_pairs, lenient=False) -> dict:
    """Signed graph over Reactome ids with split set nodes (see module doc).

    rows: fetch_reaction_rows tuples (reaction, class, input, output, catalyst,
    regulator, regulation class). comp_pairs: (complex, component), one level.
    set_pairs: (set, member), one level. Only R-HSA ids carry a signal."""
    is_set = {x for x, _ in set_pairs}
    use = lambda e: e + USE if e in is_set else e
    made = lambda e: e + MADE if e in is_set else e
    adj: dict = collections.defaultdict(set)
    for r, _k, i, o, c, g, gk in rows:
        if i: adj[use(i)].add((r, 1))
        if o: adj[r].add((made(o), 1))
        if c: adj[use(c)].add((r, 1))
        if g: adj[use(g)].add((r, -1 if "Negative" in gk else 1))
    for x in is_set:
        adj[x + MADE].add((x + USE, 1))          # a produced set can be used as the set
    for x, m in set_pairs:
        adj[use(m) if m in is_set else m].add((x + USE, 1))    # a member stands in for the set
        adj[x + MADE].add((made(m) if m in is_set else m, 1))  # a produced set hands on to members
    for x, m in comp_pairs:
        # A set that is a COMPONENT is being used by the complex: its members
        # reach it through m::use (review of PR #74: wiring it from m::made cut
        # PDGFB -> "Active PDGF dimers" -> receptor complex).
        adj[use(m)].add((x, 1))                                  # assembly
        if lenient:
            adj[x].add((made(m), 1))                             # release
    base = lambda n: n.split("::")[0]
    # Only R-HSA ids CARRY a signal; a small molecule (R-ALL) may still be a
    # READOUT, so edges into it are kept and edges out of it dropped. Dropping
    # both made every R-ALL readout unreachable by construction.
    return {a: {(b, sg) for b, sg in ts if base(b).startswith(("R-HSA-", "R-ALL-"))}
            for a, ts in adj.items() if base(a).startswith("R-HSA-")}


def fetch_one_level(pid: str):
    """(complex, component) and (set, member) pairs, one level at a time, for
    every participant of the pathway's reactions and everything nested in them."""
    from curator_oracle import _cypher
    q = (f"MATCH (p:Pathway {{stId:'{pid}'}})-[:hasEvent*]->(r:ReactionLikeEvent) WITH DISTINCT r "
         # catalystActivity->physicalEntity and regulatedBy->regulator: without
         # the second hop's labels a catalyst- or regulator-only entity was
         # never expanded (review of PR #74).
         "MATCH (r)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator*1..2]->(x) "
         "WHERE x:PhysicalEntity WITH DISTINCT x "
         "MATCH (x)-[:hasComponent|hasMember|hasCandidate*0..4]->(y)-[rel:hasComponent|hasMember|hasCandidate]->(z) "
         "WHERE y.stId IS NOT NULL AND z.stId IS NOT NULL "
         "RETURN DISTINCT y.stId + '|' + z.stId + '|' + type(rel)")
    comp, sets = set(), set()
    for row in _cypher(q):
        if row.count("|") == 2:
            y, z, t = row.split("|")
            (comp if t == "hasComponent" else sets).add((y, z))
    return comp, sets


def signed_parities(adj: dict, starts: set, goals: set, limit: int = 2_000_000) -> set:
    """Sign parities (+1 / -1) with which any goal is reachable from any start."""
    seen = {(s, 1) for s in starts}
    stack = list(seen)
    found = set()
    while stack:
        if len(seen) >= limit:
            raise RuntimeError(f"signed_parities: {limit} states exceeded; result would be partial")
        n, sg = stack.pop()
        if n in goals:
            found.add(sg)
            if len(found) == 2:
                break
        for m, s in adj.get(n, ()):
            st = (m, sg * s)
            if st not in seen:
                seen.add(st)
                stack.append(st)
    return found


def classify(strict_adj, lenient_adj, starts, goals, expected, direction) -> tuple[str, str]:
    """(reactome_route, reactome_sign)."""
    par = signed_parities(strict_adj, starts, goals)
    if par:
        if expected == "1":
            return "strict", "n/a"
        need = 1 if expected == direction else -1
        # signed_parities finds WALKS: a route touching a cycle that contains an
        # inhibition yields both parities, which is weak evidence of the sign.
        if par == {1, -1}:
            return "strict", "both_parities"
        return "strict", "matches" if need in par else "opposite_only"
    if signed_parities(lenient_adj, starts, goals):
        return "lenient", "n/a"
    return "none", "n/a"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", type=Path, required=True,
                    help="failure_structure.py --out table (has gene/output uuids and reach)")
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--split", choices=("held-out", "all"), default="held-out")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    import benchmark_vs_mpbiopath as BM
    names = {}
    for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv"):
        if line.startswith("#") or line.startswith("id\t"):
            continue
        sid, name = line.rstrip("\n").split("\t")
        names[name] = sid

    rows = [r for r in csv.DictReader(open(a.annotated, newline=""), delimiter="\t")
            if r["correct"] != "True"
            and (a.split == "all" or r["pathway"] not in TUNING_PATHWAYS)]
    graphs: dict = {}
    out = []
    for r in rows:
        pw = r["pathway"]; pid = names[pw]
        if pw not in graphs:
            rows_ = fetch_reaction_rows(pid)
            comp, sets = fetch_one_level(pid)
            strict = oracle_graph(rows_, comp, sets, lenient=False)
            lenient = oracle_graph(rows_, comp, sets, lenient=True)
            graphs[pw] = (strict, lenient, BM.neo4j_gene_to_stids, {})
        strict, lenient, g2s, cache = graphs[pw]
        gene = BM.GENE_NAME_CORRECTIONS.get((pw, r["gene"]), r["gene"])
        if gene not in cache:
            cache[gene] = set(g2s([gene]).get(gene, []))
        ko = r["key_output"]
        kkey = ("ko", ko)
        if kkey not in cache:
            cache[kkey] = BM.neo4j_dbid_to_stid({ko}).get(ko)
        starts, goal = cache[gene], cache[kkey]
        starts = {v for x in starts for v in (x, x + USE, x + MADE)}
        goals = {goal, goal + USE, goal + MADE} if goal else set()
        route, sign = ("none", "n/a") if not (starts and goals) else \
            classify(strict, lenient, starts, goals, r["expected"], r["direction"])
        out.append(dict(r, reactome_route=route, reactome_sign=sign))

    c = collections.Counter((r["error"], r["reach"], r["reactome_route"], r["reactome_sign"]) for r in out)
    print(f"{a.split}: {len(out)} failed cases")
    print(f"  {'error':<16}{'our reach':<17}{'Reactome route':<16}{'Reactome sign':<15}cases")
    for (e, rc, rt, sg), n in sorted(c.items(), key=lambda kv: -kv[1]):
        print(f"  {e:<16}{rc:<17}{rt:<16}{sg:<15}{n}")
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(out)
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
