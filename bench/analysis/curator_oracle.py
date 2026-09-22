#!/usr/bin/env python3
"""Curator oracle: read a pathway straight from Reactome and ask what the solver can reach.

Method (2026-09-19): the curators reasoned about the pathway; do the same
reading from Neo4j and check whether the network the solver runs on arrives
where that reading arrives. Where it does not, that is a construction defect on
our side regardless of what the 2019 curator wrote. It found three (specs/016):
the composition-hierarchy gap, set-member complexes shattered to leaves, and --
refuted -- dropped set-typed regulators.

Oracle graph (STRICT, i.e. what a curator counts):
  * every ReactionLikeEvent under the pathway: input, output, catalyst, signed regulator
  * composition in the ASSEMBLY direction only: component -> containing complex
    (the dissociation direction is the broadcast route measured at -15pp)
  * set membership both ways: a set IS its members
  * roots = protein-bearing entities with no incoming REACTION edge, not cofactors;
    terminals = no outgoing reaction edge

Then, against LNG's bundle for the same pathway:
  * stable-id level (sets resolved to leaves): does a route exist?
  * uuid level (the graph the solver actually runs on): does a route exist?
  * optionally, SIMULATE composition edges at uuid level (component node ->
    containing-complex node, <= N hasComponent hops) and recount

Usage:
  curator_oracle.py R-HSA-909733 --catalog /path/to/catalog [--simulate-hops 2]

The graph algorithms are pure functions over plain dicts so they are unit-tested
without Neo4j (bench/test_analysis_numbers.py). Only `fetch_*` talk to Neo4j.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Iterable

SIGNED_ADJ = dict[str, set[tuple[str, int]]]
ADJ = dict[str, set[str]]


# ----------------------------------------------------------------------------
# pure graph functions (tested)
# ----------------------------------------------------------------------------

def build_reaction_graph(rows: Iterable[tuple[str, str, str, str, str, str, str]]):
    """From (reaction, class, input, output, catalyst, regulator, regulation_class) rows.

    Returns (signed_adj, reaction_only_adj, entities, reactions, role) where role
    maps an entity to the set of roles it plays (input/output/catalyst/pos_regulator/
    neg_regulator). Empty strings mean "no such participant on this row".
    """
    adj: SIGNED_ADJ = collections.defaultdict(set)
    role: dict[str, set[str]] = collections.defaultdict(set)
    ents: set[str] = set()
    rxns: set[str] = set()
    for r, _rk, i, o, c, g, gk in rows:
        rxns.add(r)
        if i:
            adj[i].add((r, +1)); ents.add(i); role[i].add("input")
        if o:
            adj[r].add((o, +1)); ents.add(o); role[o].add("output")
        if c:
            adj[c].add((r, +1)); ents.add(c); role[c].add("catalyst")
        if g:
            neg = "Negative" in gk
            adj[g].add((r, -1 if neg else +1)); ents.add(g)
            role[g].add("neg_regulator" if neg else "pos_regulator")
    reaction_only = {k: set(v) for k, v in adj.items()}
    return adj, reaction_only, ents, rxns, role


def add_composition(adj: SIGNED_ADJ, complex_pairs: Iterable[tuple[str, str]],
                    set_pairs: Iterable[tuple[str, str]], strict: bool = True) -> None:
    """Add curator-readable composition to a signed graph, in place.

    complex_pairs are (container, member) for hasComponent; set_pairs are
    (set, member) for hasMember/hasCandidate. STRICT adds only member -> container
    for complexes (assembly direction); sets get both directions.
    """
    for x, m in complex_pairs:
        adj.setdefault(m, set()).add((x, +1))
        if not strict:
            adj.setdefault(x, set()).add((m, +1))
    for x, m in set_pairs:
        adj.setdefault(m, set()).add((x, +1))
        adj.setdefault(x, set()).add((m, +1))


def roots_and_terminals(reaction_only: SIGNED_ADJ, ents: set[str],
                        cofactors: set[str] = frozenset(), prefix: str = "R-HSA-"):
    indeg: collections.Counter = collections.Counter()
    outdeg: collections.Counter = collections.Counter()
    for a, ts in reaction_only.items():
        for t, _ in ts:
            outdeg[a] += 1; indeg[t] += 1
    roots = sorted(e for e in ents if indeg[e] == 0 and e.startswith(prefix) and e not in cofactors)
    terms = sorted(e for e in ents if outdeg[e] == 0 and e.startswith(prefix))
    return roots, terms


def drop_carriers(adj: SIGNED_ADJ, is_carrier: "callable[[str], bool]") -> SIGNED_ADJ:
    """Remove nodes that must not CARRY signal (small molecules, cofactors) from a
    signed graph: no edges into or out of them. A curator does not route a
    perturbation through GTP or ATP, and the solver pins cofactors inert, so a
    route whose only link is a released small molecule is not a curator route.
    Without this, Rho GTPases showed 3,384 'severed' routes that all died on a
    GTP dissociation sink -- an artifact of the oracle, not of the network."""
    out: SIGNED_ADJ = {}
    for a, ts in adj.items():
        if not is_carrier(a):
            continue
        kept = {(t, sg) for t, sg in ts if is_carrier(t)}
        if kept:
            out[a] = kept
    return out


def signed_reach(adj: SIGNED_ADJ, start: str) -> dict[str, int]:
    """Every node reachable from `start`, with the sign of the first path found."""
    seen = {start: +1}
    stack = [start]
    while stack:
        n = stack.pop()
        for m, s in adj.get(n, ()):
            if m not in seen:
                seen[m] = seen[n] * s
                stack.append(m)
    return seen


def reachable(adj: ADJ, starts: set[str], goals: set[str]) -> bool:
    seen = set(starts); stack = list(starts)
    while stack:
        n = stack.pop()
        for m in adj.get(n, ()):
            if m in goals:
                return True
            if m not in seen:
                seen.add(m); stack.append(m)
    return False


def shortest_path(adj: ADJ, starts: set[str], goals: set[str]) -> list[str] | None:
    prev: dict[str, str | None] = {s: None for s in starts}
    queue = collections.deque(starts)
    while queue:
        n = queue.popleft()
        if n in goals:
            path = []
            while n is not None:
                path.append(n); n = prev[n]
            return path[::-1]
        for m in adj.get(n, ()):
            if m not in prev:
                prev[m] = n; queue.append(m)
    return None


def first_break(stid_path: list[str], start_uuids: set[str],
                edge_uuids: dict[tuple[str, str], list[tuple[str, str]]]) -> tuple[int, set[str]] | None:
    """Walk a stable-id path at uuid level. Return (step, frontier) at the first step
    where the currently reachable shards of one entity have no edge into the next,
    or None if the whole path is traversable."""
    frontier = set(start_uuids)
    for i in range(1, len(stid_path)):
        hit = {b for a, b in edge_uuids.get((stid_path[i - 1], stid_path[i]), []) if a in frontier}
        if not hit:
            return i, frontier
        frontier = hit
    return None


def resolve_ids(stid: str, present: set[str], leaves_of: dict[str, set[str]]) -> set[str]:
    """Stable ids the bundle may use for `stid`: itself if present, else its expanded leaves."""
    if stid in present:
        return {stid}
    return {m for m in leaves_of.get(stid, ()) if m in present}


def simulate_composition(uadj: ADJ, direct: Iterable[tuple[str, str]], s2u: dict[str, list[str]],
                         sink_uuids: set[str]) -> tuple[ADJ, int, list[int]]:
    """Copy `uadj` with component-node -> containing-complex-node edges added for each
    (component_stid, container_stid). Sinks are never sources. Returns (adj, n_added, fan_outs)."""
    out = {k: set(v) for k, v in uadj.items()}
    added = 0; fan: list[int] = []
    for x, y in direct:
        xs = [u for u in s2u.get(x, []) if u not in sink_uuids]
        ys = s2u.get(y, [])
        if not xs or not ys:
            continue
        fan.append(len(ys))
        for ux in xs:
            for uy in ys:
                if uy not in out.setdefault(ux, set()):
                    out[ux].add(uy); added += 1
    return out, added, fan


# ----------------------------------------------------------------------------
# Neo4j fetch layer (not unit-tested; thin)
# ----------------------------------------------------------------------------

def _cypher(query: str, container: str = "reactome-neo4j") -> list[str]:
    pw = subprocess.run(["docker", "exec", container, "printenv", "NEO4J_AUTH"],
                        capture_output=True, text=True).stdout.strip().split("/")[-1]
    out = subprocess.run(["docker", "exec", container, "cypher-shell", "-u", "neo4j", "-p", pw,
                          "--format", "plain", query], capture_output=True, text=True)
    if out.returncode:
        raise RuntimeError(out.stderr[:500])
    return [line.strip().strip('"') for line in out.stdout.splitlines()[1:] if line.strip()]


def fetch_reaction_rows(pathway: str):
    q = f"""
    MATCH (p:Pathway {{stId:'{pathway}'}})-[:hasEvent*]->(r:ReactionLikeEvent) WITH DISTINCT r
    OPTIONAL MATCH (r)-[:input]->(i)  OPTIONAL MATCH (r)-[:output]->(o)
    OPTIONAL MATCH (r)-[:catalystActivity]->(ca)-[:physicalEntity]->(c)
    OPTIONAL MATCH (r)-[:regulatedBy]->(reg)-[:regulator]->(g)
    RETURN r.stId + '|' + r.schemaClass + '|' + coalesce(i.stId,'') + '|' + coalesce(o.stId,'') + '|' +
           coalesce(c.stId,'') + '|' + coalesce(g.stId,'') + '|' + coalesce(reg.schemaClass,'')
    """
    rows = []
    for row in _cypher(q):
        if row.count("|") == 6:
            rows.append(tuple(row.split("|")))
    return rows


def fetch_composition(pathway: str):
    q = (f"MATCH (p:Pathway {{stId:'{pathway}'}})-[:hasEvent*]->(r:ReactionLikeEvent) WITH DISTINCT r "
         "MATCH (r)-[:input|output|catalystActivity|regulatedBy*1..2]->(x) WITH DISTINCT x "
         "MATCH path=(x)-[:hasComponent|hasMember|hasCandidate*1..3]->(m) WHERE x.stId IS NOT NULL AND m.stId IS NOT NULL "
         "RETURN DISTINCT x.stId + '|' + m.stId + '|' + reduce(t='', rel IN relationships(path) | t + type(rel) + ',')")
    leaves_of: dict[str, set[str]] = collections.defaultdict(set)
    set_pairs, complex_pairs = set(), set()
    for row in _cypher(q):
        if row.count("|") != 2:
            continue
        x, m, rels = row.split("|")
        leaves_of[x].add(m)
        (set_pairs if rels.startswith(("hasMember", "hasCandidate")) else complex_pairs).add((x, m))
    return leaves_of, complex_pairs, set_pairs


def fetch_direct_containers(pathway: str, hops: int, sources: str = "complex"):
    """(component_stid, container_stid) pairs within `hops` hasComponent steps.

    `sources` selects which components may be sources: "complex" (complex ->
    containing complex, what LNG_COMPOSITION_EDGES emits), "protein" (a protein
    -> the complexes containing it: the assembly direction extended beyond root
    complexes), or "all". The first version of this function returned "all"
    without saying so, and a Mitotic G1 simulation built on it was reported as
    validating the complex-only design; the real bundle recovered 3 routes
    where that simulation promised 259.
    """
    label = {"complex": "x:Complex", "protein": "x:EntityWithAccessionedSequence", "all": "x"}[sources]
    q = (f"MATCH (p:Pathway {{stId:'{pathway}'}})-[:hasEvent*]->(r:ReactionLikeEvent) WITH DISTINCT r "
         f"MATCH (r)-[:input|output|catalystActivity|regulatedBy*1..2]->({label}) WITH DISTINCT x "
         f"MATCH (y:Complex)-[:hasComponent*1..{hops}]->(x) WHERE y.stId IS NOT NULL AND x.stId IS NOT NULL "
         "RETURN DISTINCT x.stId + '|' + y.stId")
    return [tuple(row.split("|")) for row in _cypher(q) if row.count("|") == 1]


# ----------------------------------------------------------------------------
# bundle loading
# ----------------------------------------------------------------------------

def load_bundle(pathway_dir: Path):
    u2s = {r["uuid"]: r["stable_id"] for r in csv.DictReader((pathway_dir / "stid_to_uuid_mapping.csv").open())}
    s2u: dict[str, list[str]] = collections.defaultdict(list)
    for u, s in u2s.items():
        s2u[s].append(u)
    kind = {}
    if (pathway_dir / "nodes.csv").exists():
        kind = {r["uuid"]: r["node_kind"] for r in csv.DictReader((pathway_dir / "nodes.csv").open())}
    ladj: ADJ = collections.defaultdict(set)
    uadj: ADJ = collections.defaultdict(set)
    edge_uuids: dict[tuple[str, str], list[tuple[str, str]]] = collections.defaultdict(list)
    for e in csv.DictReader((pathway_dir / "logic_network.csv").open()):
        a, b = e["source_id"], e["target_id"]
        sa, sb = u2s.get(a, a), u2s.get(b, b)
        ladj[sa].add(sb); uadj[a].add(b); edge_uuids[(sa, sb)].append((a, b))
    cof: set[str] = set()
    if (pathway_dir / "cofactors.csv").exists():
        for r in csv.DictReader((pathway_dir / "cofactors.csv").open()):
            cof.update(v for v in r.values() if v and v.startswith("R-"))
    return u2s, s2u, kind, ladj, uadj, edge_uuids, cof


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def run(pathway: str, catalog: Path, simulate_hops: int = 0, simulate_sources: str = "complex") -> dict:
    rows = fetch_reaction_rows(pathway)
    adj, reaction_only, ents, rxns, role = build_reaction_graph(rows)
    leaves_of, complex_pairs, set_pairs = fetch_composition(pathway)
    add_composition(adj, complex_pairs, set_pairs, strict=True)

    d = next(p for p in catalog.iterdir() if p.is_dir() and p.name.endswith(pathway))
    u2s, s2u, kind, ladj, uadj, edge_uuids, cof = load_bundle(d)
    roots, terms = roots_and_terminals(reaction_only, ents, cof)
    # Small molecules (R-ALL) and declared cofactors neither start nor carry a
    # curator's route; reactions (R-HSA reaction ids are also R-HSA-) and
    # protein-bearing entities do.
    def carries(n: str) -> bool:
        return n.startswith("R-HSA-") and n not in cof
    adj = drop_carriers(adj, carries)
    reach = {r: signed_reach(adj, r) for r in roots}
    pairs = [(r, t, reach[r][t]) for r in roots for t in terms if t in reach[r]]

    present = set(ladj) | {t for ts in ladj.values() for t in ts}
    out = {"pathway": pathway, "reactions": len(rxns), "entities": len(ents),
           "roots": len(roots), "terminals": len(terms), "pairs": len(pairs),
           "pairs_pos": sum(1 for *_, s in pairs if s > 0), "pairs_neg": sum(1 for *_, s in pairs if s < 0)}
    absent_root = collections.Counter(); absent_term = collections.Counter()
    stid_ok = 0; uuid_ok = 0; present_pairs = 0; severed_by_root: collections.Counter = collections.Counter()
    sink_uuids = {u for u, k in kind.items() if k == "dissociation_sink"}
    sim_adj = None; sim_added = 0; sim_fan: list[int] = []
    if simulate_hops:
        direct = fetch_direct_containers(pathway, simulate_hops, simulate_sources)
        sim_adj, sim_added, sim_fan = simulate_composition(uadj, direct, s2u, sink_uuids)
    sim_ok = 0
    for r, t, _s in pairs:
        R, T = resolve_ids(r, present, leaves_of), resolve_ids(t, present, leaves_of)
        if not R:
            absent_root[r] += 1; continue
        if not T:
            absent_term[t] += 1; continue
        present_pairs += 1
        if any(reachable(ladj, {rr}, {tt}) for rr in R for tt in T):
            stid_ok += 1
        RU = {u for rr in R for u in s2u.get(rr, [])}; TU = {u for tt in T for u in s2u.get(tt, [])}
        if reachable(uadj, RU, TU):
            uuid_ok += 1
        else:
            severed_by_root[r] += 1
        if sim_adj is not None and reachable(sim_adj, RU, TU):
            sim_ok += 1
    # Classify the STILL-severed pairs by the first uuid break along the bundle's
    # stable-id path: which kind of node did the signal die on?
    break_kind: collections.Counter = collections.Counter()
    break_junction: collections.Counter = collections.Counter()
    src_adj = sim_adj if sim_adj is not None else uadj
    for r, t, _s in pairs:
        R, T = resolve_ids(r, present, leaves_of), resolve_ids(t, present, leaves_of)
        if not R or not T:
            continue
        RU = {u for rr in R for u in s2u.get(rr, [])}; TU = {u for tt in T for u in s2u.get(tt, [])}
        if reachable(src_adj, RU, TU):
            continue
        path = shortest_path(ladj, set(R), set(T))
        if not path:
            break_kind["no stable-id path in bundle"] += 1; continue
        hit = first_break(path, RU, edge_uuids)
        if hit is None:
            break_kind["unlocated"] += 1; continue
        step, frontier = hit
        fk = "/".join(sorted({kind.get(u, "?") for u in frontier})) or "?"
        break_kind[fk] += 1
        break_junction[(path[step - 1], path[step], fk)] += 1
    out.update({"absent_root_pairs": sum(absent_root.values()), "absent_roots": len(absent_root),
                "severed_break_kind": break_kind.most_common(6),
                "severed_top_junctions": [(a, b, k, n) for (a, b, k), n in break_junction.most_common(4)],
                "absent_root_roles": collections.Counter("/".join(sorted(role[r])) for r in absent_root),
                "absent_term_pairs": sum(absent_term.values()), "present_pairs": present_pairs,
                "stid_connected": stid_ok, "uuid_connected": uuid_ok, "uuid_severed": present_pairs - uuid_ok,
                "severed_by_root": severed_by_root.most_common(6)})
    if sim_adj is not None:
        fs = sorted(sim_fan)
        out.update({"sim_hops": simulate_hops, "sim_sources": simulate_sources, "sim_edges": sim_added,
                    "sim_fanout_median": fs[len(fs) // 2] if fs else 0, "sim_fanout_max": fs[-1] if fs else 0,
                    "sim_uuid_connected": sim_ok, "sim_uuid_severed": present_pairs - sim_ok})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pathway", help="Reactome pathway stId, e.g. R-HSA-909733")
    ap.add_argument("--catalog", type=Path, required=True, help="catalog root holding <...>/R-HSA-xxx bundles")
    ap.add_argument("--simulate-hops", type=int, default=0,
                    help="also recount with simulated composition edges (component node -> containing complex node) within N hasComponent hops")
    ap.add_argument("--simulate-sources", choices=("complex", "protein", "all"), default="complex",
                    help="which components may be sources in the simulation (default complex = what LNG emits)")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    res = run(a.pathway, a.catalog, a.simulate_hops, a.simulate_sources)
    if a.json:
        print(json.dumps(res, default=str)); return 0   # one line, so outputs can be appended as .jsonl
    print(f"{res['pathway']}: {res['reactions']} reactions, {res['entities']} entities, "
          f"{res['roots']} roots, {res['terminals']} terminals")
    print(f"  curator-reachable root->terminal pairs (STRICT): {res['pairs']}  (+{res['pairs_pos']} / -{res['pairs_neg']})")
    print(f"  root absent from bundle: {res['absent_root_pairs']} pairs over {res['absent_roots']} roots  {dict(res['absent_root_roles'])}")
    print(f"  terminal absent from bundle: {res['absent_term_pairs']} pairs")
    pp = res["present_pairs"]
    print(f"  stable-id level connected: {res['stid_connected']}/{pp}")
    print(f"  UUID level (solver's graph) connected: {res['uuid_connected']}/{pp}   severed: {res['uuid_severed']}")
    print(f"    severed by root: {res['severed_by_root']}")
    print(f"    still-severed by break kind: {res['severed_break_kind']}")
    print(f"    top junctions: {res['severed_top_junctions']}")
    if "sim_hops" in res:
        print(f"  SIMULATED {res['sim_sources']}->complex edges (<= {res['sim_hops']} hops): {res['sim_edges']} edges, "
              f"fan-out median {res['sim_fanout_median']} max {res['sim_fanout_max']}")
        print(f"    UUID level with them: {res['sim_uuid_connected']}/{pp}   severed: {res['sim_uuid_severed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
