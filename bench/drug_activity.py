"""specs/122: a drug blocks its target's activity (DS_DRUG_ACTION=activity), instead of knocking the target out.

For target genes G, a `catalyst` or `regulator` edge from node S into reaction R is cut when S holds a protein of G
and the catalyst activity or regulation S stands for has an active unit holding a protein of G, or has none. A cut
edge is redirected to a new source node held at 0 (dropping it would let the reaction run as if it needed no
catalyst). Nothing is pinned for the drug itself. Genetic perturbations never use this.
"""
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from network_files import open_network  # noqa: E402

DRUG_ACTION = os.environ.get("DS_DRUG_ACTION", "knockout")
if DRUG_ACTION not in ("knockout", "activity"):
    raise SystemExit(f"DS_DRUG_ACTION={DRUG_ACTION!r}: must be knockout or activity (specs/122)")
CUT_TYPES = ("catalyst", "regulator")


def active_units():
    """{(reaction stId, entity stId): set of active-unit leaf stIds, or None when the activity has no active unit},
    for every human catalyst activity and regulation in the local Reactome."""
    from py2neo import Graph
    import benchmark_vs_mpbiopath as bm
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = graph.run(
        """
        MATCH (r:ReactionLikeEvent {speciesName: 'Homo sapiens'})-[:catalystActivity]->(a)-[:physicalEntity]->(pe)
        OPTIONAL MATCH (a)-[:activeUnit]->(au)
        RETURN r.stId AS r, pe.stId AS pe, COLLECT(DISTINCT au.stId) AS units
        UNION
        MATCH (r:ReactionLikeEvent {speciesName: 'Homo sapiens'})-[:regulatedBy]->(a)-[:regulator]->(pe)
        OPTIONAL MATCH (a)-[:activeUnit]->(au)
        RETURN r.stId AS r, pe.stId AS pe, COLLECT(DISTINCT au.stId) AS units
        """).data()
    units = {u for row in rows for u in row["units"]}
    import expression_gate
    below = expression_gate.entity_leaves(units)
    out = {}
    for row in rows:
        k = (row["r"], row["pe"])
        if not row["units"]:
            out.setdefault(k, None)
        else:
            leaves = set().union(*(below.get(u, {u}) for u in row["units"]))
            out[k] = leaves if out.get(k) is None else out[k] | leaves
    return out


def cut_edges(pathway_dir, target_leaves, units):
    """{(source uuid, target uuid)} of the catalyst / regulator edges a drug with these target leaves cuts."""
    pathway_dir = Path(pathway_dir)
    with open(pathway_dir / "nodes.csv") as f:
        node = {r["uuid"]: r for r in csv.DictReader(f)}
    out = set()
    with open_network(pathway_dir) as f:
        for e in csv.DictReader(f):
            if e.get("edge_type") not in CUT_TYPES:
                continue
            s = node.get(e["source_id"])
            if s is None:
                continue
            leaves = {x for x in (s.get("member_leaves") or "").split("|") if x}
            if not leaves & target_leaves:
                continue
            unit = units.get((e.get("edge_reaction_id"), s.get("diagram_entity_id")))
            if unit is None or unit & target_leaves:
                out.add((e["source_id"], e["target_id"]))
    return out


def redirected(parsed, cut, payload_extras):
    """(network payload, {block uuid: [0, 1]} observations): each cut edge's source replaced by its own block node."""
    by_uuid = {str(n["uuid"]): n for n in parsed["nodes"]}
    nodes, edges, obs = list(parsed["nodes"]), [], {}
    for e in parsed["edges"]:
        key = (str(e["parent_uuid"]), str(e["child_uuid"]))
        if key in cut:
            b = f"drugblock:{key[0]}:{key[1]}"
            if b not in obs:
                src = by_uuid[key[0]]
                nodes.append({"uuid": b, "name": "drug block", "reactome_id": None,
                              "entity_type": src.get("entity_type", ""), "baseline": src["baseline"],
                              "set_id": None})
                obs[b] = [0.0, 1.0]
            e = dict(e, parent_uuid=b)
        edges.append(e)
    return {"nodes": nodes, "edges": edges, "pathways": parsed["pathways"], **payload_extras}, obs


def payload_extras(parsed):
    """The tables that must travel with a POSTed network (benchmark_vs_mpbiopath.run_pathway's list)."""
    return {k: parsed.get(k) for k in ("containment", "cofactor_stids", "drug_stids", "pathogen_stids",
                                       "carrier_loops", "pools", "node_leaves")}


def perturb(bm, pathway_dir, parsed, genes, gene_to_stids, stid_to_uuids, units):
    """How a drug with target `genes` enters one network's solve, under DS_DRUG_ACTION.

    Returns (pins, observations, payload, reach_sources): `pins` {gene: [uuids]} as the knockout bench pins them
    (also the reachability comparator's sources in both modes, so the comparator does not move); observations for the
    solve; payload None to solve the cached network, else the redirected network to POST. Under `activity` nothing is
    pinned for the drug, so the target's own forms are scored."""
    g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
    if not any(g2u.values()):
        return {}, {}, None, set()
    pins, _ = bm.pin_occurrences(pathway_dir, g2u, gene_to_stids)
    roots = set().union(*pins.values()) if pins else set()
    if DRUG_ACTION == "knockout":
        return pins, {u: [bm.PERTURB_UI_DOWN, bm.PIN_CONFIDENCE] for u in sorted(roots)}, None, roots
    leaves = set().union(*(set(gene_to_stids.get(g, [])) for g in genes))
    cut = cut_edges(pathway_dir, leaves, units)
    if not cut:
        return pins, {}, None, roots
    payload, obs = redirected(parsed, cut, payload_extras(parsed))
    return pins, {u: [v, bm.PIN_CONFIDENCE] for u, (v, _) in obs.items()}, payload, roots
