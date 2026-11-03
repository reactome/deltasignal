"""specs/122: a drug blocks its target's activity (DS_DRUG_ACTION=activity), instead of knocking the target out.

For target genes G, a `catalyst` or `regulator` edge from node S into reaction R is cut when S holds a protein of G
and the catalyst activity or regulation S stands for has an active unit holding a protein of G, or has none. A cut
edge is redirected to a new source node held at 0 (dropping it would let the reaction run as if it needed no
catalyst). Nothing is pinned for the drug itself. Genetic perturbations never use this.

specs/123: DS_DRUG_ACTION=mechanism picks knockout or activity per drug from its documented mechanism, read from the
label table DS_DRUG_MECHANISM_TABLE (default drug_mechanism.tsv beside this file: axis, drug, label, ...); a drug
labelled `unknown` takes DS_DRUG_MECHANISM_UNKNOWN (default knockout).
"""
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from network_files import open_network  # noqa: E402

DRUG_ACTION = os.environ.get("DS_DRUG_ACTION", "knockout")
if DRUG_ACTION not in ("knockout", "activity", "mechanism"):
    raise SystemExit(f"DS_DRUG_ACTION={DRUG_ACTION!r}: must be knockout, activity (specs/122) or mechanism (specs/123)")
MECHANISM_TABLE = Path(os.environ.get("DS_DRUG_MECHANISM_TABLE", Path(__file__).resolve().with_name("drug_mechanism.tsv")))
MECHANISM_UNKNOWN = os.environ.get("DS_DRUG_MECHANISM_UNKNOWN", "knockout")
if MECHANISM_UNKNOWN not in ("knockout", "activity"):
    raise SystemExit(f"DS_DRUG_MECHANISM_UNKNOWN={MECHANISM_UNKNOWN!r}: must be knockout or activity (specs/123)")
CUT_TYPES = ("catalyst", "regulator")


def needs_units():
    """Whether some drug may be modelled as an activity inhibitor (active_units is a Neo4j query)."""
    return DRUG_ACTION != "knockout"


def describe():
    """DS_DRUG_ACTION for a Protocol line; under `mechanism`, the label table and its checksum too."""
    if DRUG_ACTION != "mechanism":
        return DRUG_ACTION
    import hashlib
    digest = hashlib.sha256(MECHANISM_TABLE.read_bytes()).hexdigest()[:12]
    return f"mechanism (labels {MECHANISM_TABLE.name} sha256 {digest}, unknown={MECHANISM_UNKNOWN})"


def drug_actions(axis, drugs, table=None):
    """{drug: knockout | activity} for one bench's drugs under DS_DRUG_ACTION; under `mechanism`, each drug's label
    from the table (rows of this axis), `unknown` as DS_DRUG_MECHANISM_UNKNOWN. A drug missing from the table is an
    error, so an arm cannot silently fall back."""
    if DRUG_ACTION != "mechanism":
        return {d: DRUG_ACTION for d in drugs}
    labels = {}
    with open(table or MECHANISM_TABLE) as f:
        for r in csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"):
            if r["axis"] == axis:
                if r["label"] not in ("knockout", "activity", "unknown"):
                    raise SystemExit(f"{table or MECHANISM_TABLE}: {axis} {r['drug']}: label {r['label']!r}")
                labels[r["drug"]] = r["label"]
    missing = sorted(set(drugs) - set(labels))
    if missing:
        raise SystemExit(f"DS_DRUG_ACTION=mechanism: no {axis} label for {missing} in {table or MECHANISM_TABLE}")
    return {d: MECHANISM_UNKNOWN if labels[d] == "unknown" else labels[d] for d in drugs}


def gene_actions(members):
    """{gene: action} for a drug set given [(action, genes)] per member drug; a gene some member knocks out is
    knocked out (a pin at 0 already removes what the cut would)."""
    out = {}
    for action, genes in members:
        for g in genes:
            out[g] = "knockout" if "knockout" in (action, out.get(g)) else action
    return out


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


def perturb(bm, pathway_dir, parsed, genes, gene_to_stids, stid_to_uuids, units, action=None):
    """How a drug with target `genes` enters one network's solve: `action` is knockout or activity for every gene, or
    {gene: action} (specs/123); None means DS_DRUG_ACTION, which must then be knockout or activity.

    Returns (pins, observations, edge_redirects, reach_sources): `pins` {gene: [uuids]} as the knockout bench pins
    them (also the reachability comparator's sources under every action, so the comparator does not move);
    observations for the solve, whose keys are the nodes actually pinned (the knocked-out genes' roots); the edges to
    cut for the activity genes, sent with the cached network's id (the server redirects each to a block node it pins
    at 0; `redirected` builds the same network client-side, which test_api_errors.jl shows solves identically).
    An activity gene is not pinned, so its own forms are scored."""
    g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
    if not any(g2u.values()):
        return {}, {}, set(), set()
    action = DRUG_ACTION if action is None else action
    by_gene = action if isinstance(action, dict) else {g: action for g in genes}
    if set(by_gene.values()) - {"knockout", "activity"} or set(genes) - set(by_gene):
        raise ValueError(f"perturb needs knockout or activity for every gene, got {by_gene}")
    pins, _ = bm.pin_occurrences(pathway_dir, g2u, gene_to_stids)
    roots = set().union(*pins.values()) if pins else set()
    held = {u for g, us in pins.items() if by_gene[g] == "knockout" for u in us}
    obs = {u: [bm.PERTURB_UI_DOWN, bm.PIN_CONFIDENCE] for u in sorted(held)}
    leaves = set().union(*(set(gene_to_stids.get(g, [])) for g in genes if by_gene[g] == "activity"))
    return pins, obs, cut_edges(pathway_dir, leaves, units) if leaves else set(), roots
