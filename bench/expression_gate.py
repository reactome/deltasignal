"""specs/121 phase 2: the nodes a cell line cannot form, from the genes it does not express.

A node is pinned at baseline in a cell line when a protein it REQUIRES is absent there:
- simple_complex, dissociation_sink: every leaf is required;
- set_variant (one member combination, specs/046): the leaves that are not members of any of its source sets, plus
  its chosen members (member_leaves lists every member of the sets, chosen or not);
- simple_entity with one leaf: that leaf; with several (a set kept as one node): only if every leaf is absent;
- set_pool, variant_pool (OR pools), reactions: never.
A leaf is absent when its gene's expression is below the threshold. Leaves of small molecules or genes without an
expression value are present.
"""
import csv
from pathlib import Path

NEVER = {"set_pool", "variant_pool", "reaction"}


def entity_leaves(stids):
    """{stId: {leaf stIds below it}} for sets and complexes (members, candidates and components, recursively), from
    the local Reactome Neo4j. A leaf (protein, small molecule) has no entry: it is its own leaf."""
    from py2neo import Graph
    import benchmark_vs_mpbiopath as bm
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $ids AS sid
        MATCH (s:PhysicalEntity {stId: sid})-[:hasMember|hasCandidate|hasComponent*1..]->(leaf)
        WHERE NOT leaf:Complex AND NOT leaf:EntitySet
        RETURN sid AS id, COLLECT(DISTINCT leaf.stId) AS leaves
        """,
        ids=sorted(stids),
    ).data()
    return {r["id"]: set(r["leaves"]) for r in rows}


def required_leaves(row, below):
    """(kind of test, leaves): ("any", L) pins when any of L is absent, ("all", L) when every one is.
    `below` is entity_leaves() for the node's source sets and chosen members."""
    kind = row.get("node_kind") or ""
    leaves = {x for x in (row.get("member_leaves") or "").split("|") if x}
    if kind in NEVER or not leaves:
        return "any", set()
    if kind == "set_variant":
        in_sets = set().union(*(below.get(s, set()) for s in (row.get("source_sets") or "").split("|") if s))
        chosen = set().union(*(below.get(c, {c}) for c in (row.get("chosen_members") or "").split("|") if c))
        return "any", (leaves - in_sets) | (leaves & chosen)
    if kind == "simple_entity" and len(leaves) > 1:
        return "all", leaves
    return "any", leaves


def gated_nodes(pathway_dir, absent_leaves, below):
    """Node uuids a cell line lacking `absent_leaves` cannot form, in one pathway."""
    out = set()
    with open(Path(pathway_dir) / "nodes.csv") as f:
        for row in csv.DictReader(f):
            how, req = required_leaves(row, below)
            if not req:
                continue
            hit = req <= absent_leaves if how == "all" else bool(req & absent_leaves)
            if hit:
                out.add(row["uuid"])
    return out


def absent_genes(expression_csv, model_id, genes, threshold=1.0):
    """Genes (of `genes`) whose log2(TPM+1) in `model_id` is below `threshold`; genes not in the table are left out."""
    with open(expression_csv) as f:
        header = f.readline().rstrip("\n").split(",")
        col = {h.split(" (")[0]: i for i, h in enumerate(header) if i}
        for line in f:
            if line.split(",", 1)[0] != model_id:
                continue
            vals = line.rstrip("\n").split(",")
            return {g for g in genes if g in col and vals[col[g]] != "" and float(vals[col[g]]) < threshold}
    return None
