#!/usr/bin/env python3
"""Is the readout reachable from the perturbed gene IN REACTOME ITSELF?

The benchmark's `no_path` bucket says the logic network has no route from the
perturbed gene to the readout. It does not say WHY, and the two possibilities
need completely different work:

  the generator lost a connection curators made   -> fixable upstream
  Reactome has no directed route either           -> not fixable by connectivity;
                                                     the curator's call rests on
                                                     something the pathway graph
                                                     does not encode

So this asks the database the same question the logic network was asked. It
walks Reactome's own bipartite structure — an entity feeds a ReactionLikeEvent
as input, catalyst or regulator; that event produces output entities; those
feed further events — restricted to events inside the pathway, which is the
same scope the logic network has.

Containment is followed in both directions (a complex's components, a set's
members, and their parents), because a gene's protein usually participates as
part of something larger and the readout is often a complex.
"""

from __future__ import annotations

import collections
import os

# Connection details come from _common, which reads the environment and
# refuses to run without NEO4J_PASSWORD rather than carrying a literal.
from _common import NEO4J_PASSWORD, NEO4J_URL, NEO4J_USER, graph  # noqa: F401


PARTICIPATION = """
MATCH (p:Pathway {stId:$pid})-[:hasEvent*1..6]->(rle:ReactionLikeEvent)
OPTIONAL MATCH (rle)-[:input]->(i:PhysicalEntity)
OPTIONAL MATCH (rle)-[:output]->(o:PhysicalEntity)
OPTIONAL MATCH (rle)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(c:PhysicalEntity)
OPTIONAL MATCH (rle)-[:regulatedBy]->(:Regulation)-[:regulator]->(g:PhysicalEntity)
RETURN rle.stId AS rle,
       COLLECT(DISTINCT i.stId) AS inputs,
       COLLECT(DISTINCT o.stId) AS outputs,
       COLLECT(DISTINCT c.stId) AS catalysts,
       COLLECT(DISTINCT g.stId) AS regulators
"""

# Containment among the entities the pathway actually touches.
CONTAINMENT = """
MATCH (p:Pathway {stId:$pid})-[:hasEvent*1..6]->(:ReactionLikeEvent)
      -[:input|output|catalystActivity|regulatedBy|physicalEntity|regulator*1..2]->(pe:PhysicalEntity)
OPTIONAL MATCH (pe)-[:hasComponent|hasMember|hasCandidate*1..3]->(child:PhysicalEntity)
RETURN DISTINCT pe.stId AS parent, COLLECT(DISTINCT child.stId) AS children
"""


def build_pathway_graph(pathway_stid: str):
    """entity -> entities reachable downstream, through the pathway's events."""
    rows = graph().run(PARTICIPATION, pid=pathway_stid).data()
    feeds: dict[str, set[str]] = collections.defaultdict(set)   # entity -> reactions
    produces: dict[str, set[str]] = collections.defaultdict(set)  # reaction -> entities
    entities: set[str] = set()
    for r in rows:
        rle = r["rle"]
        for group in ("inputs", "catalysts", "regulators"):
            for e in r[group] or []:
                if e:
                    feeds[e].add(rle)
                    entities.add(e)
        for e in r["outputs"] or []:
            if e:
                produces[rle].add(e)
                entities.add(e)

    contains: dict[str, set[str]] = collections.defaultdict(set)
    for r in graph().run(CONTAINMENT, pid=pathway_stid).data():
        parent = r["parent"]
        if not parent:
            continue
        entities.add(parent)
        for child in r["children"] or []:
            if child:
                contains[parent].add(child)
                entities.add(child)
    return feeds, produces, contains, entities


def reachable(feeds, produces, contains, sources: set[str],
              max_steps: int = 200000, containment: str = "both") -> set[str]:
    """Entities downstream of `sources`, following containment both ways.

    Containment is bidirectional on purpose: a perturbed protein acts through
    the complex that contains it, and a complex's behaviour is attributed to
    its members. That is exactly what the logic network's decomposition does.
    """
    # `containment` controls how permissive the walk is, and it matters:
    #   "none"  pure reaction chaining. Strictest.
    #   "down"  a complex's state is attributed to its members.
    #   "both"  also lets a member act through its parent — which additionally
    #           permits parent-then-sibling, i.e. CO-MEMBERSHIP rather than
    #           causation. That is the "shared ancestor" pattern, and counting
    #           it as a route would overstate how much the generator lost.
    up = collections.defaultdict(set)
    if containment == "both":
        for parent, children in contains.items():
            for child in children:
                up[child].add(parent)
    if containment == "none":
        contains = {}

    seen = set(sources)
    queue = collections.deque(sources)
    steps = 0
    while queue and steps < max_steps:
        steps += 1
        entity = queue.popleft()
        nxt = set(contains.get(entity, ())) | set(up.get(entity, ()))
        for rle in feeds.get(entity, ()):
            nxt |= produces.get(rle, set())
        for e in nxt:
            if e not in seen:
                seen.add(e)
                queue.append(e)
    return seen
