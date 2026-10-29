"""specs/123 preparation: the part of a network that can carry a drug's effect to its readouts.

S = nodes reachable from the drug's sources (its pinned roots, or the sources of its cut edges) that also reach a
readout, plus the sources and the reached readouts. Every edge into S is kept; a parent outside S cannot change
under the drug (no route from the sources), so it is a boundary node pinned at its own baseline. The slim network
should solve to the same readout values as the full one; `compare` checks that, since loops that rail away from
baseline without any perturbation would break the assumption.
"""
import collections


def reach(adj, start):
    seen, stack = set(start), list(start)
    while stack:
        for v in adj.get(stack.pop(), ()):
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def slim(parsed, sources, readout_nodes):
    """(payload nodes, payload edges, boundary uuids) of the slim network."""
    fwd, back = collections.defaultdict(set), collections.defaultdict(set)
    for e in parsed["edges"]:
        p, c = str(e["parent_uuid"]), str(e["child_uuid"])
        fwd[p].add(c)
        back[c].add(p)
    down = reach(fwd, sources)
    up = reach(back, {r for r in readout_nodes if r in down})
    keep = (down & up) | set(sources)
    edges = [e for e in parsed["edges"] if str(e["child_uuid"]) in keep]
    boundary = {str(e["parent_uuid"]) for e in edges} - keep
    nodes = [n for n in parsed["nodes"] if str(n["uuid"]) in keep | boundary]
    return nodes, edges, boundary


def boundary_observations(parsed, boundary, confidence=1.0):
    """Each boundary node pinned at its baseline, on the UI scale (baseline x 100)."""
    base = {str(n["uuid"]): n["baseline"] for n in parsed["nodes"]}
    return {u: [100.0 * base[u], confidence] for u in sorted(boundary)}
