#!/usr/bin/env python3
"""How often do the benchmark's own pins meet the two pool traps of specs/052 R2?

The cycle-balance pool (specs/039, DS_CYCLE_MODE=balance) takes its supply
from the producers of its states that the pool does not itself reach, and
that reach search stops at every pinned node. So a pin on a pool's own
recycling path turns a recycling step into "external supply" (trap B,
the PTK6 case). Separately, a pinned pool intermediate or step copy sends
the whole pool to the iteration fallback (trap A).

For every scored case in a results directory, this reads the pinned uuids the
bench logged (`gene_uuids`) and asks, structurally, whether they hit either
trap. Nothing is solved: a cut producer reading baseline changes nothing, so
trap B counts are an upper bound on cases whose prediction can move.

Mirrors the solver (src/solvers/steady_state.jl, pool setup):
- reach follows activator (positive) edges only, inside the pool's strongly
  connected component, and never enters a pinned node;
- the pins of a solve include the inert holds: cofactors (DS_COFACTOR_MODE=
  inert) and pathogen nodes (DS_PATHOGEN_MODE=inert), so those are barriers
  in the control too, and a gene pin only counts for what it adds;
- a pool is solved only if all its nodes sit in one cyclic component.
Supply producers are the positive in-edges of a state from outside the pool's
own nodes and step copies, which is how index_pools builds supply_rxns.

Usage: pool_pin_census.py BUILD_DIR RESULTS_DIR
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from network_files import open_network  # noqa: E402


def scc(nodes, adj):
    """Iterative Tarjan: node -> component id."""
    index, low, comp, on, stack, out = {}, {}, {}, set(), [], 0
    counter = 0
    for root in nodes:
        if root in index:
            continue
        work = [(root, iter(adj.get(root, ())))]
        index[root] = low[root] = counter; counter += 1
        stack.append(root); on.add(root)
        while work:
            v, it = work[-1]
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter; counter += 1
                    stack.append(w); on.add(w)
                    work.append((w, iter(adj.get(w, ()))))
                    break
                if w in on:
                    low[v] = min(low[v], index[w])
            else:
                work.pop()
                if work:
                    low[work[-1][0]] = min(low[work[-1][0]], low[v])
                if low[v] == index[v]:
                    while True:
                        w = stack.pop(); on.discard(w); comp[w] = out
                        if w == v:
                            break
                    out += 1
    return comp


def load_pathway(d):
    adj, inn, nodes = collections.defaultdict(set), collections.defaultdict(set), set()
    with open_network(d) as f:
        for r in csv.DictReader(f):
            s, t = r["source_id"], r["target_id"]
            nodes.update((s, t))
            if r["pos_neg"] == "pos":
                adj[s].add(t); inn[t].add(s)
    comp = scc(sorted(nodes), adj)
    size = collections.Counter(comp.values())
    stid_uuids = collections.defaultdict(set)
    with open(d / "stid_to_uuid_mapping.csv") as f:
        rd = csv.reader(f); next(rd)
        for row in rd:
            stid_uuids[row[1]].add(row[0])
    holds = set()
    for name in ("cofactors.csv", "pathogens.csv"):
        p = d / name
        if p.exists():
            with open(p) as f:
                for r in csv.DictReader(f):
                    if name != "cofactors.csv" or r.get("in_network", "1") != "0":
                        holds |= stid_uuids.get(r["stable_id"], set())
    pools = collections.defaultdict(lambda: {"states": set(), "inter": set(), "steps": set()})
    with open(d / "pools.csv") as f:
        for r in csv.DictReader(f):
            pools[r["pool_id"]]["states" if r["role"] == "state" else "inter"].add(r["node_uuid"])
    with open(d / "pool_transitions.csv") as f:
        for r in csv.DictReader(f):
            pools[r["pool_id"]]["steps"].add(r["reaction_uuid"])
    solved = {}
    for pid, p in pools.items():
        allp = p["states"] | p["inter"] | p["steps"]
        cs = {comp.get(u) for u in allp}
        if len(cs) != 1 or None in cs or size[next(iter(cs))] < 2:
            continue
        c = next(iter(cs))
        prods = {a for s in p["states"] for a in inn[s] if a not in allp}
        solved[pid] = dict(p, comp=c, prods=prods)
    return adj, comp, holds, solved


def reach(p, adj, comp, barrier):
    seen, stack = set(p["states"]), list(p["states"])
    while stack:
        v = stack.pop()
        for w in adj.get(v, ()):
            if w in seen or w in barrier or comp.get(w) != p["comp"]:
                continue
            seen.add(w); stack.append(w)
    return seen


def main(build, results):
    build, results = Path(build), Path(results)
    owner = {}
    for d in sorted(p for p in build.iterdir() if (p / "pools.csv").exists()):
        with open(d / "nodes.csv", newline="") as f:
            for r in csv.DictReader(f):
                owner[r["uuid"]] = d
    cache = {}
    tally = collections.Counter()
    per_path = collections.defaultdict(collections.Counter)
    for axis in ("curator", "experimental"):
        with open(results / f"{axis}_cases.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["valid"] != "1":
                    continue
                pins = set(filter(None, r["gene_uuids"].split("|")))
                dirs = {owner[u] for u in pins if u in owner}
                tally[(axis, "cases")] += 1
                if not dirs:
                    continue
                d = dirs.pop()
                if d not in cache:
                    cache[d] = load_pathway(d)
                adj, comp, holds, solved = cache[d]
                if not solved:
                    continue
                tally[(axis, "cases in a pathway with a solved pool")] += 1
                trap_a = trap_b = state_pin = False
                for p in solved.values():
                    if pins & (p["inter"] | p["steps"]):
                        trap_a = True
                    if pins & p["states"]:
                        state_pin = True
                    r0 = reach(p, adj, comp, holds)
                    r1 = reach(p, adj, comp, holds | pins)
                    if (p["prods"] & r0) - r1:
                        trap_b = True
                ok = r["predicted"] == r["expected"]
                for flag, name in ((trap_a, "A: pin on a pool intermediate/step (fallback)"),
                                   (trap_b, "B: pin cuts a pool recycling path"),
                                   (state_pin, "pin on a pool state (handled)")):
                    if flag:
                        tally[(axis, name)] += 1
                        tally[(axis, name + ", correct")] += ok
                        per_path[(axis, name)][r["pathway"]] += 1
    for axis in ("curator", "experimental"):
        print(f"== {axis}")
        for k in sorted(k for k in tally if k[0] == axis):
            print(f"  {k[1]:60s} {tally[k]}")
        for k in sorted(k for k in per_path if k[0] == axis):
            top = ", ".join(f"{p} {n}" for p, n in per_path[k].most_common(6))
            print(f"  [{k[1][:2]}] by pathway: {top}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
