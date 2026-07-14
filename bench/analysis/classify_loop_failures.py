#!/usr/bin/env python3
"""Classify the loop topology behind specific readouts: is the readout inside
(or downstream of) a strongly-connected component, and is that SCC a NEGATIVE
feedback loop, a POSITIVE feedback loop, or a recycling/catalytic cycle?

Signed feedback: the net sign of a cycle = product of its edge signs
(pos=+1, neg=-1). Odd number of negative edges => net-negative feedback
(stabilizing / can invert direction on convergence); even => net-positive
(amplifying / bistable). We find one representative cycle through each SCC and
also report the internal edge-type / sign census.

Usage:
  python3 bench/analysis/classify_loop_failures.py <catalog_dir> <stid> [<stid>...]
  # stids given WITHOUT the R-HSA- prefix are auto-prefixed.
"""
import csv, os, sys
from collections import defaultdict
from pathlib import Path

CATALOG = Path(os.environ.get(
    "PATHWAY_CATALOG", str(Path.home() / "gitroot" / "logic-network-generator" / "output")))


def load(dirname):
    d = CATALOG / dirname
    stid = {}
    for r in csv.DictReader(open(d / "stid_to_uuid_mapping.csv")):
        stid[r["uuid"]] = r["stable_id"]
    adj = defaultdict(list)         # u -> [(v, sign, edge_type)]
    radj = defaultdict(list)        # reverse
    for r in csv.DictReader(open(d / "logic_network.csv")):
        s, t = r["source_id"], r["target_id"]
        sign = -1 if r["pos_neg"] == "neg" else 1
        adj[s].append((t, sign, r["edge_type"]))
        radj[t].append(s)
    return stid, adj, radj


def tarjan(adj, nodes):
    idx = {}; low = {}; on = set(); st = []; comp = {}; c = [0]; nc = [0]
    for s in nodes:
        if s in idx: continue
        work = [(s, iter(adj.get(s, [])))]
        idx[s] = low[s] = c[0]; c[0] += 1; st.append(s); on.add(s)
        while work:
            v, it = work[-1]; adv = False
            for (w, _sgn, _et) in it:
                if w not in idx:
                    idx[w] = low[w] = c[0]; c[0] += 1; st.append(w); on.add(w)
                    work.append((w, iter(adj.get(w, [])))); adv = True; break
                elif w in on:
                    low[v] = min(low[v], idx[w])
            if adv: continue
            work.pop()
            if work: low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == idx[v]:
                nc[0] += 1
                while True:
                    w = st.pop(); on.discard(w); comp[w] = nc[0]
                    if w == v: break
    return comp


def find_cycle_sign(adj, members):
    """DFS within `members` to find one cycle; return its net sign (product of
    edge signs) or None if no cycle found in the sampled DFS."""
    members = set(members)
    start = next(iter(members))
    stack = [(start, 1)]
    onpath = {start: 1}      # node -> cumulative sign from start
    visited = set()
    # iterative DFS tracking path sign
    def dfs(u, acc):
        visited.add(u)
        for (w, sgn, _et) in adj.get(u, []):
            if w not in members: continue
            if w in onpath:
                return acc * sgn // onpath[w] * onpath[w]  # cycle closed; net = product along cycle
            if w not in visited:
                onpath[w] = acc * sgn
                r = dfs(w, acc * sgn)
                if r is not None: return r
        onpath.pop(u, None)
        return None
    # simpler: since any cycle's sign is what we want, walk edges and detect back-edge sign product
    # recompute cleanly:
    onpath.clear(); visited.clear()
    sys.setrecursionlimit(100000)
    def dfs2(u, acc):
        visited.add(u); onpath[u] = acc
        for (w, sgn, _et) in adj.get(u, []):
            if w not in members: continue
            if w in onpath:                       # back edge -> cycle
                return (acc * sgn) * onpath[w]     # sign of loop = product around; onpath[w]=±1 prefix
            if w not in visited:
                r = dfs2(w, acc * sgn)
                if r is not None: return r
        onpath.pop(u, None)
        return None
    try:
        return dfs2(start, 1)
    except RecursionError:
        return None


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    dirname = sys.argv[1]
    stids = [s if s.startswith("R-HSA-") else f"R-HSA-{s}" for s in sys.argv[2:]]
    stid, adj, radj = load(dirname)
    nodes = set(adj) | {v for u in adj for (v, _s, _e) in adj[u]}
    comp = tarjan(adj, nodes)
    comp_size = defaultdict(int)
    comp_members = defaultdict(list)
    for u, cid in comp.items():
        comp_size[cid] += 1; comp_members[cid].append(u)
    by_stid = defaultdict(list)
    for u, s in stid.items(): by_stid[s].append(u)

    print(f"# {dirname.split('_R-HSA')[0]}")
    for s in stids:
        us = by_stid.get(s, [])
        print(f"\nreadout {s}: {len(us)} UUID(s)")
        seen_c = set()
        in_loop = False
        for u in us:
            cid = comp.get(u); sz = comp_size.get(cid, 1)
            if sz > 1 and cid not in seen_c:
                seen_c.add(cid); in_loop = True
                mem = comp_members[cid]
                memset = set(mem)
                negc = catc = tot = 0
                for x in mem:
                    for (v, sgn, et) in adj.get(x, []):
                        if v in memset:
                            tot += 1
                            if sgn < 0: negc += 1
                            if et == "catalyst": catc += 1
                cyc = find_cycle_sign(adj, mem)
                kind = ("NEGATIVE feedback" if cyc is not None and cyc < 0
                        else "positive/recycling" if cyc is not None
                        else "unknown")
                print(f"  in SCC #{cid}: size={sz}, internal edges={tot} "
                      f"(neg={negc}, catalyst={catc}); sample-cycle net sign="
                      f"{'-' if (cyc is not None and cyc<0) else '+' if cyc is not None else '?'} => {kind}")
        if not in_loop:
            # not in a loop; is it downstream of one? check 1-2 hops upstream
            up_in_loop = False
            for u in us:
                for p in radj.get(u, []):
                    if comp_size.get(comp.get(p), 1) > 1:
                        up_in_loop = True
            print(f"  NOT in a nontrivial SCC; immediate-upstream-in-loop={up_in_loop}")


if __name__ == "__main__":
    main()
