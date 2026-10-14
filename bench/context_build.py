#!/usr/bin/env python3
"""Build a "Reactome context" catalog from a built catalog (specs/056).

Each pathway P becomes P plus its neighbours, joined by directed bridge edges:
- an entity is identified by its full `nodes.csv` identity (every column but
  the uuid), and is a HUB when its diagram entity is in more than `--hub`
  of the build's pathways;
- cofactors, drugs, pathogens, reactions and set pools are never bridged;
- a BRIDGE runs from a copy that has an in-edge in its own pathway (produced
  there) to a copy of the same identity, in a different pathway, that has no
  in-edge in its own pathway (a root there). Only bridges between P and one
  neighbour are added (a star around P, never neighbour to neighbour);
- the neighbours of P are the pathways with at least one such bridge to or
  from P;
- a neighbour keeps only what can reach P: the ancestors, inside it, of its
  bridges into P, closed over its pools and carrier loops (a pool or loop
  with one kept node is kept whole, with its ancestors). P's values depend
  only on their ancestors, so the pruning is exact; `--no-prune` keeps whole
  neighbours, to check that.
A bridge is `pos,or,context_bridge`: a root fed by several producers reads
their mean. No uuid is merged and nothing inside a pathway is bridged.

Every side table is the union of the included pathways' tables, with pool and
loop ids prefixed by their pathway. `context_own.csv` lists P's own uuids, so
the bench scores P's readouts on P's own nodes only.

Usage: context_build.py BASE_BUILD OUT_BUILD [--hub 5] [--census]
"""
import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

KEY = ("node_kind", "diagram_entity_id", "compartment", "member_leaves", "source_sets", "chosen_members")
NEVER = {"reaction", "set_pool"}
PREFIXED = {"pools.csv": ("pool_id",), "pool_carriers.csv": ("pool_id",),
            "pool_transitions.csv": ("pool_id", "path_id"), "carrier_loops.csv": ("loop_id",)}
UNION = ("nodes.csv", "stid_to_uuid_mapping.csv", "entity_reaction_proxy_mapping.csv",
         "containment.csv", "cofactors.csv", "drugs.csv", "pathogens.csv", "carrier_loops.csv",
         "pools.csv", "pool_carriers.csv", "pool_transitions.csv", "boundary_edges.csv")
EDGE_HEAD = ["source_id", "target_id", "pos_neg", "and_or", "edge_type", "stoichiometry", "edge_reaction_id"]


def rows(p):
    if not p.exists():
        return None, []
    with open(p, newline="") as f:
        r = csv.DictReader(f)
        return r.fieldnames, list(r)


def load(d):
    _, nodes = rows(d / "nodes.csv")
    held = set()
    for t in ("cofactors.csv", "drugs.csv", "pathogens.csv"):
        held |= {r["stable_id"] for r in rows(d / t)[1]}
    indeg = collections.Counter()
    for t in ("logic_network.csv", "boundary_edges.csv"):
        for r in rows(d / t)[1]:
            indeg[r["target_id"]] += 1
    ident = {}
    for r in nodes:
        if r["node_kind"] in NEVER or r["diagram_entity_id"] in held:
            continue
        ident[r["uuid"]] = tuple(r.get(k, "") for k in KEY)
    return {"ident": ident, "indeg": indeg}


def plan(base, hub):
    paths = sorted(p for p in base.iterdir() if (p / "nodes.csv").exists())
    info = {p.name: load(p) for p in paths}
    n_paths = collections.Counter()
    for v in info.values():
        for de in {k[1] for k in v["ident"].values()}:
            n_paths[de] += 1
    by_key = collections.defaultdict(lambda: collections.defaultdict(lambda: ([], [])))
    for name, v in info.items():
        for u, k in v["ident"].items():
            if n_paths[k[1]] > hub:
                continue
            by_key[k][name][0 if v["indeg"][u] else 1].append(u)
    bridges = collections.defaultdict(list)          # (A, B) -> [(src in A, dst in B)]
    for k, per in by_key.items():
        if len(per) < 2:
            continue
        for a, (prod, _) in per.items():
            for b, (_, root) in per.items():
                if a != b:
                    bridges[(a, b)] += [(s, t) for s in prod for t in root]
    return [p.name for p in paths], bridges, n_paths


def kept_nodes(qd, srcs):
    """The ancestors of `srcs` in pathway dir `qd`, closed over its pools and
    carrier loops."""
    inc = collections.defaultdict(set)
    for t in ("logic_network.csv", "boundary_edges.csv"):
        for r in rows(qd / t)[1]:
            inc[r["target_id"]].add(r["source_id"])
    groups = collections.defaultdict(set)
    for t, cols in (("pools.csv", ("node_uuid",)), ("pool_carriers.csv", ("carrier_uuid", "release_reaction_uuid")),
                    ("pool_transitions.csv", ("source_uuid", "target_uuid", "reaction_uuid")),
                    ("carrier_loops.csv", ("free_uuid", "loop_source_uuid"))):
        for r in rows(qd / t)[1]:
            gid = (t.split("_")[0], r.get("pool_id") or r.get("loop_id"))
            groups[gid] |= {r[c] for c in cols if r.get(c)}
    member_of = collections.defaultdict(set)
    for g, us in groups.items():
        for u in us:
            member_of[u].add(g)
    seen, stack, done = set(), list(srcs), set()
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u)
        stack.extend(inc.get(u, ()))
        for g in member_of.get(u, ()):
            if g not in done:
                done.add(g); stack.extend(groups[g])
    return seen


def keep_row(t, r, kept):
    if t in ("nodes.csv", "stid_to_uuid_mapping.csv"):
        return r["uuid"] in kept
    if t == "entity_reaction_proxy_mapping.csv":
        return r["proxy_uuid"] in kept
    if t in ("logic_network.csv", "boundary_edges.csv"):
        return r["target_id"] in kept
    if t == "pools.csv":
        return r["node_uuid"] in kept
    if t == "pool_carriers.csv":
        return r["carrier_uuid"] in kept
    if t == "pool_transitions.csv":
        return r["source_uuid"] in kept
    if t == "carrier_loops.csv":
        return r["free_uuid"] in kept
    return True


def build_one(base, out, P, nbrs, bridges, prune=True):
    d = out / P
    d.mkdir(parents=True)
    members = [P] + nbrs
    kept = {P: None}
    for q in nbrs:
        kept[q] = kept_nodes(base / q, [s for s, _ in bridges.get((q, P), ())]) if prune else None
    for t in UNION:
        head, allr = None, []
        for q in members:
            h, rs = rows(base / q / t)
            if h is None:
                continue
            head = head or h
            for r in rs:
                if kept[q] is not None and not keep_row(t, r, kept[q]):
                    continue
                for c in PREFIXED.get(t, ()):
                    r[c] = f"{q}:{r[c]}"
                allr.append(r)
        if head is None:
            continue
        if t in ("cofactors.csv", "drugs.csv", "pathogens.csv", "containment.csv"):
            allr = list({tuple(r[c] for c in head): r for r in allr}.values())
        with open(d / t, "w", newline="") as f:
            w = csv.DictWriter(f, head, lineterminator="\n")
            w.writeheader(); w.writerows(allr)
    with open(d / "logic_network.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(EDGE_HEAD)
        for q in members:
            for r in rows(base / q / "logic_network.csv")[1]:
                if kept[q] is None or keep_row("logic_network.csv", r, kept[q]):
                    w.writerow([r[c] for c in EDGE_HEAD])
        n = 0
        for q in nbrs:
            for a, b in ((P, q), (q, P)):
                for s, t in bridges.get((a, b), ()):
                    if kept[q] is not None and q == b and t not in kept[q]:
                        continue
                    w.writerow([s, t, "pos", "or", "context_bridge", "1", ""])
                    n += 1
    with open(d / "context_own.csv", "w", newline="") as f:
        f.write("uuid\n" + "".join(r["uuid"] + "\n" for r in rows(base / P / "nodes.csv")[1]))
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("base"); ap.add_argument("out", nargs="?")
    ap.add_argument("--hub", type=int, default=5)
    ap.add_argument("--census", action="store_true")
    ap.add_argument("--no-prune", dest="prune", action="store_false")
    ap.add_argument("--only", action="append", default=[], help="build only these pathways (checks)")
    a = ap.parse_args(argv)
    base = Path(a.base).resolve()
    names, bridges, n_paths = plan(base, a.hub)
    nbr = collections.defaultdict(set)
    for (x, y), bs in bridges.items():
        if bs:
            nbr[x].add(y); nbr[y].add(x)
    size = {p: sum(1 for _ in open(base / p / "nodes.csv")) - 1 for p in names}
    print(f"pathways {len(names)}; hub cutoff > {a.hub} pathways; "
          f"hub entities {sum(1 for n in n_paths.values() if n > a.hub)} of {len(n_paths)}")
    print(f"bridge pairs (directed) {sum(1 for b in bridges.values() if b)}; "
          f"bridges {sum(len(b) for b in bridges.values())}")
    print("pathway\town_nodes\tneighbours\tcontext_nodes\tbridges_in_star")
    for p in names:
        ns = sorted(nbr[p])
        nb = sum(len(bridges.get((p, q), ())) + len(bridges.get((q, p), ())) for q in ns)
        print(f"{p}\t{size[p]}\t{len(ns)}\t{size[p] + sum(size[q] for q in ns)}\t{nb}")
    if a.census:
        return
    out = Path(a.out).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    out.mkdir(parents=True)
    tot = {}
    for p in a.only or names:
        tot[p] = {"neighbours": sorted(nbr[p]),
                  "bridges": build_one(base, out, p, sorted(nbr[p]), bridges, a.prune)}
    b = json.load(open(base / "BUILD.json"))
    repo = Path(__file__).resolve().parent.parent
    try:
        sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True).stdout.strip() or "unknown"
        dirty = bool(subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--",
                                     "bench/context_build.py"], capture_output=True, text=True).stdout.strip())
    except OSError:                      # no git (the test container)
        sha, dirty = "unknown", False
    b.update({"build_id": out.name, "context_builder_commit": sha + ("-dirty" if dirty else ""), "context_of": base.name, "context_hub": a.hub,
              "context_pruned": a.prune, "context": tot,
              "pathways_built": len(tot), "status": "complete" if not a.only else "partial"})
    json.dump(b, open(out / "BUILD.json", "w"), indent=1)


if __name__ == "__main__":
    main()
