#!/usr/bin/env python3
"""Build a "Reactome context" catalog from a built catalog (specs/056).

Each pathway P becomes P plus its neighbours, joined by directed bridge edges:
- an entity is identified by its full `nodes.csv` identity (every column but
  the uuid), and is a HUB when its diagram entity is in more than `--hub`
  of the build's pathways;
- cofactors, drugs, pathogens, reactions and set pools are never bridged;
- a neighbour Q first loses every reaction P already has (its reaction
  nodes, their edges and the edges attributed to it), so it contributes only
  what P lacks; a parent pathway and its sub-pathway do not re-link P's own
  reactions to their copies;
- a BRIDGE runs from a copy that has an in-edge in its own (reduced) pathway
  (produced there) to a copy of the same identity, in the other pathway, that
  has no in-edge there (a root). Only bridges between P and one neighbour
  are added (a star around P, never neighbour to neighbour);
- the neighbours of P are the pathways with at least one such bridge INTO
  P (one that only receives from P cannot change P's values);
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

Under `--feed-in`, no bridge leaves P: neighbours only feed P's roots, and
the context cannot close a cycle through P.

Under `--curated LINKS S_Q` (specs/132) bridges come only from Reactome's curated handoffs, never from identity:
LINKS (bench/analysis/curated_links_132.py) lists pairs (r in P, s outside P, entity e an output of s and an input,
catalyst or regulator of r); S_Q names the pathway holding s (its lowest ancestor with a diagram). A bridge runs from
Q's copy of e that s outputs to P's copy of e that feeds r and is a root in P (a root as the bench pins it,
`benchmark_vs_mpbiopath.root_indegree`: depletion and recycle in-edges do not count). Feed-in only, no hub rule, centres
limited by `--centres` (the scored pathways); `--l2` keeps only pairs whose s lies in a pathway P's diagrams
encapsulate. With an empty LINKS it builds each centre alone (the control). `--whole` (specs/133 amendment 2): a source
pathway is neither reduced nor pruned, so it solves as it does alone.

Under `--curated LINKS S_Q --transitive` (specs/133) handoffs are followed back as far as they go, as a tree per centre
(`plan_transitive`); LINKS comes from bench/analysis/closure_133.py, whose `extracellular` handoffs are never bridged
(under either planner). `--l2` applies to both planners. BUILD.json records both filters and the root rule.

Usage: context_build.py BASE_BUILD OUT_BUILD [--hub 5] [--feed-in] [--census]
       context_build.py BASE_BUILD OUT_BUILD --curated LINKS S_Q --centres LIST [--l2] [--transitive] [--whole]
"""
import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402  (the bench's root definition)

KEY = ("node_kind", "diagram_entity_id", "compartment", "member_leaves", "source_sets", "chosen_members")
NEVER = {"reaction", "set_pool"}
PREFIXED = {"pools.csv": ("pool_id",), "pool_carriers.csv": ("pool_id",),
            "pool_transitions.csv": ("pool_id", "path_id"), "carrier_loops.csv": ("loop_id",)}
UNION = ("nodes.csv", "stid_to_uuid_mapping.csv", "entity_reaction_proxy_mapping.csv",
         "containment.csv", "cofactors.csv", "drugs.csv", "pathogens.csv", "carrier_loops.csv",
         "pools.csv", "pool_carriers.csv", "pool_transitions.csv", "boundary_edges.csv", "rnas.csv")
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
    both = [r for t in ("logic_network.csv", "boundary_edges.csv") for r in rows(d / t)[1]]
    edges = [(r["source_id"], r["target_id"], r.get("edge_reaction_id") or "") for r in both]
    etyped = [(r["source_id"], r["target_id"], r.get("edge_reaction_id") or "", r.get("edge_type") or "") for r in both]
    rnode = {r["uuid"]: r["diagram_entity_id"] for r in nodes if r["node_kind"] == "reaction"}
    ident = {}
    for r in nodes:
        if r["node_kind"] in NEVER or r["diagram_entity_id"] in held:
            continue
        ident[r["uuid"]] = tuple(r.get(k, "") for k in KEY)
    return {"ident": ident, "edges": edges, "etyped": etyped, "rnode": rnode, "reactions": set(rnode.values())}


def bench_roots(v, edges, nodes):
    """The members of `nodes` that the bench counts as roots over `edges` (a subset of pathway `v`'s edges): in-degree 0
    by bench.benchmark_vs_mpbiopath.root_indegree, so a depletion in-edge (under the default DS_PIN_ROOT_IGNORE) and a
    recycle edge do not stop a copy being a root. A curated bridge may enter only a copy the bench would pin (review
    2026-11-05: the builders had counted every in-edge)."""
    keep = set(edges)
    indeg = bm.root_indegree({"source_id": a, "target_id": b, "edge_type": t}
                             for a, b, rid, t in v["etyped"] if (a, b, rid) in keep)
    return {u for u in nodes if not indeg[u]}


def reduced_edges(q, without):
    """Q's edges with every reaction in `without` removed: its reaction nodes,
    their edges, and edges attributed to it."""
    return [(a, b, rid) for a, b, rid in q["edges"]
            if rid not in without and q["rnode"].get(a) not in without and q["rnode"].get(b) not in without]


def plan(base, hub):
    """Bridges for every ordered pair (A, B): {(A, B): [(src in A, dst in B)]},
    where the side that is not the centre pathway is reduced by the centre's
    reactions. Keyed (P, Q, direction) as {(P, Q): {"in": [...], "out": [...]}}."""
    paths = sorted(p for p in base.iterdir() if (p / "nodes.csv").exists())
    info = {p.name: load(p) for p in paths}
    n_paths = collections.Counter()
    for v in info.values():
        for de in {k[1] for k in v["ident"].values()}:
            n_paths[de] += 1
    by_key = collections.defaultdict(lambda: collections.defaultdict(list))
    for name, v in info.items():
        for u, k in v["ident"].items():
            if n_paths[k[1]] <= hub:
                by_key[k][name].append(u)
    shared = collections.defaultdict(list)        # (P, Q) -> [(uuids in P, uuids in Q)]
    for k, per in by_key.items():
        for a in per:
            for b in per:
                if a != b:
                    shared[(a, b)].append((per[a], per[b]))
    indeg = {n: collections.Counter(t for _, t, _ in v["edges"]) for n, v in info.items()}
    bridges = {}
    for (P, Q), pairs in shared.items():
        qe = reduced_edges(info[Q], info[P]["reactions"])
        qin = collections.Counter(t for _, t, _ in qe)
        into = [(s, t) for ps, qs in pairs for s in qs if qin[s] for t in ps if not indeg[P][t]]
        out = [(s, t) for ps, qs in pairs for s in ps if indeg[P][s] for t in qs if not qin[t]]
        if into or out:
            bridges[(P, Q)] = {"in": into, "out": out, "edges": qe}
    return [p.name for p in paths], bridges, n_paths


def plan_curated(base, links, s_q, centres, l2_only=False, members=None, whole=False):
    """Bridges from curated handoffs only: {(P, Q): {"in": [(src, dst)], "out": [], "edges": Q's reduced edges}},
    plus one census row per handoff pair and entity (why it did or did not make a bridge). `members` (a TSV of
    entity, member) lists each handed-off set's members: the generator expands a set into its members, so a set
    handoff is a handoff of each member, matched member to same member. Handoffs flagged `extracellular` (specs/133;
    a LINKS file from closure_133.py) are never bridged, as under plan_transitive."""
    mem = collections.defaultdict(set)
    if members:
        for r in csv.DictReader(open(members), delimiter="\t"):
            mem[r["entity"]].add(r["member"])
    q_of = {r["s"]: r["q_diagram"] for r in csv.DictReader(open(s_q), delimiter="\t")}
    links = [r for r in csv.DictReader(open(links), delimiter="\t") if r["pathway"] in centres
             and (r["l2"] == "True" or not l2_only)]
    info, ents, roots = {}, {}, {}

    def get(name):
        if name not in info:
            d = base / name
            v = load(d)
            v["typed"] = [(r["source_id"], r["target_id"], r["edge_type"], r.get("edge_reaction_id") or "")
                          for r in rows(d / "logic_network.csv")[1]]
            info[name] = v
            ents[name] = {r["uuid"]: r["diagram_entity_id"] for r in rows(d / "nodes.csv")[1]}
            roots[name] = bench_roots(v, v["edges"], ents[name])
        return info[name]

    bridges, census = {}, []
    for L in links:
        P, Q = L["pathway"], q_of.get(L["s"])
        row = dict(L, q=Q or "", status="")
        if L.get("extracellular") == "True":
            row["status"] = "extracellular"; census.append(row); continue
        if not Q or not (base / Q / "nodes.csv").exists():
            row["status"] = "no_source_network"; census.append(row); continue
        if Q == P:
            row["status"] = "source_is_centre"; census.append(row); continue
        vp, vq = get(P), get(Q)
        row["status"], pairs = match_link(L, base / P, vp, vq, ents[P], ents[Q], roots[P], mem)
        if pairs:
            if (P, Q) not in bridges:
                bridges[(P, Q)] = {"in": [], "out": [],
                                   "edges": vq["edges"] if whole else reduced_edges(vq, vp["reactions"])}
            bridges[(P, Q)]["in"] += pairs
        census.append(row)
    for b in bridges.values():
        b["in"] = sorted(set(b["in"]))
    return bridges, census


def match_link(L, xdir, vx, vq, ex, eq, roots_x, mem):
    """(census status, [(src in Q, dst in X)]) for one curated handoff (s in Q precedes r in X through entity e):
    Q's copies of e that s outputs, bridged to X's copies of e that feed r and are in `roots_x`. A set entity is matched
    member to same member: its members (from `mem`) count as e, and a set POOL node at r (specs/033: the set's own id,
    fed by `set_member` edges, as a set-valued catalyst or regulator is built) is replaced by the members feeding it.
    `no_same_member_pair` when copies exist on both sides and at a root but no two are the same entity."""
    e = L["entity"]
    held = set()
    for t in ("cofactors.csv", "drugs.csv", "pathogens.csv"):
        held |= {r["stable_id"] for r in rows(xdir / t)[1]}
    if e in held:
        return "held_entity", []
    ids = ({e} | mem[e]) - held
    role = "regulator" if L["role"].startswith("regulator") else L["role"]
    # an input edge carries its reaction's stId; a catalyst or regulator edge carries none and points at the
    # reaction node, whose diagram entity is the reaction
    dst = {a for a, b, t, rid in vx["typed"] if t == role and ex.get(a) in ids
           and (rid == L["r"] if role == "input" else ex.get(b) == L["r"])}
    feeds = collections.defaultdict(set)
    for a, b, t, _ in vx["typed"]:
        if t == "set_member" and b in dst and ex.get(a) not in held:
            feeds[b].add(a)
    for pool, ms in feeds.items():
        dst = (dst - {pool}) | ms
        ids |= {ex[m] for m in ms}
    dst = sorted(dst)
    src = sorted({b for a, b, t, rid in vq["typed"] if rid == L["s"] and t == "output" and eq.get(b) in ids})
    if not dst:
        return "entity_not_at_r_in_centre", []
    if not src:
        return "entity_not_output_of_s_in_source", []
    if not any(u in roots_x for u in dst):
        return "centre_copy_not_root", []
    pairs = [(a, b) for a in src for b in dst if b in roots_x and eq[a] == ex[b]]
    return ("bridged", pairs) if pairs else ("no_same_member_pair", [])


def plan_transitive(base, links, s_q, centres, members=None, whole=False, l2_only=False):
    """specs/133: curated handoffs followed back as far as they go, as a tree per centre P. Depth 0 is P; a pathway Q
    joins at depth d+1 when a handoff runs from Q into a root copy that a depth-d pathway X keeps (every root of P;
    for X != P a root among X's kept nodes, the ancestors of its own bridge sources). Each pathway joins once, at its
    least depth, and bridges run only from depth d+1 into depth d, so the context has no cycle. A joining pathway loses
    the reactions of every pathway at a lower depth. Handoffs flagged `extracellular` (cell-to-cell signals) are never
    bridged. Returns ({(P, Q): {"in": [(src in Q, dst)], "out": [], "edges": Q's reduced edges}}, census rows with
    the centre, the receiving pathway and Q's depth). Under `whole` (amendment 2) no pathway is reduced. Under `l2_only`
    a handoff counts only if its source lies in a pathway the receiving pathway's diagrams encapsulate (closure_133's
    `l2`, computed for the receiving pathway), as plan_curated's --l2."""
    mem = collections.defaultdict(set)
    if members:
        for r in csv.DictReader(open(members), delimiter="\t"):
            mem[r["entity"]].add(r["member"])
    q_of = {r["s"]: r["q_diagram"] for r in csv.DictReader(open(s_q), delimiter="\t")}
    by_x = collections.defaultdict(list)
    for r in csv.DictReader(open(links), delimiter="\t"):
        if r["l2"] == "True" or not l2_only:
            by_x[r["pathway"]].append(r)
    info, ents = {}, {}

    def get(name):
        if name not in info:
            d = base / name
            v = load(d)
            v["typed"] = [(r["source_id"], r["target_id"], r["edge_type"], r.get("edge_reaction_id") or "")
                          for r in rows(d / "logic_network.csv")[1]]
            info[name] = v
            ents[name] = {r["uuid"]: r["diagram_entity_id"] for r in rows(d / "nodes.csv")[1]}
        return info[name]

    bridges, census = {}, []
    for P in sorted(centres):
        vp = get(P)
        depth, edges = {P: 0}, {P: vp["edges"]}
        seen_rx = set(vp["reactions"])           # reactions of every pathway joined so far
        allowed = {P: bench_roots(vp, vp["edges"], ents[P])}
        level = [P]
        while level:
            joined = collections.defaultdict(list)       # Q -> its bridges this round
            for X in level:
                vx = get(X)
                for L in sorted(by_x.get(X, []), key=lambda r: (r["r"], r["s"], r["entity"])):
                    Q = q_of.get(L["s"])
                    row = dict(L, centre=P, q=Q or "", depth=depth[X] + 1, status="")
                    if L.get("extracellular") == "True":
                        row["status"] = "extracellular"
                    elif not Q or not (base / Q / "nodes.csv").exists():
                        row["status"] = "no_source_network"
                    elif Q in depth and depth[Q] <= depth[X]:
                        row["status"] = "source_already_joined"
                    else:
                        row["status"], pairs = match_link(L, base / X, vx, get(Q), ents[X], ents[Q], allowed[X], mem)
                        if pairs:
                            joined[Q] += pairs
                    census.append(row)
            nxt = []
            for Q in sorted(q for q, ps in joined.items() if ps):
                vq = get(Q)
                qe = vq["edges"] if whole else reduced_edges(vq, seen_rx)
                depth[Q], edges[Q] = depth[level[0]] + 1, qe
                srcs = sorted({a for a, _ in joined[Q]})
                kept = kept_nodes(base / Q, srcs, qe)
                allowed[Q] = bench_roots(vq, qe, kept)
                bridges[(P, Q)] = {"in": sorted(set(joined[Q])), "out": [], "edges": qe}
                nxt.append(Q)
            for Q in nxt:
                seen_rx |= get(Q)["reactions"]
            level = nxt
    return bridges, census


def kept_nodes(qd, srcs, qe):
    """The ancestors of `srcs` over edges `qe` of pathway dir `qd`, closed over
    its pools and carrier loops."""
    inc = collections.defaultdict(set)
    for a, b, _ in qe:
        inc[b].add(a)
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


def build_one(base, out, P, nbrs, bridges, prune=True, feedback=True):
    d = out / P
    d.mkdir(parents=True)
    members = [P] + nbrs
    kept, edges = {P: None}, {}
    for q in nbrs:
        bq = bridges[(P, q)]
        edges[q] = {(a, b, rid) for a, b, rid in bq["edges"]}
        # unpruned: the whole reduced neighbour, still without P's reactions
        kept[q] = kept_nodes(base / q, [s for s, _ in bq["in"]], bq["edges"]) if prune else \
            {x for a, b, _ in bq["edges"] for x in (a, b)} | {s for s, _ in bq["in"]} | {t for _, t in bq["out"]}
    for t in UNION:
        head, allr = None, []
        for q in members:
            h, rs = rows(base / q / t)
            if h is None:
                continue
            head = head or h
            for r in rs:
                if q != P and t in ("logic_network.csv", "boundary_edges.csv") and \
                        (r["source_id"], r["target_id"], r.get("edge_reaction_id") or "") not in edges[q]:
                    continue
                if kept[q] is not None and not keep_row(t, r, kept[q]):
                    continue
                for c in PREFIXED.get(t, ()):
                    r[c] = f"{q}:{r[c]}"
                allr.append(r)
        if head is None:
            continue
        if t in ("cofactors.csv", "drugs.csv", "pathogens.csv", "containment.csv", "rnas.csv"):
            allr = list({tuple(r[c] for c in head): r for r in allr}.values())
        with open(d / t, "w", newline="") as f:
            w = csv.DictWriter(f, head, lineterminator="\n")
            w.writeheader(); w.writerows(allr)
    with open(d / "logic_network.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(EDGE_HEAD)
        for q in members:
            for r in rows(base / q / "logic_network.csv")[1]:
                if q != P and ((r["source_id"], r["target_id"], r.get("edge_reaction_id") or "") not in edges[q]
                               or not keep_row("logic_network.csv", r, kept[q])):
                    continue
                w.writerow([r[c] for c in EDGE_HEAD])
        n = 0
        for q in nbrs:
            bq = bridges[(P, q)]
            outs = [(s_, t_) for s_, t_ in bq["out"] if t_ in kept[q]] if feedback else []
            for s_, t_ in bq["in"] + outs:
                w.writerow([s_, t_, "pos", "or", "context_bridge", "1", ""])
                n += 1
    with open(d / "context_own.csv", "w", newline="") as f:
        f.write("uuid\n" + "".join(r["uuid"] + "\n" for r in rows(base / P / "nodes.csv")[1]))
    # specs/133 amendment 1: each node's pathway, so the solver's identity rules stay inside one pathway
    with open(d / "context_origin.csv", "w", newline="") as f:
        f.write("uuid,pathway\n")
        for q in members:
            for r in rows(base / q / "nodes.csv")[1]:
                if kept[q] is None or keep_row("nodes.csv", r, kept[q]):
                    f.write(f"{r['uuid']},{q}\n")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("base"); ap.add_argument("out", nargs="?")
    ap.add_argument("--hub", type=int, default=5)
    ap.add_argument("--census", action="store_true")
    ap.add_argument("--no-prune", dest="prune", action="store_false")
    ap.add_argument("--feed-in", dest="feedback", action="store_false",
                    help="no bridge out of P: neighbours only feed P's roots, so no cross-pathway cycle")
    ap.add_argument("--only", action="append", default=[], help="build only these pathways (checks)")
    ap.add_argument("--curated", nargs=2, metavar=("LINKS", "S_Q"), help="specs/132: bridge curated handoffs only")
    ap.add_argument("--centres", help="specs/132: pathway list (bench/catalog_pathways.tsv format) to build")
    ap.add_argument("--members", help="specs/132: entity, member TSV (analysis/set_members_132.py)")
    ap.add_argument("--l2", action="store_true", help="specs/132: only handoffs from pathways P's diagrams encapsulate")
    ap.add_argument("--whole", action="store_true",
                    help="specs/133 amendment 2: with --curated, source pathways enter whole (not reduced, not pruned)")
    ap.add_argument("--transitive", action="store_true",
                    help="specs/133: with --curated, follow handoffs back as a tree (LINKS from closure_133.py)")
    a = ap.parse_args(argv)
    base = Path(a.base).resolve()
    if a.curated:
        return main_curated(a, base)
    names, bridges, n_paths = plan(base, a.hub)
    nbr = collections.defaultdict(set)
    for (x, y), bs in bridges.items():
        if bs["in"]:                     # only a bridge INTO P can move P's values
            nbr[x].add(y)
    size = {p: sum(1 for _ in open(base / p / "nodes.csv")) - 1 for p in names}
    print(f"pathways {len(names)}; hub cutoff > {a.hub} pathways; "
          f"hub entities {sum(1 for n in n_paths.values() if n > a.hub)} of {len(n_paths)}")
    print(f"neighbour pairs {len(bridges)}; bridges into the centre "
          f"{sum(len(b['in']) for b in bridges.values())}, out of it {sum(len(b['out']) for b in bridges.values())} (before pruning)")
    print("pathway\town_nodes\tneighbours\tcontext_nodes\tbridges_in_star")
    for p in names:
        ns = sorted(nbr[p])
        nb = sum(len(bridges[(p, q)]["in"]) + len(bridges[(p, q)]["out"]) for q in ns)
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
                  "bridges": build_one(base, out, p, sorted(nbr[p]), bridges, a.prune, a.feedback)}
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
              "context_pruned": a.prune, "context_feedback": a.feedback, "context": tot,
              "pathways_built": len(tot), "status": "complete" if not a.only else "partial"})
    json.dump(b, open(out / "BUILD.json", "w"), indent=1)


def main_curated(a, base):
    if not a.centres:
        sys.exit("--curated needs --centres")
    centres = [l.split("\t")[0] for l in open(a.centres) if l.startswith("R-HSA")]
    with open(a.curated[0]) as f:
        has_extra = "extracellular" in f.readline().rstrip("\n").split("\t")
    if a.transitive:
        bridges, census = plan_transitive(base, a.curated[0], a.curated[1], set(centres), a.members, a.whole, a.l2)
    else:
        bridges, census = plan_curated(base, a.curated[0], a.curated[1], set(centres), a.l2, a.members, a.whole)
    nbr = collections.defaultdict(set)
    for (p, q), b in bridges.items():
        nbr[p].add(q)
    print(f"centres {len(centres)}; handoff rows {len(census)}; "
          + "; ".join(f"{k} {v}" for k, v in sorted(collections.Counter(r['status'] for r in census).items())))
    print(f"centres with a bridge {len(nbr)}; bridges {sum(len(b['in']) for b in bridges.values())}")
    if a.census or not a.out:
        return
    out = Path(a.out).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    out.mkdir(parents=True)
    with open(out / "curated_census.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, list(census[0]) if census else ["status"], delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(census)
    tot = {}
    for p in centres:
        tot[p] = {"neighbours": sorted(nbr[p]),
                  "bridges": build_one(base, out, p, sorted(nbr[p]), bridges, not a.whole, False)}
    b = json.load(open(base / "BUILD.json"))
    repo = Path(__file__).resolve().parent.parent
    sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip() or "unknown"
    dirty = bool(subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", "bench/context_build.py"],
                                capture_output=True, text=True).stdout.strip())
    b.update({"build_id": out.name, "context_builder_commit": sha + ("-dirty" if dirty else ""),
              "context_of": base.name,
              "context_rule": ("curated_transitive" + ("_l2" if a.l2 else "")) if a.transitive else
                              "curated_l2" if a.l2 else "curated_l1",
              "context_l2": a.l2,
              # specs/133: extracellular handoffs are never bridged; a LINKS file without the column flags none
              "context_extracellular": "never_bridged" if has_extra else "no column",
              "context_root_ignore": bm.PIN_ROOT_IGNORE or "none",
              "context_links": a.curated[0], "context_members": a.members, "context_pruned": not a.whole,
              "context_whole": a.whole, "context_feedback": False, "context": tot,
              "pathway_list": a.centres, "pathways_requested": len(centres), "pathways_built": len(tot),
              "status": "complete"})
    json.dump(b, open(out / "BUILD.json", "w"), indent=1)


if __name__ == "__main__":
    main()
