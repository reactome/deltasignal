#!/usr/bin/env python3
"""Detect "siloed UUID" cases in a logic network.

When the logic-network generator positionally decomposes a single Reactome
stable_id into multiple UUIDs but does not bridge them, signal cannot flow
from upstream-reaction-output UUIDs to downstream-reaction-input UUIDs
even though they represent the SAME biological entity. The result is
gene→readout cascades that dead-end inside a single PE.

Two detection modes:

  STRICT silo (default count): a stable_id has at least one pure
  producer UUID (incoming, no outgoing), at least one pure consumer UUID
  (outgoing, no incoming), and ZERO bridge UUIDs. Guaranteed no flow.

  EFFECTIVE silo (with --effective): per-stid forward BFS from each
  producer-only UUID through the full network — flag the stid if at
  least one producer's output never reaches any of the stid's
  consumer-only UUIDs. Catches additional cases where bridge UUIDs
  exist but don't connect the producer half to the consumer half.

Known limitation
----------------
Neither mode catches "internally-fragmented all-bridge" stids — when
every UUID has both incoming AND outgoing edges, but the UUIDs fall
into disconnected sub-clusters in the broader network. Example:
R-HSA-3209194 (TP53 Tetramer) has 13 UUIDs, all bridges, but 9 of them
only feed a dissociation reaction while the signaling hub UUID has its
own dedicated upstream that the mainstream producers can't reach. This
pattern needs case-by-case investigation via `trace_paths.py`.

Report
------
  Per-pathway totals (multi-UUID stids, siloed stids, percent siloed)
  Top-N worst offenders by UUID count (with Reactome name lookup if py2neo
  is available)

Usage
-----
  python check_silo_bug.py                 # all 9 experimental pathways
  python check_silo_bug.py PATHWAY_DIR     # one specific pathway folder
  python check_silo_bug.py --top 30        # show more worst-offenders

Re-run after a generator change; the silo% should drop. The TP53 Tetramer
(R-HSA-3209194) is the headline bellwether — its silo state is a
single-number proxy for whether the fix worked.
"""
import csv, sys
from collections import defaultdict, deque
from pathlib import Path
from _common import CATALOG, EXP_PATHWAYS, pw_dir


def _bfs_reach(out_map, sources, sinks, max_depth=15):
    """Forward BFS in `out_map` (uuid → [target_uuids]) from any source until
    hitting any sink. Returns True if any sink is reachable."""
    seen = set(sources)
    q = deque((s, 0) for s in sources)
    while q:
        u, d = q.popleft()
        if u in sinks: return True
        if d >= max_depth: continue
        for t in out_map.get(u, ()):
            if t not in seen:
                seen.add(t); q.append((t, d + 1))
    return False


def _effective_silo(stid_uuids, inc, out, max_depth=15):
    """A stid is effectively siloed when at least one producer-only UUID
    cannot reach ANY of the stid's consumer-only UUIDs through the full
    network (within max_depth hops). Returns (is_silo, n_prod, n_cons, n_both)."""
    producers, consumers, both = [], [], []
    for u in stid_uuids:
        has_in, has_out = bool(inc.get(u)), bool(out.get(u))
        if has_in and has_out: both.append(u)
        elif has_in:           producers.append(u)
        elif has_out:          consumers.append(u)
    if not producers or not consumers:
        return False, len(producers), len(consumers), len(both)
    # Build forward target map from outgoing edges
    forward = {u: [t for t, *_ in out.get(u, [])] for u in out}
    # If ANY producer can't reach ANY consumer, the stid is effectively siloed.
    for p in producers:
        if not _bfs_reach(forward, {p}, set(consumers), max_depth):
            return True, len(producers), len(consumers), len(both)
    return False, len(producers), len(consumers), len(both)


def silo_analysis(dirname, effective=False):
    """Returns (totals, worst) for one pathway directory.
    totals = dict with multi, silo, prod_orphans, cons_orphans, bridges, all_uuids.
    worst  = list of (stid, n_uuids, n_producers, n_consumers, n_bridges) for
             siloed stids only, sorted by descending n_uuids.
    If effective=True, uses BFS to also catch stids whose bridges exist
    but don't connect producers to consumers."""
    inc = defaultdict(list)
    out = defaultdict(list)  # uuid → [(target, pn, ao, et)]
    fp = CATALOG / dirname / "logic_network.csv"
    for r in csv.reader(open(fp)):
        if r[0] == "source_id":
            continue
        src, tgt, pn, ao, et, _ = r
        inc[tgt].append(src)
        out[src].append((tgt, pn, ao, et))
    rev = defaultdict(list)
    for r in csv.DictReader(open(CATALOG / dirname / "stid_to_uuid_mapping.csv")):
        rev[r["stable_id"]].append(r["uuid"])

    multi = silo = 0
    worst = []
    all_uuids = set()
    prod_orphans = cons_orphans = bridges = orphans = 0
    for stid, uuids in rev.items():
        for u in uuids:
            all_uuids.add(u)
            has_in, has_out = bool(inc.get(u)), bool(out.get(u))
            if has_in and has_out: bridges += 1
            elif has_in:           prod_orphans += 1
            elif has_out:          cons_orphans += 1
            else:                  orphans += 1
        if len(uuids) < 2: continue
        multi += 1
        n_prod = n_cons = n_both = 0
        for u in uuids:
            has_in, has_out = bool(inc.get(u)), bool(out.get(u))
            if has_in and has_out: n_both += 1
            elif has_in:           n_prod += 1
            elif has_out:          n_cons += 1
        if effective:
            is_silo, _, _, _ = _effective_silo(uuids, inc, out)
        else:
            is_silo = bool(n_prod and n_cons and not n_both)
        if is_silo:
            silo += 1
            worst.append((stid, len(uuids), n_prod, n_cons, n_both))
    worst.sort(key=lambda x: -x[1])
    return ({
        "multi": multi, "silo": silo,
        "prod_orphans": prod_orphans, "cons_orphans": cons_orphans,
        "bridges": bridges, "orphans": orphans,
        "all_uuids": len(all_uuids),
    }, worst)


def name_lookup(stids):
    """Best-effort Reactome name lookup. Returns dict stid → displayName."""
    try:
        from py2neo import Graph
        g = Graph("bolt://localhost:7687", auth=("neo4j", "reactome"))
        Q = """UNWIND $stids AS s
        MATCH (p {stId: s})
        OPTIONAL MATCH (p)-[:referenceEntity]->(re)
        RETURN p.stId AS stid, p.displayName AS name,
               coalesce(re.geneName[0], '') AS gene"""
        return {r["stid"]: (r["name"], r["gene"])
                for r in g.run(Q, stids=stids).data()}
    except Exception:
        return {}


def main(target=None, top=10, effective=False):
    if target:
        dirs = [target]
    else:
        dirs = [pw_dir(p) for p in EXP_PATHWAYS]

    mode = "EFFECTIVE (BFS)" if effective else "STRICT (no-bridge)"
    print(f"Mode: {mode}")
    print(f"{'pathway':45s} {'#stids':>7s} {'#multi':>7s} "
          f"{'#siloed':>8s} {'silo%':>6s}   "
          f"{'prod-only':>10s} {'cons-only':>10s} {'bridges':>8s}")
    grand = defaultdict(int)
    all_worst = []
    for d in dirs:
        if not (CATALOG / d).exists():
            print(f"  {d[:43]:45s} (missing)"); continue
        totals, worst = silo_analysis(d, effective=effective)
        name = d.split("_R-HSA-")[0][:43]
        n_stids = totals["multi"] + (totals["bridges"] + totals["prod_orphans"]
                                     + totals["cons_orphans"] + totals["orphans"])
        # Note: n_stids above is approximate; multi already excludes singletons.
        # Reporting "multi" + raw counts is more honest.
        sp = totals["silo"] * 100 / totals["multi"] if totals["multi"] else 0
        print(f"  {name:45s} {len(worst) + (totals['multi'] - totals['silo']):>7d} "
              f"{totals['multi']:>7d} {totals['silo']:>8d} {sp:5.1f}%   "
              f"{totals['prod_orphans']:>10d} {totals['cons_orphans']:>10d} "
              f"{totals['bridges']:>8d}")
        for k in ("multi", "silo", "prod_orphans", "cons_orphans", "bridges"):
            grand[k] += totals[k]
        for stid, n, p, c, b in worst:
            all_worst.append((d, stid, n, p, c, b))

    sp = grand["silo"] * 100 / grand["multi"] if grand["multi"] else 0
    print(f"\n  {'TOTAL':45s} {' '*7} {grand['multi']:>7d} "
          f"{grand['silo']:>8d} {sp:5.1f}%   "
          f"{grand['prod_orphans']:>10d} {grand['cons_orphans']:>10d} "
          f"{grand['bridges']:>8d}")

    # Top worst offenders by UUID count
    all_worst.sort(key=lambda x: -x[2])
    show = all_worst[:int(top)]
    names = name_lookup([w[1] for w in show])
    print(f"\n=== Top {len(show)} siloed stids by UUID count ===")
    print(f"  {'pathway':30s} {'stid':16s} {'#uuids':>7s} {'prod':>5s} "
          f"{'cons':>5s}  name")
    for d, stid, n, p, c, _ in show:
        nm, gene = names.get(stid, ("?", ""))
        pw = d.split("_R-HSA-")[0][:28]
        gene_s = f"[{gene}]" if gene else ""
        print(f"  {pw:30s} {stid:16s} {n:7d} {p:5d} {c:5d}  "
              f"{nm[:55] if nm else '?':55s} {gene_s}")

    # Headline indicator — only meaningful under effective mode
    print(f"\n=== Bellwether: R-HSA-3209194 (TP53 Tetramer) ===")
    for d in dirs:
        if not (CATALOG / d).exists(): continue
        _, worst = silo_analysis(d, effective=True)
        for stid, n, p, c, b in worst:
            if stid == "R-HSA-3209194":
                pw = d.split("_R-HSA-")[0]
                print(f"  {pw}: SILOED ({n} UUIDs, {p} producers, "
                      f"{c} consumers, {b} bridges that don't connect halves)")
                break
        else:
            continue
        break
    else:
        print("  R-HSA-3209194 either not in 9-pathway scope or no longer siloed.")


if __name__ == "__main__":
    args = sys.argv[1:]
    top = 10
    effective = False
    if "--top" in args:
        i = args.index("--top")
        top = int(args[i + 1])
        args = args[:i] + args[i + 2:]
    if "--effective" in args:
        effective = True
        args.remove("--effective")
    target = args[0] if args else None
    main(target, top, effective)
