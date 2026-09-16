#!/usr/bin/env python3
"""What do the diagram's `links` connect that the generated network does not?

`diagram_connectivity.py` reads the layout JSON but uses only its `edges` —
the reaction glyphs — to bridge a producer reaction to a consumer reaction
across a shared glyph. It never reads `layout["links"]`, and those are a
separate population of connections drawn on the diagram:

    Interaction                  2,995
    EntitySetAndMemberLink       2,693
    EntitySetAndEntitySetLink      841
    FlowLine                        78

`EntitySetAndMemberLink` and `EntitySetAndEntitySetLink` largely restate
set membership, which Neo4j already has as hasMember/hasCandidate. `Interaction`
and `FlowLine` do not correspond to curated reactions at all, so nothing else
in the pipeline can recover them.

This reports, per pathway: how many links there are, how many join two entities
that BOTH appear as nodes in the generated network, and how many of those pairs
are not already connected there. The last number is the size of the gap.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
from pathlib import Path


def network_reachability(pathway_dir: Path):
    """stable_id -> set of stable_ids reachable from it in the logic network."""
    stid = {}
    f = pathway_dir / "stid_to_uuid_mapping.csv"
    if not f.exists():
        return None, None
    with f.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("stable_id"):
                stid[r["uuid"]] = r["stable_id"]
    adj = collections.defaultdict(set)
    lf = pathway_dir / "logic_network.csv"
    if not lf.exists():
        return None, None
    with lf.open() as fh:
        for e in csv.DictReader(fh):
            adj[e["source_id"]].add(e["target_id"])
    by_stid = collections.defaultdict(set)
    for uuid, s in stid.items():
        by_stid[s].add(uuid)
    return adj, by_stid


def reach(adj, sources, cap=120000):
    seen = set(sources)
    stack = list(sources)
    while stack:
        u = stack.pop()
        for v in adj.get(u, ()):
            if v not in seen:
                seen.add(v)
                stack.append(v)
                if len(seen) > cap:
                    return seen
    return seen


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--diagrams", type=Path,
                    default=Path.home() / "reactome-diagrams" / "97")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()

    totals = collections.Counter()
    per_pathway = []

    for d in sorted(p for p in args.catalog.iterdir() if p.is_dir()):
        if "_R-HSA-" not in d.name:
            continue
        pid = "R-HSA-" + d.name.rsplit("_R-HSA-", 1)[1]
        layout = args.diagrams / f"{pid}.json"
        if not layout.exists():
            totals["pathway_without_own_diagram"] += 1
            continue
        try:
            lay = json.loads(layout.read_text())
        except Exception:
            totals["layout_unreadable"] += 1
            continue

        glyph_to_entity = {n["id"]: n.get("reactomeId")
                           for n in lay.get("nodes", []) if n.get("reactomeId")}
        adj, by_stid = network_reachability(d)
        if adj is None:
            continue
        # entity dbId -> stable id is not in the layout, so match on dbId via
        # the companion graph.json node list.
        gfile = args.diagrams / f"{pid}.graph.json"
        dbid_to_stid = {}
        if gfile.exists():
            try:
                gr = json.loads(gfile.read_text())
                dbid_to_stid = {n["dbId"]: n["stId"]
                                for n in gr.get("nodes", []) if n.get("stId")}
            except Exception:
                pass

        present = joinable = unconnected = 0
        kinds = collections.Counter()
        cache = {}
        for link in lay.get("links", []) or []:
            cls = link.get("renderableClass") or "?"
            totals["link:" + cls] += 1
            kinds[cls] += 1
            src = [glyph_to_entity.get(x.get("id")) for x in (link.get("inputs") or [])]
            dst = [glyph_to_entity.get(x.get("id")) for x in (link.get("outputs") or [])]
            s_st = {dbid_to_stid.get(x) for x in src if x}
            t_st = {dbid_to_stid.get(x) for x in dst if x}
            s_st.discard(None)
            t_st.discard(None)
            if not s_st or not t_st:
                continue
            present += 1
            s_uu = {u for s in s_st for u in by_stid.get(s, ())}
            t_uu = {u for s in t_st for u in by_stid.get(s, ())}
            if not s_uu or not t_uu:
                continue
            joinable += 1
            key = tuple(sorted(s_uu))
            if key not in cache:
                cache[key] = reach(adj, s_uu)
            if not (t_uu & cache[key]):
                unconnected += 1
                totals["unconnected:" + cls] += 1

        totals["links_total"] += sum(kinds.values())
        totals["links_both_entities_known"] += present
        totals["links_both_in_network"] += joinable
        totals["links_not_already_connected"] += unconnected
        if unconnected:
            per_pathway.append((unconnected, joinable, d.name))

    print("diagram `links` across the catalog:")
    for k in sorted(totals):
        if k.startswith("link:"):
            print(f"   {k:<44} {totals[k]:>6}")
    print()
    for k in ("links_total", "links_both_entities_known", "links_both_in_network",
              "links_not_already_connected", "pathway_without_own_diagram"):
        print(f"   {k:<44} {totals[k]:>6}")
    print("\n   not already connected, by link class:")
    for k in sorted(totals):
        if k.startswith("unconnected:"):
            print(f"      {k[12:]:<41} {totals[k]:>6}")
    per_pathway.sort(reverse=True)
    print(f"\ntop pathways by unconnected diagram links:")
    for n, j, name in per_pathway[:args.top]:
        print(f"   {n:>5} of {j:>5} joinable   {name[:52]}")


if __name__ == "__main__":
    main()
