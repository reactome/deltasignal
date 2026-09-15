#!/usr/bin/env python3
"""Count stable ids split into a pure SINK and a pure SOURCE uuid.

The GPVI regression exposed a very specific shape of the UUID silo. The
catalyst `R-HSA-442307` occurs as two uuids that share a stable id and are
disconnected from each other:

    0f995db5…  in-degree 0, out-degree 5   feeds the reaction; nothing feeds it
    22fd0285…  in-degree 1, out-degree 0   receives SYK:p-VAV; goes nowhere

The signal arrives at one copy of the entity and the reaction reads the other,
so the perturbation has no route at all. This is narrower than "the silo bug":
it is specifically a sink/source pair for one curated entity, where a single
edge from the sink to the source would restore exactly the route the curation
implies and nothing else.

Three broader fixes have already failed — full merge (-204pp), all-pairs
bridges (2.3M edges), connectivity-aware one-bridge-per-silo (-4pp). So measure
the prevalence of THIS shape, AND what closing it would connect, before
building anything.

MEASURED 2026-09-15 over the 92-pathway catalog at Release97:

  4,415 stable ids carry a disconnected sink/source pair, in 91 of 92 pathways.
  Bridging all of them makes 160,090 nodes newly reachable from signals that
  currently die at a sink (that figure adds every bridge at once and walks from
  all sinks, so overlapping reach is counted once; the per-bridge figures below
  are each measured against the unbridged graph and do not sum to it).

That is not a targeted repair, and it explains why one-bridge-per-silo lost:
the marginal downstream reach of a single bridge is wildly skewed —

  min 2   p25 24   median 88   p75 452   p90 2,779   max 39,374

In `Transcriptional_regulation_by_RUNX1` eight bridges each connect roughly
39,300 of the network's 39,773 nodes. That is the hub-flooding signature that also cost
macro-F1 0.663 -> 0.479 in the cross-pathway stitch. But 38.9% of bridges reach
50 nodes or fewer and 52.8% reach 100 or fewer, and the GPVI catalyst that
surfaced all of this (`R-HSA-442307`) reaches 36 — a genuinely local repair.

So the untested hypothesis is a REACH-CAPPED bridge: restore the route only
where doing so does not connect an entity's arrival point to the whole network.
That is distinct from all three failed attempts, which bridged unconditionally.
Testing it needs a catalog regeneration plus a wide-set A/B per arm.
"""

from __future__ import annotations

import argparse
import csv
import collections
from pathlib import Path


def scan(pathway_dir: Path):
    stid_file = pathway_dir / "stid_to_uuid_mapping.csv"
    logic_file = pathway_dir / "logic_network.csv"
    if not stid_file.exists() or not logic_file.exists():
        return None

    stid = {}
    with stid_file.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("stable_id"):
                stid[r["uuid"]] = r["stable_id"]

    indeg: collections.Counter = collections.Counter()
    outdeg: collections.Counter = collections.Counter()
    with logic_file.open() as fh:
        for e in csv.DictReader(fh):
            outdeg[e["source_id"]] += 1
            indeg[e["target_id"]] += 1

    by_stid = collections.defaultdict(list)
    for uuid, s in stid.items():
        by_stid[s].append(uuid)

    pairs = []
    for s, uuids in by_stid.items():
        if len(uuids) < 2:
            continue
        sinks = [u for u in uuids if indeg[u] > 0 and outdeg[u] == 0]
        sources = [u for u in uuids if indeg[u] == 0 and outdeg[u] > 0]
        if sinks and sources:
            pairs.append((s, len(sinks), len(sources),
                          sum(outdeg[u] for u in sources)))
    return {
        "nodes": len(stid),
        "edges": sum(outdeg.values()),
        "split_stids": sum(1 for v in by_stid.values() if len(v) > 1),
        "sink_source_stids": len(pairs),
        "pairs": pairs,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--reach", action="store_true",
                    help="also report the marginal downstream reach of each bridge")
    args = ap.parse_args()

    rows = []
    total_pairs = 0
    for d in sorted(p for p in args.catalog.iterdir() if p.is_dir()):
        got = scan(d)
        if not got:
            continue
        rows.append((d.name, got))
        total_pairs += got["sink_source_stids"]

    print(f"pathways scanned: {len(rows)}")
    print(f"stable ids with a disconnected sink/source pair: {total_pairs}\n")
    rows.sort(key=lambda r: -r[1]["sink_source_stids"])
    print(f"{'pathway':<58}{'nodes':>7}{'split':>7}{'sink/src':>9}")
    for name, got in rows[:args.top]:
        print(f"{name[:56]:<58}{got['nodes']:>7}{got['split_stids']:>7}"
              f"{got['sink_source_stids']:>9}")
    affected = sum(1 for _, g in rows if g["sink_source_stids"])
    print(f"\npathways with at least one: {affected} of {len(rows)}")
    if args.reach:
        report_reach(rows, args.catalog)


def report_reach(rows, catalog: Path) -> None:
    """Marginal downstream reach of each bridge, taken one at a time.

    The aggregate is misleading: a handful of bridges connect an entire
    network while most are local. This is the distribution that decides
    whether a reach cap is a principled rule or just a tuned threshold.
    """
    import statistics

    gains = []
    for name, _ in rows:
        d = catalog / name
        stid = {}
        with (d / "stid_to_uuid_mapping.csv").open() as fh:
            for r in csv.DictReader(fh):
                if r.get("stable_id"):
                    stid[r["uuid"]] = r["stable_id"]
        adj = collections.defaultdict(list)
        indeg: collections.Counter = collections.Counter()
        outdeg: collections.Counter = collections.Counter()
        with (d / "logic_network.csv").open() as fh:
            for e in csv.DictReader(fh):
                adj[e["source_id"]].append(e["target_id"])
                outdeg[e["source_id"]] += 1
                indeg[e["target_id"]] += 1
        by = collections.defaultdict(list)
        for u, s in stid.items():
            by[s].append(u)
        for s, us in by.items():
            if len(us) < 2:
                continue
            sinks = [u for u in us if indeg[u] > 0 and outdeg[u] == 0]
            sources = [u for u in us if indeg[u] == 0 and outdeg[u] > 0]
            if not (sinks and sources):
                continue
            # TRUE marginal gain: what the SINK can newly reach once bridged,
            # not everything downstream of one source. The first version
            # measured the latter and overstated the local-repair fraction by
            # about a point (40.1% -> 38.9% at <=50 nodes); it also ignored all
            # but one source per stable id.
            sink = max(sinks, key=lambda u: indeg[u])

            def walk(start):
                seen, q = {start}, [start]
                while q:
                    u = q.pop()
                    for v in adj.get(u, ()):
                        if v not in seen:
                            seen.add(v)
                            q.append(v)
                return seen

            before = walk(sink)
            adj[sink].extend(sources)
            after = walk(sink)
            del adj[sink][-len(sources):]
            gains.append(len(after - before))

    if not gains:
        return
    print(f"\nmarginal downstream reach per bridge ({len(gains)} bridges):")
    print(f"  min {min(gains)}  p25 {int(statistics.quantiles(gains, n=4)[0])}"
          f"  median {int(statistics.median(gains))}"
          f"  p75 {int(statistics.quantiles(gains, n=4)[2])}"
          f"  p90 {int(statistics.quantiles(gains, n=10)[8])}  max {max(gains)}")
    for thr in (10, 50, 100, 500):
        n = sum(1 for v in gains if v <= thr)
        print(f"  reach <= {thr:>4}: {n:>5} ({100 * n / len(gains):.1f}%)")


if __name__ == "__main__":
    main()
