#!/usr/bin/env python3
"""specs/147 P2: every curated reaction and participant vs the generated bundle.

For each benchmark pathway, every ReactionLikeEvent Neo4j reaches by hasEvent*
(fetch_reaction_rows — the oracle's own fetch, FR: imported, never re-queried
differently) is checked against the bundle:

  reaction_absent        no reaction node carries its stable id
                         (subclassed by node_exclusions.csv: _excluded means
                          the generator RECORDED dropping it, with its reason)
  input_dropped          a curated input with no input edge into any copy
  catalyst_dropped       a curated catalyst with no catalyst edge into any copy
  regulator_dropped      a curated regulator with no regulator edge into any copy
  regulator_sign_flipped a Negative* regulation wired pos, or Positive* wired neg
  output_dropped         a curated output no copy produces

A participant counts as present when ANY node whose stable id matches it
directly (nodes.csv diagram_entity_id), or that CARRIES it (member_leaves),
or that carries any of its Neo4j composition leaves (hasMember, hasCandidate,
hasComponent — the generator's recorded set-expansion and decomposition
policies mean a set's or decomposed complex's own id never survives), or any
set_pool fed by such a node, has the right edge type into/out of any copy of
the reaction. So a finding is a participant none of whose forms, members or
leaves touches the reaction at all. Cofactor policy does not drop edges (DS_COFACTOR_MODE
pins them inert), so a dropped cofactor input still counts as a finding.

Writes per-class counts and a traced-examples TSV.

Usage: structural_diff.py CATALOG [--out DIR] [--pathways N]
"""

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / "bench" / "analysis"))
from curator_oracle import fetch_composition, fetch_reaction_rows  # noqa: E402


def load_bundle(pdir: Path):
    nodes = list(csv.DictReader(open(pdir / "nodes.csv")))
    edges = list(csv.DictReader(open(pdir / "logic_network.csv")))
    carriers = collections.defaultdict(set)   # stid -> uuids that are or carry it
    reaction_uuids = collections.defaultdict(list)
    for n in nodes:
        if n["node_kind"] == "reaction":
            reaction_uuids[n["diagram_entity_id"]].append(n["uuid"])
            continue
        if n["diagram_entity_id"]:
            carriers[n["diagram_entity_id"]].add(n["uuid"])
        for leaf in (n["member_leaves"] or "").split("|"):
            if leaf:
                carriers[leaf].add(n["uuid"])
    # a set_pool is fed by set_member edges from carriers; it stands in for them
    pool_feeds = collections.defaultdict(set)  # pool uuid -> source uuids
    for e in edges:
        if e["edge_type"] == "set_member":
            pool_feeds[e["target_id"]].add(e["source_id"])
    into = collections.defaultdict(list)   # reaction uuid -> (source, type, sign)
    outof = collections.defaultdict(list)  # reaction uuid -> target
    for e in edges:
        into[e["target_id"]].append((e["source_id"], e["edge_type"], e["pos_neg"]))
        outof[e["source_id"]].append((e["target_id"], e["edge_type"]))
    exclusions = {}
    xp = pdir / "node_exclusions.csv"
    if xp.exists():
        exclusions = {r["stable_id"]: r["reason"] for r in csv.DictReader(open(xp))}
    return carriers, reaction_uuids, pool_feeds, into, outof, exclusions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catalog", type=Path)
    ap.add_argument("--out", type=Path, default=HERE.parent)
    ap.add_argument("--pathways", type=int, default=0, help="first N only (0 = all)")
    args = ap.parse_args()

    pdirs = sorted(p for p in args.catalog.iterdir()
                   if p.name.startswith("R-HSA-") and (p / "nodes.csv").exists())
    if args.pathways:
        pdirs = pdirs[:args.pathways]

    counts = collections.Counter()
    per_pathway = collections.defaultdict(collections.Counter)
    examples = []
    totals = collections.Counter()

    for pdir in pdirs:
        pw = pdir.name
        rows = fetch_reaction_rows(pw)
        reactions = collections.defaultdict(lambda: {"inputs": set(), "outputs": set(),
                                                     "catalysts": set(), "regulators": {}})
        for stid, _cls, i, o, c, g, regcls in rows:
            r = reactions[stid]
            if i: r["inputs"].add(i)
            if o: r["outputs"].add(o)
            if c: r["catalysts"].add(c)
            if g: r["regulators"][g] = regcls
        carriers, reaction_uuids, pool_feeds, into, outof, exclusions = load_bundle(pdir)
        leaves_of, _comp_pairs, _set_pairs = fetch_composition(pw)

        def carrier_set(stid):
            cset = set(carriers.get(stid, ()))
            for leaf in leaves_of.get(stid, ()):
                cset |= carriers.get(leaf, set())
            return cset

        def present(stid, rx_uuids, want_types, want_sign=None):
            """Any copy of the reaction has an edge of a wanted type from a
            node that is/carries stid, or from a pool fed by one."""
            cset = carrier_set(stid)
            for ru in rx_uuids:
                for src, etype, sign in into[ru]:
                    if etype not in want_types:
                        continue
                    if want_sign is not None and sign != want_sign:
                        continue
                    if src in cset or (pool_feeds.get(src, set()) & cset):
                        return True
            return False

        for stid, r in sorted(reactions.items()):
            totals["reactions"] += 1
            totals["participants"] += (len(r["inputs"]) + len(r["outputs"])
                                       + len(r["catalysts"]) + len(r["regulators"]))
            rx = reaction_uuids.get(stid)
            if not rx:
                cls = "reaction_absent_excluded" if stid in exclusions else "reaction_absent"
                counts[cls] += 1
                per_pathway[pw][cls] += 1
                if len([e for e in examples if e["class"] == cls]) < 25:
                    examples.append({"class": cls, "pathway": pw, "reaction": stid,
                                     "detail": exclusions.get(stid, "")})
                continue
            for i in r["inputs"]:
                if not present(i, rx, {"input"}):
                    counts["input_dropped"] += 1
                    per_pathway[pw]["input_dropped"] += 1
                    if len([e for e in examples if e["class"] == "input_dropped"]) < 25:
                        examples.append({"class": "input_dropped", "pathway": pw,
                                         "reaction": stid, "detail": i})
            for c in r["catalysts"]:
                if not present(c, rx, {"catalyst"}):
                    counts["catalyst_dropped"] += 1
                    per_pathway[pw]["catalyst_dropped"] += 1
                    if len([e for e in examples if e["class"] == "catalyst_dropped"]) < 25:
                        examples.append({"class": "catalyst_dropped", "pathway": pw,
                                         "reaction": stid, "detail": c})
            for g, regcls in r["regulators"].items():
                negative = "Negative" in regcls
                ok = present(g, rx, {"regulator"}, "neg" if negative else "pos")
                if not ok and present(g, rx, {"regulator"}):
                    counts["regulator_sign_flipped"] += 1
                    per_pathway[pw]["regulator_sign_flipped"] += 1
                    examples.append({"class": "regulator_sign_flipped", "pathway": pw,
                                     "reaction": stid, "detail": f"{g} ({regcls})"})
                elif not ok:
                    counts["regulator_dropped"] += 1
                    per_pathway[pw]["regulator_dropped"] += 1
                    if len([e for e in examples if e["class"] == "regulator_dropped"]) < 25:
                        examples.append({"class": "regulator_dropped", "pathway": pw,
                                         "reaction": stid, "detail": f"{g} ({regcls})"})
            for o in r["outputs"]:
                cset = carrier_set(o)
                made = any(tgt in cset or (pool_feeds.get(tgt, set()) & cset)
                           for ru in rx for tgt, etype in outof[ru] if etype == "output")
                if not made:
                    counts["output_dropped"] += 1
                    per_pathway[pw]["output_dropped"] += 1
                    if len([e for e in examples if e["class"] == "output_dropped"]) < 25:
                        examples.append({"class": "output_dropped", "pathway": pw,
                                         "reaction": stid, "detail": o})
        print(f"{pw}: {len(reactions)} curated reactions, "
              f"{sum(per_pathway[pw].values())} findings", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / "examples.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, ["class", "pathway", "reaction", "detail"], delimiter="\t")
        w.writeheader()
        w.writerows(examples)
    summary = {"pathways": len(pdirs), "curated": dict(totals), "findings": dict(counts),
               "worst_pathways": {cls: sorted(((pp[cls], pw) for pw, pp in per_pathway.items()
                                               if pp[cls]), reverse=True)[:5]
                                  for cls in counts}}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["findings"], indent=1))
    print("curated:", dict(totals))


if __name__ == "__main__":
    main()
