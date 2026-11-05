#!/usr/bin/env python3
"""specs/135: the total change calls left under DS_TOTAL_READOUT=abundance, by how their abundance copies are made, and
what a basal expression level and a 24 h protein stock would do to them.

A copy is `starting` (no in-edge but depletion), `made` (by the reactions in the network's `output` edges into it) or
`other` (an in-edge from something else, e.g. assembled from partners only).
The producing reactions of made copies are printed with their counts, for a Neo4j check of what they are.

What-if (single-copy readouts only, where the readout is the copy): a basal level beside the expression regulators
reads v -> (1 + v) / 2; a stock of which a fraction f is replaced over the treatment reads v -> (1 - f) + f v, with
f = 1 - 2^(-T / t_half). Calls use the bench's cutoffs.

Usage: census_135.py KORKUT_CASES [T_HOURS T_HALF_HOURS]   (needs NEO4J_PASSWORD and DS_CATALOG_ROOT)
"""
import collections
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# the census is of abundance copies: refuse another readout rather than silently classify its copies
if not os.environ.setdefault("DS_TOTAL_READOUT", "abundance").startswith("abundance"):
    raise SystemExit(f"census_135 counts abundance copies; DS_TOTAL_READOUT={os.environ['DS_TOTAL_READOUT']!r}")
import benchmark_vs_mpbiopath as bm  # noqa: E402
import korkut_bench as kb  # noqa: E402


def copy_producers(copies, network_rows):
    """{copy: (kind, sorted producing reaction ids)}, kind `made`, `other` or `starting`."""
    prod, other = collections.defaultdict(set), set()
    for e in network_rows:
        if e["edge_type"] == "output":
            prod[e["target_id"]].add(e["source_id"])
        elif e["edge_type"] != "depletion":
            other.add(e["target_id"])
    return {u: ("made" if prod[u] else "other" if u in other else "starting", sorted(prod[u])) for u in copies}


def stock_fraction(t_hours, t_half_hours):
    """The share of a protein replaced over t_hours at a turnover half-life of t_half_hours."""
    return 1.0 - 2.0 ** (-t_hours / t_half_hours)


def what_if(v, basal=False, f=1.0):
    """A synthesis copy's fold under a basal expression level and a stock of which f is new."""
    if basal:
        v = (1.0 + v) / 2.0
    return (1.0 - f) + f * v


def call(ui):
    return {bm.DOWN: "down", bm.NORMAL: "not", bm.UP: "up"}[bm.classify(ui)]


def main():
    cases_path = sys.argv[1]
    t, th = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (24.0, 46.0)
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] == "total"]
    genes = sorted({g for r in readouts for g in r["sites"].split(",")})
    r_st = kb.readout_stids(readouts, kb.total_gene_stids(readouts, bm.neo4j_gene_to_stids(genes)), {})
    cases = [r for r in csv.DictReader(open(cases_path), delimiter="\t")
             if r["status"] == "scored" and r["readout"] in r_st]
    changes = collections.Counter((c["pathway_id"], c["readout"]) for c in cases if c["predicted"] != "not")
    cls, rx = {}, collections.Counter()
    for pid in sorted({c["pathway_id"] for c in cases}):
        pdir = bm.CATALOG_ROOT / pid
        with bm.open_network(pdir) as f:
            rows = list(csv.DictReader(f))
        diagram = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(open(pdir / "nodes.csv"))}
        for name, nodes in kb.readout_nodes(r_st, bm.load_stid_to_uuids(pdir)).items():
            prods = copy_producers(kb.abundance_copies(nodes, rows, kb.TOTAL_READOUT), rows)
            kinds = {k for k, _ in prods.values()}
            cls[(pid, name)] = "+".join(sorted(kinds)) or "none"
            for p in {diagram.get(r, r) for _, ps in prods.values() for r in ps}:
                rx[p] += changes[(pid, name)]
    out = collections.Counter()
    for c in cases:
        k = cls.get((c["pathway_id"], c["readout"]), "?")
        out[(k, "cases")] += 1
        if c["predicted"] != "not":
            out[(k, "change")] += 1
            out[(k, "right")] += c["predicted"] == c["measured"]
    print("copies           cases  change calls  right")
    for k in sorted({k[0] for k in out}):
        print(f"{k:15s} {out[(k, 'cases')]:6d} {out[(k, 'change')]:13d} {out[(k, 'right')]:6d}")
    print("\nproducing reactions of made copies (weighted by change calls):")
    for k, v in rx.most_common():
        if v:
            print(f"  {k} {v}")
    single = [c for c in cases if c["predicted"] != "not" and c["readout_nodes"] == "1"
              and cls.get((c["pathway_id"], c["readout"])) == "made"]
    f = stock_fraction(t, th)
    print(f"\nsingle-copy made readouts with a change call: {len(single)}; stock f = {f:.3f} (T {t:g} h, t1/2 {th:g} h)")
    for label, basal, ff in (("basal", True, 1.0), ("stock", False, f), ("basal + stock", True, f)):
        moves = collections.Counter()
        near = 0
        for c in single:
            u = what_if(float(c["pred_ui"]), basal, ff)
            moves[f"{c['predicted']}->{call(u)}"] += 1
            near += abs(u - bm.DOWN_CUTOFF) < 0.01
        print(f"  {label:14s} {dict(sorted(moves.items()))}  within 0.01 of the DOWN cutoff: {near}")


if __name__ == "__main__":
    main()
