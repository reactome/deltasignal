#!/usr/bin/env python3
"""How much of the decryptM reachability gap is the pathway boundary?

For each scored decryptM case, is the measured site reachable from a node built
from one of the drug's targets?
- inside the case's own pathway, and
- in the STITCHED graph: all pathways of the build merged by node identity,
  the key the 2026-07-16 cross-pathway measurement used (node_kind,
  diagram_entity_id, compartment, member_leaves, chosen_members).

Both use the same generous definition (any target-containing node, any edge).
The identity key also merges copies WITHIN a pathway (variant and silo twins),
so a third count, the keyed graph of the case's own pathway alone, separates
that from crossing the boundary. Graph only; nothing is solved.

Usage: decryptm_stitched_reach.py BUILD CASES_TSV TARGET_EWAS_CSV PHOSPHO_EWAS_CSV
  CASES_TSV is a decryptm_bench.py case log (status, split, measured, targets,
  uniprot, position, pathway_id, predicted, reach_predicted).
"""
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from network_files import open_network  # noqa: E402

RESIDUE = {"O-phospho-L-serine", "O-phospho-L-threonine", "O4'-phospho-L-tyrosine"}


def neo4j_rows(path):
    for line in open(path):
        if line.startswith('"'):
            yield [x.strip().strip('"') for x in line.split(",")]


def load_pathway(pdir):
    """(uuid -> identity key, uuid -> leaves, edges as (src, dst) uuids)."""
    key, leaves = {}, {}
    with open(pdir / "nodes.csv", newline="") as f:
        for r in csv.DictReader(f):
            key[r["uuid"]] = (r["node_kind"], r["diagram_entity_id"], r.get("compartment", ""),
                              r.get("member_leaves", ""), r.get("chosen_members", ""))
            leaves[r["uuid"]] = set(filter(None, (r.get("member_leaves") or "").split("|")))
            leaves[r["uuid"]].add(r["diagram_entity_id"])
    with open_network(pdir) as f:
        edges = [(r["source_id"], r["target_id"]) for r in csv.DictReader(f)]
    return key, leaves, edges


def reachable(adj, sources):
    seen, stack = set(sources), list(sources)
    while stack:
        for v in adj.get(stack.pop(), ()):
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def main(argv):
    build, cases_tsv, target_csv, phospho_csv = map(Path, argv[1:5])
    gene_ewas = collections.defaultdict(set)
    for st, g in neo4j_rows(target_csv):
        gene_ewas[g].add(st)
    site_ewas = collections.defaultdict(set)
    for st, acc, _g, coord, mod in neo4j_rows(phospho_csv):
        if mod in RESIDUE:
            site_ewas[(acc, int(coord))].add(st)
    cases = [r for r in csv.DictReader(open(cases_tsv), delimiter="\t") if r["status"] == "scored"]
    pathways = {}
    s_adj = collections.defaultdict(set)
    s_leaves = collections.defaultdict(set)
    for pdir in sorted(p for p in build.iterdir() if (p / "nodes.csv").exists()):
        key, leaves, edges = load_pathway(pdir)
        pathways[pdir.name] = (key, leaves, edges)
        for u, k in key.items():
            s_leaves[k] |= leaves[u]
        for a, b in edges:
            if a in key and b in key:
                s_adj[key[a]].add(key[b])
    s_has = collections.defaultdict(set)          # leaf -> stitched keys built from it
    for k, ls in s_leaves.items():
        for l in ls:
            s_has[l].add(k)
    s_reach_cache = {}
    tally = collections.Counter()
    keyed_tally = collections.Counter()
    for c in cases:
        key, leaves, edges = pathways[c["pathway_id"]]
        t_ewas = set().union(*(gene_ewas.get(g, set()) for g in c["targets"].split(";")))
        r_ewas = site_ewas.get((c["uniprot"], int(c["position"])), set())
        adj = collections.defaultdict(set)
        for a, b in edges:
            adj[a].add(b)
        src = {u for u, ls in leaves.items() if ls & t_ewas}
        dst = {u for u, ls in leaves.items() if ls & r_ewas}
        own = bool(reachable(adj, src) & dst) if src and dst else False
        tk = frozenset(t_ewas)
        if tk not in s_reach_cache:
            s_reach_cache[tk] = reachable(s_adj, set().union(*(s_has[e] for e in t_ewas)) if t_ewas else set())
        s_dst = set().union(*(s_has[e] for e in r_ewas)) if r_ewas else set()
        stitched = bool(s_reach_cache[tk] & s_dst)
        k_adj = collections.defaultdict(set)          # keyed, own pathway only
        for a, b in edges:
            if a in key and b in key:
                k_adj[key[a]].add(key[b])
        k_src = {key[u] for u in src}
        keyed = bool(reachable(k_adj, k_src) & {key[u] for u in dst}) if k_src and dst else False
        reg = "regulated" if c["measured"] in ("down", "up") else "unchanged"
        tally[(c["split"], reg, own, stitched)] += 1
        keyed_tally[(c["split"], reg)] += keyed
    for split in ("held-out", "tuning"):
        for reg in ("regulated", "unchanged"):
            n = sum(v for k, v in tally.items() if k[:2] == (split, reg))
            if not n:
                continue
            own = sum(v for k, v in tally.items() if k[:2] == (split, reg) and k[2])
            st = sum(v for k, v in tally.items() if k[:2] == (split, reg) and k[3])
            gain = sum(v for k, v in tally.items() if k[:2] == (split, reg) and k[3] and not k[2])
            print(f"{split:9s} {reg:10s} n={n:5d}  reachable in pathway {own:5d}  "
                  f"keyed own pathway {keyed_tally[(split, reg)]:5d}  stitched {st:5d}  "
                  f"(only when stitched {gain})")
    return tally


if __name__ == "__main__":
    main(sys.argv)
