#!/usr/bin/env python3
"""specs/174 registered check: do K562's curated regulator -> target pairs fail because the knockdown moves nothing?

Pairs: the knocked-down gene's protein is in a regulator of a gene-expression reaction (a reaction with a gene or
transcript leaf, `ensembl:`) whose output carries the target's protein, in the served Release 98 networks. Symbols map
to UniProt by the reviewed human proteome FASTA (GN=), so no Neo4j is needed.

Usage: k562_174.py BUILD FASTA TRUTH OUT_JSON
"""
import csv
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def symbol_to_acc(fasta):
    out = {}
    for ln in open(fasta):
        if ln.startswith(">sp|"):
            m = re.search(r"GN=(\S+)", ln)
            if m:
                out.setdefault(m.group(1), ln.split("|")[1])
    return out


def regulator_pairs(build):
    """{(regulator accession, target accession)}: the regulator is a leaf of a regulator node of an expression reaction
    whose output node has the target as a leaf."""
    pairs = set()
    for pdir in sorted(p for p in Path(build).iterdir() if (p / "logic_network.csv").exists()):
        ref = {r["leaf_stid"]: r["reference_stid"] for r in csv.DictReader(open(pdir / "leaf_refs.csv"))}
        leaves = {r["uuid"]: [x for x in r["member_leaves"].split("|") if x] for r in csv.DictReader(open(pdir / "nodes.csv"))}
        ins, outs, regs = defaultdict(set), defaultdict(set), defaultdict(set)
        for e in csv.DictReader(open(pdir / "logic_network.csv")):
            if e["edge_type"] == "output":
                outs[e["source_id"]].add(e["target_id"])
            elif e["edge_type"] == "input":
                ins[e["target_id"]].add(e["source_id"])
            elif e["edge_type"] == "regulator":
                regs[e["target_id"]].add(e["source_id"])

        def accs(u, prefix):
            return {ref[x].split(":", 1)[1] for x in leaves.get(u, ()) if ref.get(x, "").startswith(prefix)}
        for rx, rs in regs.items():
            if not any(accs(u, "ensembl:") for u in ins[rx] | outs[rx]):
                continue
            targets = set().union(*(accs(u, "uniprot:") for u in outs[rx])) if outs[rx] else set()
            for u in rs:
                for a in accs(u, "uniprot:"):
                    pairs.update((a, t) for t in targets if t != a)
    return pairs


def main():
    import dataset_roles
    dataset_roles.require_role("replogle2022_a", "score")
    build, fasta, truth, out = sys.argv[1:5]
    s2a = symbol_to_acc(fasta)
    pairs = regulator_pairs(build)
    rows = list(csv.DictReader(open(truth), delimiter="\t"))
    effect = defaultdict(list)
    for r in rows:
        effect[r["gene"]].append(r["measured_50"] != "not")
    hit = [r for r in rows if (s2a.get(r["gene"]), s2a.get(r["target"])) in pairs]
    changed = [r for r in hit if r["measured_50"] != "not"]
    failed = [r for r in hit if r["measured_50"] == "not"]

    def gene_effect(rs):
        return [sum(effect[r["gene"]]) / len(effect[r["gene"]]) for r in rs]

    fe = gene_effect(failed)
    res = {
        "regulator_pairs_in_build": len(pairs), "scored_pairs": len(hit), "changed": len(changed),
        "genes": len({r["gene"] for r in hit}),
        "failed_gene_moves_nothing": sum(x == 0 for x in fe),
        "failed_gene_effect_median": statistics.median(fe) if fe else None,
        "changed_gene_effect_median": statistics.median(gene_effect(changed)) if changed else None,
        "all_genes_effect_median": statistics.median(sum(v) / len(v) for v in effect.values()),
        "failed_knockdown_median": statistics.median(float(r["knockdown"]) for r in failed) if failed else None,
        "changed_knockdown_median": statistics.median(float(r["knockdown"]) for r in changed) if changed else None,
        "by_gene": {g: {"pairs": sum(r["gene"] == g for r in hit), "changed": sum(r["gene"] == g for r in changed),
                        "effect": round(sum(effect[g]) / len(effect[g]), 3), "targets_scored": len(effect[g])}
                    for g in sorted({r["gene"] for r in hit})},
    }
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "by_gene"}, indent=1))


if __name__ == "__main__":
    main()
