#!/usr/bin/env python3
"""specs/182: does a Lun protein's starting form act by itself in Reactome? Structure only; no Lun value is read.

Usage: classes_182.py BUILD > activation_classes.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
sys.path.insert(0, str(HERE.parents[1] / "174"))
import lun_a_173 as A  # noqa: E402
import trace_174 as T  # noqa: E402

IDENTITY = {"variant_pool", "variant_split", "specialization", "set_member"}   # the same material under another node
ACTS = {"catalyst", "regulator"}
ORDER = "ABCN"


def classify(net, acc):
    """(class, starting forms, reactions entered) for one protein in one network."""
    want = f"uniprot:{acc}"

    def has(u):
        return any(net.ref.get(leaf) == want for leaf in net.leaves.get(u, ()))
    starts = [u for u in net.nodes_with(acc) if net.indeg[u] == 0]
    seen, stack = set(starts), list(starts)
    while stack:                                   # a starting form and the pool / variant nodes that stand for it
        u = stack.pop()
        for v, _, et in net.fwd.get(u, ()):
            if et in IDENTITY and v not in seen:
                seen.add(v)
                stack.append(v)
    entered = set()
    for u in seen:
        for v, _, et in net.fwd.get(u, ()):
            if et in ACTS:
                return "A", len(starts), 0
            if et == "input":
                entered.add(v)
    if not entered:
        return "N", len(starts), 0
    cats = defaultdict(list)
    for u, es in net.fwd.items():
        for v, _, et in es:
            if et == "catalyst" and v in entered:
                cats[v].append(u)
    free = [r for r in entered if not any(not has(c) for c in cats[r])]   # no catalyst free of the protein
    return ("B" if free else "C"), len(starts), len(entered)


def pooled(classes):
    return min(classes, key=ORDER.index) if classes else "N"


def main():
    build = Path(sys.argv[1])
    pois = [r for r in A.read_tsv(A.MAPS / "pois.tsv")]
    per = defaultdict(dict)
    for pdir in sorted(p for p in build.iterdir() if (p / "logic_network.csv").exists()):
        net = T.Net(pdir)
        present = {ref.split(":")[-1] for ref in net.ref.values()}
        for r in pois:
            if r["accession"] in present:
                per[r["gene"]][pdir.name] = classify(net, r["accession"])[0]
    print(f"# specs/182: activation class of each Lun protein from build {build.name} (structure only; no Lun value read).")
    print("# A the starting form is a catalyst or regulator; B it first binds or modifies itself; C another enzyme first; N none.")
    print("gene\tclass\tnetworks\tper_network")
    for r in pois:
        d = per.get(r["gene"], {})
        print(f"{r['gene']}\t{pooled(list(d.values()))}\t{len(d)}\t{';'.join(f'{k}:{v}' for k, v in sorted(d.items()))}")


if __name__ == "__main__":
    main()
