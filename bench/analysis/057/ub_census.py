"""specs/057 census: depletion edges that a ubiquitination aims at a partner
passing through it unchanged.

Inputs: the output of ub_reactions.cypher (cypher-shell --format plain), and a
catalog build. For each ubiquitination (Ub an input) in a catalog pathway, a
substrate under the current generator rule is an input leaf whose gene is in
the output. Under LNG_UB_SUBSTRATE_CHANGED it must also leave in a changed
form: a leaf whose stId appears among the output's leaves is passed through.
Prints the depletion edges whose target is such a passed-through leaf (and not
a changed substrate of another ubiquitination in the same pathway).

    python3 bench/analysis/057/ub_census.py ub.txt ~/deltasignal-catalogs/builds/<id>
"""
import collections
import csv
import json
import os
import sys


def main(ub_txt, build):
    passed, changed = collections.defaultdict(set), collections.defaultdict(set)
    nrx = nrx_passed = 0
    for line in open(ub_txt).read().splitlines()[1:]:
        pid, _rx, ins, outleaves, outgenes = json.loads("[" + line.replace("NULL", "null") + "]")
        og = {g for gs in outgenes for g in (gs or [])}
        ol = set(outleaves)
        nrx += 1
        hit = False
        for leaf, genes in ins:
            if not genes or not og & set(genes):
                continue
            if leaf in ol:
                passed[pid].add(leaf)
                hit = True
            else:
                changed[pid].add(leaf)
        nrx_passed += hit
    print(f"ubiquitinations: {nrx}; with a passed-through partner counted as substrate: {nrx_passed}")
    hits, total = [], 0
    for pid in sorted(set(passed) | set(changed)):
        d = os.path.join(build, pid)
        if not os.path.isdir(d):
            continue
        st = {r[0]: r[1] for r in csv.reader(open(os.path.join(d, "stid_to_uuid_mapping.csv")))}
        for f in ("logic_network.csv", "boundary_edges.csv"):
            p = os.path.join(d, f)
            if not os.path.exists(p):
                continue
            for r in csv.DictReader(open(p)):
                if r["edge_type"] != "depletion":
                    continue
                t = st.get(r["target_id"])
                if t in passed[pid] and t not in changed[pid]:
                    hits.append((pid, st.get(r["source_id"]), t))
    per = collections.Counter(h[0] for h in hits)
    print(f"depletion edges onto a passed-through partner: {len(hits)} in {len(per)} pathways")
    for p, n in per.most_common():
        print(f"  {n:4d} {p}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
