"""specs/099: is each unhandled interconversion cycle a carrier loop the generator
listed (and the solver kept, e.g. a displacer-driven open release under fed_single),
or listed nowhere (the pool and carrier-loop finders both missed it)?
Usage: python3 listed.py BUILD cycles_X_classified.tsv"""
import csv, os, sys
build, f = sys.argv[1:]
from collections import Counter
c = Counter(); ex = {}
for r in csv.DictReader(open(f), delimiter="\t"):
    if r["class"] != "interconversion":
        continue
    d = os.path.join(build, r["pathway"])
    loops = list(csv.DictReader(open(os.path.join(d, "carrier_loops.csv")))) if os.path.exists(os.path.join(d, "carrier_loops.csv")) else []
    free = {l["free_stid"] for l in loops}
    closed = {l["free_stid"] for l in loops if l["closed"] == "True"}
    ents = set(r["entities"].split("|"))
    k = "listed closed loop" if ents & closed else "listed open loop" if ents & free else "not listed"
    c[k] += 1; ex.setdefault(k, set()).add(r["pathway"])
print(f.split("/")[-1], {k: f"{v} cycles, {len(ex[k])} pathways: {sorted(ex[k])}" for k, v in c.items()})
