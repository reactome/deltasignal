"""specs/099: classify each shortest material cycle (material_cycles.jl) by what
moves around it.
  self_loop       - one entity: a reaction outputs one of its own inputs unchanged
  interconversion - a protein leaf is in every entity on the cycle (one protein
                    changing state: a GTPase cycle, a modification cycle the
                    generator's pool finder did not take)
  metabolic       - every leaf on the cycle is a small molecule
  exchange        - anything else: material passed between different entities
Leaf types from Neo4j (EntityWithAccessionedSequence = protein; SimpleEntity and
ChemicalDrug = small molecule). Run from the generator's poetry env with
NEO4J_PASSWORD set. Usage: python3 classify_cycles.py cycles_X.tsv [...]
"""
import csv, sys
from collections import Counter
sys.path.insert(0, "/home/awright/gitroot/logic-network-generator")
from src.neo4j_connector import get_graph

files = sys.argv[1:]
rows = {f: list(csv.DictReader(open(f), delimiter="\t")) for f in files}
leaves = {l for rs in rows.values() for r in rs for part in r["leaves"].split(";") for l in part.split("|") if l}
kinds = {x["s"]: x["k"] for x in get_graph().run(
    "UNWIND $l AS s MATCH (e:PhysicalEntity {stId: s}) RETURN s, e.schemaClass AS k", l=sorted(leaves)).data()}
SMALL = {"SimpleEntity", "ChemicalDrug", "ProteinDrug"}
for f, rs in rows.items():
    out = []
    for r in rs:
        ents = [e for e in r["entities"].split("|") if e]
        ls = [set(x for x in part.split("|") if x) for part in r["leaves"].split(";")]
        if len(set(ents)) <= 1:
            k = "self_loop"
        else:
            common = set.intersection(*ls) if ls else set()
            if any(kinds.get(l) == "EntityWithAccessionedSequence" for l in common):
                k = "interconversion"
            elif all(kinds.get(l) in SMALL for s in ls for l in s):
                k = "metabolic"
            else:
                k = "exchange"
        out.append({**r, "class": k})
    with open(f.replace(".tsv", "_classified.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(out)
    c = Counter(x["class"] for x in out)
    nodes = Counter(); pws = {}
    for x in out:
        nodes[x["class"]] += int(x["comp_size"]); pws.setdefault(x["class"], set()).add(x["pathway"])
    print(f, {k: f"{v} cycles, {len(pws[k])} pathways" for k, v in c.items()})
