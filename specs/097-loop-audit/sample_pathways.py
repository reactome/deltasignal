"""specs/097: a deterministic sample of Reactome human pathways outside the 92
benchmark pathways, for the loop audit (no ground truth exists for them).

Candidates: pathways under a human top-level pathway that have their own diagram
and at least 5 direct reactions, not among the 92, outside Disease and Drug ADME
(mutant and drug variants; drugs are held inert by the solver). Each candidate is
assigned to the first of its top-level pathways by name. N pathways are allocated
to the areas in proportion to their candidates (at least 1 each), and drawn
evenly spaced through each area's candidates sorted by stId.
Run from the generator's poetry env with NEO4J_PASSWORD set.
Usage: python3 sample_pathways.py <catalog_pathways.tsv> <N> > sample.tsv
"""
import csv, sys
sys.path.insert(0, "/home/awright/gitroot/logic-network-generator")
from src.neo4j_connector import get_graph

bench = {r["id"] for r in csv.DictReader((l for l in open(sys.argv[1]) if not l.startswith("#")), delimiter="\t")}
n_total = int(sys.argv[2])
rows = get_graph().run("""
MATCH (t:TopLevelPathway {speciesName:'Homo sapiens'})-[:hasEvent*]->(p:Pathway)
WHERE p.hasDiagram = true AND NOT t.displayName IN ['Disease', 'Drug ADME']
WITH t, p, size([(p)-[:hasEvent]->(r:ReactionLikeEvent) | r]) AS nr
WHERE nr >= 5
RETURN p.stId AS id, p.displayName AS name, t.displayName AS area""").data()
area_of = {}
for r in sorted(rows, key=lambda r: (r["id"], r["area"])):
    if r["id"] not in bench and r["id"] not in area_of:
        area_of[r["id"]] = (r["area"], r["name"])
by_area = {}
for pid, (area, name) in area_of.items():
    by_area.setdefault(area, []).append((pid, name))
total = sum(len(v) for v in by_area.values())
alloc = {a: max(1, round(n_total * len(v) / total)) for a, v in by_area.items()}
print("id\tpathway_name")
for a in sorted(by_area):
    c = sorted(by_area[a]); k = min(alloc[a], len(c))
    for i in range(k):
        pid, name = c[(i * len(c)) // k]
        print(f"{pid}\t{name.replace(' ', '_').replace('/', '_')}")
print(f"# {total} candidates; {sum(min(alloc[a], len(by_area[a])) for a in by_area)} drawn", file=sys.stderr)
