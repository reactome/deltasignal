"""specs/070: curated hand-offs the network does not make.

A curated hand-off: within one catalog pathway, reaction A outputs entity E and
reaction B uses E (input, catalyst or regulator), exactly (not via sets).
Cofactors excluded. Each is looked up in handoff_audit.py's output; if absent
the diagram status is reported.

  missing_handoffs.py HANDOFFS_TSV REACTIONS_TSV BUILD_DIR DIAGRAM_DIR OUT_TSV
"""
import collections, csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from handoff_audit import load_reactions, load_diagrams, diagram_status, closure, order

def main(handoffs, reactions_tsv, build, diagram_dir, out_tsv):
    have = set()
    with open(handoffs) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            have.add((r["pid"], r["r_from"], r["r_to"], r["role"]))
    rx = load_reactions(reactions_tsv)
    by_rxn, links = load_diagrams(diagram_dir)
    rows = []
    for pid, R in sorted(rx.items()):
        cof = set()
        cf = Path(build) / pid / "cofactors.csv"
        if cf.exists():
            cof = {r["stable_id"] for r in csv.DictReader(open(cf))}
        anc = closure(R)
        users = collections.defaultdict(set)
        for b, d in R.items():
            for role in ("input", "catalyst", "regulator+", "regulator-"):
                for e in d[role]:
                    users[e].add((b, role))
        for a, d in R.items():
            for e in d["outs"] - cof:
                for b, role in users[e]:
                    if b == a:
                        continue
                    present = (pid, a, b, role) in have
                    rows.append((pid, a, b, role, e, "present" if present else "MISSING",
                                 order(a, b, anc), diagram_status(a, b, role, by_rxn, links, None)))
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pid", "r_from", "r_to", "role", "entity", "network", "order", "diagram"])
        w.writerows(rows)
    print(collections.Counter((r[5], r[7]) for r in rows).most_common())

if __name__ == "__main__":
    main(*sys.argv[1:])
