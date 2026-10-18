"""specs/070: why is each curated hand-off missing? (reads curated_handoffs_r97.tsv)

For a missing (A, E, B): A's output nodes of E (any node A outputs whose
entity is E) and B's nodes of E (any node joined to B's reaction nodes whose
entity is E). Mechanism:
  seam   B reads E through a variant_pool / set_variant node (specs/046
         member variants) while A outputs E as one node
  copy   both are plain nodes, different copies (positional silo)
  other  anything else (printed)

  missing_mechanism.py BUILD_DIR CURATED_HANDOFFS_TSV OUT_TSV
"""
import collections, csv, sys
from pathlib import Path

def main(build, curated, out_tsv):
    miss = collections.defaultdict(list)
    with open(curated) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["network"] == "MISSING":
                miss[r["pid"]].append(r)
    rows = []
    for pid, rs in sorted(miss.items()):
        P = Path(build) / pid
        ent, kind = {}, {}
        for r in csv.DictReader(open(P / "nodes.csv")):
            ent[r["uuid"]], kind[r["uuid"]] = r["diagram_entity_id"], r["node_kind"]
        outs = collections.defaultdict(set)   # (reaction stid, entity) -> nodes
        uses = collections.defaultdict(set)
        outdeg = collections.Counter()
        for name in ("logic_network.csv", "boundary_edges.csv"):
            for r in csv.DictReader(open(P / name)):
                s, t, et = r["source_id"], r["target_id"], r["edge_type"]
                outdeg[s] += 1
                if et == "output":
                    outs[(ent[s], ent[t])].add(t)
                elif et in ("input", "catalyst", "regulator"):
                    uses[(ent[t], ent[s])].add(s)
        for r in rs:
            a_nodes = outs.get((r["r_from"], r["entity"]), set())
            b_nodes = uses.get((r["r_to"], r["entity"]), set())
            bk = {kind[n] for n in b_nodes}
            if not a_nodes:
                mech = "a_output_absent"
            elif not b_nodes:
                mech = "b_use_absent"
            elif bk & {"variant_pool", "set_variant", "set_pool"}:
                mech = "seam"
            else:
                mech = "copy"
            dead = all(outdeg[n] == 0 for n in a_nodes) if a_nodes else ""
            rows.append((pid, r["r_from"], r["r_to"], r["role"], r["entity"], r["order"],
                         r["diagram"], mech, "|".join(sorted(bk)), dead))
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pid", "r_from", "r_to", "role", "entity", "order", "diagram",
                    "mechanism", "b_node_kinds", "a_output_dead_end"])
        w.writerows(rows)
    print(collections.Counter((r[6], r[7]) for r in rows).most_common())
    print("A output dead end:", collections.Counter((r[7], r[9]) for r in rows).most_common())

if __name__ == "__main__":
    main(*sys.argv[1:])
