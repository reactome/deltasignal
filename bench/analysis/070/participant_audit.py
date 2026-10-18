"""specs/070: does each reaction node read and write exactly its curated participants?

For every reaction in a catalog pathway, the network's participants (the
diagram_entity_id of every node joined to any copy of the reaction by an
input / output / catalyst / regulator edge) are compared with Neo4j's, role
by role. A network participant counts as curated when it is the curated
entity, a set member of it, or (for derived decomposition) a component of it.

  participant_audit.py BUILD_DIR REACTIONS_TSV ENTITIES_TSV OUT_TSV
"""
import collections, csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from handoff_audit import load_reactions, load_entities

ROLE = {"input": "input", "output": "outs", "catalyst": "catalyst"}


def main(build, reactions_tsv, entities_tsv, out_tsv):
    build = Path(build)
    rx = load_reactions(reactions_tsv)
    members, components = load_entities(entities_tsv)
    rows = []
    for pid, R in sorted(rx.items()):
        P = build / pid
        ent = {}
        with open(P / "nodes.csv") as f:
            for r in csv.DictReader(f):
                ent[r["uuid"]] = r["diagram_entity_id"]
        net = collections.defaultdict(lambda: collections.defaultdict(set))
        with open(P / "logic_network.csv") as f:
            for r in csv.DictReader(f):
                et = r["edge_type"]
                if et == "input":
                    net[ent[r["target_id"]]]["input"].add(ent[r["source_id"]])
                elif et == "output":
                    net[ent[r["source_id"]]]["outs"].add(ent[r["target_id"]])
                elif et == "catalyst":
                    net[ent[r["target_id"]]]["catalyst"].add(ent[r["source_id"]])
                elif et == "regulator":
                    k = "regulator+" if r["pos_neg"] == "pos" else "regulator-"
                    net[ent[r["target_id"]]][k].add(ent[r["source_id"]])
        for rid, d in sorted(R.items()):
            if rid not in net:
                rows.append((pid, rid, "ALL", "reaction_absent", "", ""))
                continue
            for role in ("input", "outs", "catalyst", "regulator+", "regulator-"):
                cur, got = d[role], net[rid][role]
                def covers(x, pool):
                    for y in pool:
                        if x == y or x in members.get(y, ()) or y in members.get(x, ()) \
                                or (members.get(x, set()) & members.get(y, set())):
                            return y
                    return None
                for x in sorted(got):
                    y = covers(x, cur)
                    if y is None:
                        comp = any(x in components.get(c, ()) for c in cur)
                        rows.append((pid, rid, role, "extra_component" if comp else "extra", x, ""))
                for y in sorted(cur):
                    if not any(covers(x, {y}) for x in got):
                        rows.append((pid, rid, role, "missing", "", y))
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pid", "reaction", "role", "kind", "network_entity", "curated_entity"])
        w.writerows(rows)
    print(collections.Counter((r[2], r[3]) for r in rows).most_common())


if __name__ == "__main__":
    main(*sys.argv[1:])
