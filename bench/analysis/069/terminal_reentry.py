"""specs/069 rule T (terminal re-entry): census and catalog transform.

A material edge re-enters when a reaction with no follower in its pathway
(nothing in the pathway lists it as precedingEvent) outputs an entity that a
reaction PRECEDING it (transitively, within the pathway) consumes as an input.
The pathway's end then feeds its own start, a cycle the curated event order
does not contain. Exempt: the generator's carrier-loop free forms (closed by
DS_CARRIER_LOOP_MODE, specs/049) and cofactors (held inert, specs/007).

  census  REACTIONS_TSV BUILD_DIR OUT_TSV
      REACTIONS_TSV is reactions.cypher's output for the catalog pathways
      (pid, r.stId, cls, ins, outs, prev; lists joined by ';').
  apply   SRC_BUILD DST_BUILD CUTS_TSV [output MEMBERS_TSV]
      copies SRC_BUILD (without results/) and drops edges in
      logic_network.csv and boundary_edges.csv. Default (rule T as
      registered, refuted): every `input` edge whose source node's diagram
      entity is the cut entity and whose edge_reaction_id is the earlier
      reaction. `output` (amendment 1): every `output` edge of the terminal
      reaction (edge_reaction_id = r_late) into a node whose diagram entity
      is the cut entity or one of its set members (MEMBERS_TSV, from
      members.cypher). Prints the per-pathway count of dropped edges.
"""
import collections
import csv
import json
import shutil
import sys
from pathlib import Path


def _split(s):
    return {x for x in s.split(";") if x}


def load_reactions(path):
    by_pathway = collections.defaultdict(dict)
    with open(path) as f:
        for row in csv.reader(f, skipinitialspace=True):
            if row[0] == "pid":
                continue
            pid, r, cls, ins, outs, prev = row
            by_pathway[pid][r] = dict(cls=cls, ins=_split(ins), outs=_split(outs), prev=_split(prev))
    return by_pathway


def _column(path, col):
    if not path.exists():
        return set()
    with open(path) as f:
        return {row[col] for row in csv.DictReader(f) if row.get(col)}


def census(reactions_tsv, build, out_tsv):
    build = Path(build)
    rows = []
    for pid, R in sorted(load_reactions(reactions_tsv).items()):
        ancestors = {}

        def anc(r, stack=()):
            if r not in ancestors:
                s = set()
                for q in R[r]["prev"]:
                    if q in R and q not in stack:
                        s.add(q)
                        s |= anc(q, stack + (r,))
                ancestors[r] = s
            return ancestors[r]

        followers = collections.Counter(q for r in R for q in R[r]["prev"] if q in R)
        consumers = collections.defaultdict(set)
        for r, d in R.items():
            for e in d["ins"]:
                consumers[e].add(r)
        exempt = (_column(build / pid / "carrier_loops.csv", "free_stid")
                  | _column(build / pid / "cofactors.csv", "stable_id"))
        seen = set()
        for late in sorted(R):
            if followers[late]:
                continue
            for e in sorted(R[late]["outs"] - exempt):
                for early in sorted(consumers[e]):
                    if early != late and early in anc(late) and late not in anc(early) \
                            and (e, early) not in seen:
                        seen.add((e, early))
                        rows.append((pid, e, early, late))
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pid", "entity", "r_early", "r_late"])
        w.writerows(rows)
    per = collections.Counter(r[0] for r in rows)
    print(f"{len(rows)} (entity, earlier reaction) pairs in {len(per)} pathways")


def apply(src, dst, cuts_tsv, mode="input", members_tsv=None):
    src, dst = Path(src), Path(dst)
    if dst.exists():
        sys.exit(f"{dst} exists")
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("results"))
    if mode not in ("input", "output"):
        sys.exit(f"unknown mode {mode}")
    members = collections.defaultdict(set)
    if mode == "output":
        with open(members_tsv) as f:
            for row in csv.reader(f, skipinitialspace=True):
                if row[0] != "id":
                    members[row[0]] = _split(row[1])
    cuts = collections.defaultdict(set)
    with open(cuts_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if mode == "input":
                cuts[row["pid"]].add((row["entity"], row["r_early"], "input"))
            else:
                for e in {row["entity"]} | members[row["entity"]]:
                    cuts[row["pid"]].add((e, row["r_late"], "output"))
    dropped = {}
    for pid, pairs in sorted(cuts.items()):
        P = dst / pid
        with open(P / "nodes.csv") as f:
            ent = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(f)}
        n = 0
        for name in ("logic_network.csv", "boundary_edges.csv"):
            with open(P / name) as f:
                reader = csv.DictReader(f)
                fields, rows = reader.fieldnames, list(reader)
            end = "source_id" if mode == "input" else "target_id"
            keep = [r for r in rows
                    if (ent.get(r[end]), r["edge_reaction_id"], r["edge_type"]) not in pairs]
            n += len(rows) - len(keep)
            with open(P / name, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields)
                w.writeheader()
                w.writerows(keep)
        dropped[pid] = n
    meta = json.loads((dst / "BUILD.json").read_text())
    meta.update(build_id=dst.name, variant=f"terminal_reentry_{mode}",
                variant_note=f"specs/069 rule T ({mode} cut) applied to {src.name} by bench/analysis/069/terminal_reentry.py",
                variant_edges_dropped=dropped)
    (dst / "BUILD.json").write_text(json.dumps(meta, indent=1))
    print(f"{sum(dropped.values())} edges dropped in {sum(1 for v in dropped.values() if v)} pathways")
    for pid, n in sorted(dropped.items(), key=lambda kv: -kv[1]):
        print(pid, n)


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"census": census, "apply": apply}[cmd](*args)
