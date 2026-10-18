"""specs/070: is every hand-off in the generated networks one the curator intended?

A hand-off is the network passing a reaction's product to another reaction:
reaction A -output-> node n -role-> reaction B (role: input, catalyst,
regulator+/-), directly, or through generator-derived edges between n and the
node B reads (assembly, variant_pool, variant_split, specialization,
set_member). Each hand-off is checked, at the level of stable ids, against three
curator records:

  participants  Neo4j: some output of A and some role-participant of B are the
                same entity ("same"), or related by set membership ("set") or
                composition ("component"), or not related ("none").
  order         precedingEvent within the pathway: B after A ("forward"), A
                after B ("reverse"), both ("cycle"), neither ("unordered").
  diagram       Release layout JSON: some diagram draws A's output and B's
                participant as ONE glyph ("shared"); joins two glyphs by a
                dashed link ("link"); draws both reactions but on separate
                glyphs ("separate"); never draws A and B together ("apart");
                or draws neither ("undrawn").

  handoff_audit.py BUILD_DIR REACTIONS_TSV ENTITIES_TSV DIAGRAM_DIR OUT_TSV

Cofactor nodes (cofactors.csv) are skipped: the solver holds them inert.
"""
import collections
import csv
import glob
import json
import sys
from pathlib import Path

DERIVED = {"assembly", "variant_pool", "variant_split", "specialization", "set_member"}
ROLE_KEYS = {"input": ("inputs",), "catalyst": ("catalysts",),
             "regulator+": ("activators", "catalysts"), "regulator-": ("inhibitors",)}
MAX_DERIVED_HOPS = 3


def _split(s):
    return [x for x in s.split(";") if x]


def _num(stid):
    return int(stid.rsplit("-", 1)[1])


def load_reactions(path):
    by_pathway = collections.defaultdict(dict)
    with open(path) as f:
        for row in csv.reader(f, skipinitialspace=True):
            if row[0] == "pid":
                continue
            pid, rid, cls, ins, outs, cats, pos, neg, prev = row
            by_pathway[pid][rid] = {"cls": cls, "outs": set(_split(outs)), "prev": set(_split(prev)),
                                    "input": set(_split(ins)), "catalyst": set(_split(cats)),
                                    "regulator+": set(_split(pos)), "regulator-": set(_split(neg))}
    return by_pathway


def load_entities(path):
    members, components = {}, {}
    with open(path) as f:
        for row in csv.reader(f, skipinitialspace=True):
            if row[0] == "id":
                continue
            members[row[0]] = set(_split(row[2]))
            components[row[0]] = set(_split(row[3]))
    return members, components


def relation(outs, parts, members, components):
    """Strongest curated relation between any output of A and any participant of B."""
    if outs & parts:
        return "same"
    best = "none"
    for x in outs:
        mx = members.get(x, set()) | {x}
        for y in parts:
            my = members.get(y, set()) | {y}
            if (members.get(x) or members.get(y)) and mx & my:
                return "set"
            if x in components.get(y, ()) or y in components.get(x, ()):
                best = "component"
    return best


def load_diagrams(diagram_dir):
    """reaction dbId -> list of (diagram, {role key: set of glyph ids}, outputs glyph ids);
    per diagram: glyph id -> entity dbId, and dashed links as glyph pairs."""
    by_rxn = collections.defaultdict(list)
    links = {}
    for f in sorted(glob.glob(f"{diagram_dir}/R-HSA-*.json")):
        if f.endswith(".graph.json"):
            continue
        d = json.load(open(f))
        name = Path(f).stem
        glyph = {n["id"]: n.get("reactomeId") for n in d["nodes"]}
        pairs = set()
        for l in d.get("links", []):
            for a in l.get("inputs", []):
                for b in l.get("outputs", []):
                    pairs.add((a["id"], b["id"]))
                    pairs.add((b["id"], a["id"]))
        links[name] = (glyph, pairs)
        for e in d["edges"]:
            if e.get("reactomeId") is None:
                continue
            roles = {k: {x["id"] for x in e.get(k, [])} for k in
                     ("inputs", "outputs", "catalysts", "activators", "inhibitors")}
            by_rxn[e["reactomeId"]].append((name, roles))
    return by_rxn, links


def diagram_status(a, b, role, by_rxn, links, entity_ok):
    """entity_ok(glyph_entity_dbid) -> True when that glyph can carry this hand-off."""
    da = {n: r for n, r in by_rxn.get(_num(a), [])}
    db = {n: r for n, r in by_rxn.get(_num(b), [])}
    if not da and not db:
        return "undrawn"
    common = set(da) & set(db)
    if not common:
        return "apart"
    best = "separate"
    for name in common:
        outs = da[name]["outputs"]
        uses = set().union(*(db[name][k] for k in ROLE_KEYS[role]))
        if outs & uses:
            return "shared"
        glyph, pairs = links[name]
        if any((g, h) in pairs for g in outs for h in uses):
            best = "link"
    return best


def closure(R):
    anc = {}

    def go(r, stack=()):
        if r not in anc:
            s = set()
            for q in R[r]["prev"]:
                if q in R and q not in stack:
                    s.add(q)
                    s |= go(q, stack + (r,))
            anc[r] = s
        return anc[r]
    for r in R:
        go(r)
    return anc


def order(a, b, anc):
    f, r = a in anc.get(b, ()), b in anc.get(a, ())
    return "cycle" if f and r else "forward" if f else "reverse" if r else "unordered"


def pathway_handoffs(P):
    """(A, B, role, via) at stable-id level, with the entity stable id of n."""
    ent, kind = {}, {}
    with open(P / "nodes.csv") as f:
        for r in csv.DictReader(f):
            ent[r["uuid"]] = r["diagram_entity_id"]
            kind[r["uuid"]] = r["node_kind"]
    cof = set()
    if (P / "cofactors.csv").exists():
        with open(P / "cofactors.csv") as f:
            cof = {r["stable_id"] for r in csv.DictReader(f)}
    producers = collections.defaultdict(set)
    uses = collections.defaultdict(set)     # node -> (reaction uuid, role)
    derived = collections.defaultdict(set)  # node -> (node, edge type)
    for name in ("logic_network.csv", "boundary_edges.csv"):
        with open(P / name) as f:
            for r in csv.DictReader(f):
                s, t, et = r["source_id"], r["target_id"], r["edge_type"]
                if et == "output":
                    producers[t].add(s)
                elif et in ("input", "catalyst"):
                    uses[s].add((t, et))
                elif et == "regulator":
                    uses[s].add((t, "regulator+" if r["pos_neg"] == "pos" else "regulator-"))
                elif et in DERIVED:
                    derived[s].add((t, et))
    out = collections.Counter()
    for n, prods in producers.items():
        if ent.get(n) in cof:
            continue
        # walk derived edges from n
        frontier, seen = [(n, ())], {n}
        reach = []
        while frontier:
            m, via = frontier.pop()
            reach.append((m, via))
            if len(via) < MAX_DERIVED_HOPS:
                for t, et in derived.get(m, ()):
                    if t not in seen:
                        seen.add(t)
                        frontier.append((t, via + (et,)))
        for m, via in reach:
            for rb, role in uses.get(m, ()):
                for ra in prods:
                    a, b = ent[ra], ent[rb]
                    if a == b:
                        continue
                    out[(a, b, role, "+".join(sorted(set(via))) or "direct", ent[n], ent[m])] += 1
    return out


def main(build, reactions_tsv, entities_tsv, diagram_dir, out_tsv):
    build = Path(build)
    rx = load_reactions(reactions_tsv)
    members, components = load_entities(entities_tsv)
    by_rxn, links = load_diagrams(diagram_dir)
    rows = []
    for pid in sorted(rx):
        R = rx[pid]
        anc = closure(R)
        for (a, b, role, via, en, em), copies in sorted(pathway_handoffs(build / pid).items()):
            if a not in R or b not in R:
                rows.append((pid, a, b, role, via, en, em, copies, "not_in_pathway", "", ""))
                continue
            rel = relation(R[a]["outs"], R[b][role], members, components)
            rows.append((pid, a, b, role, via, en, em, copies, rel, order(a, b, anc),
                         diagram_status(a, b, role, by_rxn, links, None)))
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pid", "r_from", "r_to", "role", "via", "entity_out", "entity_used", "node_copies",
                    "participants", "order", "diagram"])
        w.writerows(rows)
    print(f"{len(rows)} hand-offs")


if __name__ == "__main__":
    main(*sys.argv[1:])
