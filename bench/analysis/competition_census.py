"""specs/050 research R1 census: would competitive allocation of a shared consumed species
predict the curated direction? Analysis only.

Motif: entity S is an `input` to >= 2 distinct reactions (by edge_reaction_id);
S is not a cofactor, not a carrier-loop free form, not a pool node. The
perturbed gene's pinned node is a partner of consumer Ri (another input, or a
catalyst, of Ri's reaction node). The readout is reachable from a sibling
consumer Rj and NOT from Ri's own reaction node. Competition then predicts the
readout moves opposite to the gene: knockout -> UP, overexpression -> DOWN.

Usage: competition_census.py BUILD_DIR CASES_TSV
"""
import csv, sys, os, collections

build, cases_tsv = sys.argv[1], sys.argv[2]
UP, NORMAL, DOWN = "2", "1", "0"
SINK = os.environ.get("SINK") == "1"
NEGREG = None
if os.environ.get("NEGREG"):
    NEGREG = {l.split(",")[0].strip().strip('"') for l in open(os.environ["NEGREG"]) if l.startswith('"')}

def load(pdir):
    edges = list(csv.DictReader(open(f"{pdir}/logic_network.csv")))
    out = collections.defaultdict(set)
    for e in edges:
        out[e["source_id"]].add(e["target_id"])
    excluded = set()
    for f, col in (("carrier_loops.csv", "free_uuid"), ("pools.csv", "node_uuid")):
        p = f"{pdir}/{f}"
        if os.path.exists(p):
            excluded |= {r[col] for r in csv.DictReader(open(p))}
    cof = set()
    p = f"{pdir}/cofactors.csv"
    if os.path.exists(p):
        cofst = {list(r.values())[0] for r in csv.DictReader(open(p))}
        for r in csv.DictReader(open(f"{pdir}/stid_to_uuid_mapping.csv")):
            if r["stable_id"].split("::")[0] in cofst:
                cof.add(r["uuid"])
    # consumers: S -> {reaction stid: set(reaction nodes)}
    cons = collections.defaultdict(lambda: collections.defaultdict(set))
    partners = collections.defaultdict(set)   # reaction node -> partner entity nodes
    for e in edges:
        if e["edge_type"] in ("input", "catalyst"):
            partners[e["target_id"]].add(e["source_id"])
        if e["edge_type"] == "input":
            cons[e["source_id"]][e["edge_reaction_id"]].add(e["target_id"])
    acts = {e["source_id"] for e in edges if e["edge_type"] in ("catalyst", "regulator")}
    shared = {s: rs for s, rs in cons.items()
              if len(rs) >= 2 and s not in excluded and s not in cof}
    # closed_twin: a closed carrier free form follows its gene's root copy, so it
    # counts as the gene (same reactome_id)
    stid = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(f"{pdir}/stid_to_uuid_mapping.csv"))}
    closed_free = collections.defaultdict(set)
    p = f"{pdir}/carrier_loops.csv"
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            if r["closed"] == "True":
                closed_free[stid.get(r["free_uuid"])].add(r["free_uuid"])
    twin_of = {u: closed_free[st] for u, st in stid.items() if st in closed_free}
    return out, shared, partners, cof, acts, twin_of

def reach(out, starts, cache={}):
    seen, stack = set(starts), list(starts)
    while stack:
        u = stack.pop()
        for v in out.get(u, ()):
            if v not in seen:
                seen.add(v); stack.append(v)
    return seen

ids = {}
for r in csv.DictReader(open(f"{build}/pathways.tsv"), delimiter="\t"):
    ids[list(r.values())[1]] = list(r.values())[0]

tally = collections.Counter()
examples = collections.defaultdict(list)
nets = {}
for c in csv.DictReader(open(cases_tsv), delimiter="\t"):
    if c["valid"] != "1" or c["direction"] == NORMAL:
        continue
    pid = ids.get(c["pathway"])
    if not pid or not os.path.isdir(f"{build}/{pid}"):
        continue
    if pid not in nets:
        nets[pid] = load(f"{build}/{pid}")
        nets[pid] += ({},)
    out, shared, partners, cof, acts, twin_of, rcache = nets[pid]
    genes = set(filter(None, c["gene_uuids"].split("|")))
    genes |= {f for g in list(genes) for f in twin_of.get(g, ())}
    outs = set(filter(None, c["output_uuids"].split("|")))
    hit = False
    for s, rs in shared.items():
        own = [st for st, nodes in rs.items() if any(genes & partners[n] for n in nodes)
               and (NEGREG is None or st in NEGREG)]
        if not own or s in genes:
            continue
        for st in own:
            def r_of(st_):
                k = (s, st_)
                if k not in rcache:
                    rcache[k] = reach(out, rs[st_])
                return rcache[k]
            mine = r_of(st)
            if outs & mine:
                continue
            # SINK=1: the gene's branch is inert -- nothing downstream of it acts
            # (catalyst or regulator edge out) on anything
            if SINK and mine & acts:
                continue
            if any(outs & r_of(o) for o in rs if o != st):
                hit = True
    if not hit:
        continue
    # does the gene reach the readout at all other than through a competition?
    gk = ("gene", frozenset(genes))
    if gk not in rcache:
        rcache[gk] = reach(out, genes)
    direct = bool(outs & rcache[gk])
    exp_kind = {"1": "exp NORMAL", c["direction"]: "exp same as gene"}.get(c["expected"], "exp opposite")
    tally[f"  detail: {'direct route' if direct else 'no direct route'}, {exp_kind}, {'right' if c['predicted']==c['expected'] else 'wrong'}"] += 1
    comp = UP if c["direction"] == DOWN else DOWN
    right = c["predicted"] == c["expected"]
    agree = comp == c["expected"]
    key = ("now right" if right else "now wrong") + (", competition agrees" if agree else ", competition disagrees")
    tally[key] += 1
    tally[c["pathway"] + " | " + key] += 1
    if not right and agree and len(examples[c["pathway"]]) < 3:
        examples[c["pathway"]].append(f'{c["gene"]} {c["direction"]} -> {c["key_output"]}')
for k in sorted(k for k in tally if "|" not in k):
    print(f"{tally[k]:6d}  {k}")
print("\nby pathway (now wrong & agrees / now right & disagrees):")
pw = collections.Counter(); pb = collections.Counter()
for k, v in tally.items():
    if "|" in k:
        p, kk = k.split(" | ")
        if kk == "now wrong, competition agrees": pw[p] += v
        if kk == "now right, competition disagrees": pb[p] += v
for p in sorted(set(pw) | set(pb), key=lambda p: -(pw[p] - pb[p])):
    print(f"  {pw[p]:4d} / {pb[p]:4d}  {p}  {examples.get(p, [])[:2]}")
