"""specs/050 step 2: motif census. Analysis only; it changes no prediction.

Per pathway of a catalog build, how often each registered motif occurs, as the
solver sees it (one unperturbed solve through the running API, so the counts
are what the default model actually applied), next to the canonical accuracy
on both axes and the specs that name the pathway (a proxy for "inspected when
a rule was written"). The point is to find where a motif is common and nothing
was traced: the motif-level hold-out of specs/050.

Usage:
  motif_census.py BUILD_DIR RESULTS_DIR OUT_TSV [--negreg NEGREG_CSV]

BUILD_DIR    a catalog build (must be what the API at DS_API serves)
RESULTS_DIR  builds/<id>/results/<solver-sha> (curator_cases.tsv, experimental_cases.tsv)
NEGREG_CSV   reaction stIds annotated GO "negative regulation of ..." (specs/050 R1)
"""
import argparse
import collections
import csv
import json
import os
import re
import sys
from pathlib import Path
from urllib.request import Request, urlopen

DS_API = os.environ.get("DS_API", "http://localhost:8080")
REPO = Path(__file__).resolve().parents[2]

SOLVE_COUNTS = ["cycle_pools_solved", "cycle_pools_multistep", "carrier_loop_forms",
                "carrier_loop_edges_dropped", "carrier_loop_twin_edges", "self_inhibitors",
                "self_inhibitor_leaf_pairs", "pathogens_held"]


def post(path, payload):
    req = Request(f"{DS_API}{path}", data=json.dumps(payload).encode(),
                  headers={"Content-Type": "application/json"}, method="POST")
    with urlopen(req, timeout=600) as r:
        return json.loads(r.read())


def largest_scc(edges):
    """Size of the largest strongly connected component (iterative Tarjan)."""
    adj = collections.defaultdict(list)
    nodes = set()
    for s, t in edges:
        adj[s].append(t)
        nodes.update((s, t))
    index, low, on, stack, best, counter = {}, {}, set(), [], 0, [0]
    for root in nodes:
        if root in index:
            continue
        work = [(root, 0)]
        while work:
            v, i = work.pop()
            if i == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on.add(v)
            nbrs = adj[v]
            if i < len(nbrs):
                work.append((v, i + 1))
                w = nbrs[i]
                if w not in index:
                    work.append((w, 0))
                elif w in on:
                    low[v] = min(low[v], index[w])
                continue
            if low[v] == index[v]:
                size = 0
                while True:
                    w = stack.pop()
                    on.discard(w)
                    size += 1
                    if w == v:
                        break
                best = max(best, size)
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[v])
    return best


def competition_sites(edges, excluded, negreg):
    """specs/050 R1 motif: a species consumed (`input`) by >= 2 distinct reactions,
    at least one of them GO-annotated negative regulation; not a carrier free form
    or pool node. Counts species."""
    cons = collections.defaultdict(set)
    for e in edges:
        if e["edge_type"] == "input":
            cons[e["source_id"]].add(e["edge_reaction_id"])
    return sum(1 for s, rs in cons.items()
               if len(rs) >= 2 and s not in excluded and negreg and rs & negreg)


def case_stats(path):
    stats = collections.defaultdict(lambda: [0, 0])
    if not os.path.exists(path):
        return stats
    for c in csv.DictReader(open(path), delimiter="\t"):
        if c["valid"] != "1":
            continue
        stats[c["pathway"]][0] += 1
        stats[c["pathway"]][1] += c["predicted"] == c["expected"]
    return stats


def spec_mentions(names):
    """Which specs/NNN name each pathway (by display name or stId)."""
    texts = {}
    for d in sorted((REPO / "specs").iterdir()):
        if d.is_dir() and re.match(r"\d{3}-", d.name):
            texts[d.name[:3]] = " ".join(p.read_text(errors="ignore") for p in d.glob("*.md"))
    out = {}
    for pid, name in names.items():
        human = name.replace("_", " ")
        out[pid] = [n for n, t in texts.items() if pid in t or human in t or name in t]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("results")
    ap.add_argument("out")
    ap.add_argument("--negreg")
    a = ap.parse_args()
    build = Path(a.build)
    health = json.loads(urlopen(f"{DS_API}/api/health", timeout=30).read())
    served = (health.get("catalog") or {}).get("build_id")
    if served != build.resolve().name:
        sys.exit(f"the API serves {served!r}, not {build.resolve().name!r}")
    negreg = set()
    if a.negreg:
        negreg = {l.split(",")[0].strip().strip('"') for l in open(a.negreg) if l.startswith('"')}
    names = {}
    for r in csv.DictReader(open(build / "pathways.tsv"), delimiter="\t"):
        v = list(r.values())
        names[v[0]] = v[1]
    cur = case_stats(f"{a.results}/curator_cases.tsv")
    exp = case_stats(f"{a.results}/experimental_cases.tsv")
    specs = spec_mentions(names)
    cols = (["pathway_id", "pathway", "nodes", "edges", "largest_scc"] + SOLVE_COUNTS +
            ["drug_stids", "cofactor_stids", "set_pool_edges", "competition_sites",
             "curator_cases", "curator_acc", "experimental_cases", "experimental_acc", "specs"])
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for pid, name in sorted(names.items(), key=lambda x: x[1]):
            parsed = post("/api/parse", {"pathway_id": pid})
            solved = post("/api/solve", {"network_id": parsed["network_id"], "observations": {}})
            if solved.get("status") != "success":
                sys.exit(f"{pid}: unperturbed solve failed: {solved.get('message')}")
            edges = parsed["edges"]
            pdir = build / pid
            excluded = set()
            for fn, col in (("carrier_loops.csv", "free_uuid"), ("pools.csv", "node_uuid")):
                if (pdir / fn).exists():
                    excluded |= {r[col] for r in csv.DictReader(open(pdir / fn))}
            raw = list(csv.DictReader(open(pdir / "logic_network.csv")))
            et = collections.Counter(e.get("edge_type") for e in edges)
            c, x = cur.get(name, [0, 0]), exp.get(name, [0, 0])
            w.writerow([pid, name, len(parsed["nodes"]), len(edges),
                        largest_scc((str(e["parent_uuid"]), str(e["child_uuid"])) for e in edges)] +
                       [solved.get(k, 0) for k in SOLVE_COUNTS] +
                       [len(parsed.get("drug_stids") or []), len(parsed.get("cofactor_stids") or []),
                        et.get("set_member", 0) + et.get("variant_pool", 0),
                        competition_sites(raw, excluded, negreg),
                        c[0], f"{c[1] / c[0]:.3f}" if c[0] else "",
                        x[0], f"{x[1] / x[0]:.3f}" if x[0] else "",
                        ",".join(specs[pid])])
            print(pid, name, flush=True)


if __name__ == "__main__":
    main()
