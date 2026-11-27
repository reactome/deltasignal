#!/usr/bin/env python3
"""specs/181 part M: every Lun pair with a curated route, by outcome and shortest route length (reaction steps).

Usage: routes_181.py PREDICTIONS BUILD OUT_JSON A_FILE B_FILE [CLASS_DIR]
"""
import csv
import json
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
for d in ("173", "174", "176", "177"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
import score_177 as S  # noqa: E402
import trace_174 as T  # noqa: E402

BANDS = ((0, 0), (1, 5), (6, 10), (11, 20), (21, 10**6))


def outcome(c):
    if not B.changed(c["measured"]):
        return "false alarm" if B.changed(c["model"]) else "no change, none called"
    if not B.changed(c["model"]):
        return "silent"
    return "right" if c["model"] == c["measured"] else "wrong direction"


def band(steps):
    if steps is None:
        return "no path found"
    return next(str(lo) if lo == hi else f"{lo}-{hi}" if hi < 10**6 else f"{lo}+" for lo, hi in BANDS if lo <= steps <= hi)


def distances(net, sources):
    """Reaction steps (reaction nodes passed) on the fewest-edges path from any source to every node."""
    dist = {s: (0, 0) for s in sources}
    q = deque(sources)
    while q:
        u = q.popleft()
        e, r = dist[u]
        for v, _, _ in net.fwd.get(u, ()):
            if v not in dist:
                dist[v] = (e + 1, r + (net.kind.get(v) == "reaction"))
                q.append(v)
    return dist


def route_lengths(pred_file, build, cases):
    """{(gene, readout): reaction steps on the shortest curated route} for `cases` ({(gene, readout): case})."""
    by_pw = defaultdict(lambda: defaultdict(set))
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        k = (r["gene"], r["readout"])
        if k in cases and r["self"] != "1" and r["status"] == "scored" and r["reach_predicted"] == "up":
            by_pw[r["pathway_id"]][r["gene"]].add(r["readout"])
    acc = {r["gene"]: r["accession"] for r in A.read_tsv(A.MAPS / "pois.tsv")}
    ro = {r["readout"]: r for r in A.read_tsv(A.MAPS / "readouts.tsv")}
    best = {}
    for pw in sorted(by_pw):
        net = T.Net(Path(build) / pw)
        targets = {}
        for g, rs in by_pw[pw].items():
            src = net.nodes_with(acc[g])
            dist = distances(net, [u for u in src if net.indeg[u] == 0] or src)
            for r in rs:
                if r not in targets:
                    targets[r] = T.readout_targets(net, ro[r])
                ds = [dist[t][1] for t in targets[r] if t in dist]
                if ds and ((g, r) not in best or min(ds) < best[(g, r)]):
                    best[(g, r)] = min(ds)
    return best


def main():
    import dataset_roles
    pred_file, build, out, a_file, b_file = sys.argv[1:6]
    tk, tyr = S.read_classes(sys.argv[6]) if len(sys.argv) > 6 else S.read_classes()
    genes = A.marker_genes()
    m, p = {}, {}
    for half, path in (("A", a_file), ("B", b_file)):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        m.update(S.measured(path, genes, half))
        p.update(S.predictions(pred_file, half))
    cases = {(c["gene"], c["readout"]): c for c in S.build_cases(m, p, tk, tyr) if c["reach"] == "up"}
    best = route_lengths(pred_file, build, cases)
    tab = defaultdict(Counter)
    for k, c in cases.items():
        tab[band(best.get(k))][outcome(c)] += 1
    res = {"pairs_with_route": len(cases), "by_route_length": {b: dict(v) for b, v in sorted(tab.items())}}
    for b, v in tab.items():
        calls = v["right"] + v["wrong direction"] + v["false alarm"]
        real = v["right"] + v["wrong direction"] + v["silent"]
        res["by_route_length"][b].update({"pairs": sum(v.values()), "share_measured_changed": real / sum(v.values()),
                                          "precision_of_calls": v["right"] / calls if calls else None})
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
