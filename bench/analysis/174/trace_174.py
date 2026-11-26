#!/usr/bin/env python3
"""specs/174: why Lun 2019 half-A cases pass or fail (registered classes), with network traces of the largest groups.

  cases  A_FILE PREDICTIONS OUT_DIR       the specs/173 measured call and frozen predictions per pair, classified;
                                          writes cases.tsv and classes.json (outside the repo); prints counts only
  trace  OUT_DIR BUILD                    shortest signed path, pinned protein -> readout, for every case in the
                                          registered groups; writes traces.tsv and groups.json
  trace-all OUT_DIR BUILD                 the same for every case with a route (post hoc); writes traces_all.tsv

Readout nodes are matched from the catalog's own leaf_refs.csv (UniProt accession and phosphosite position), not Neo4j.
Half B is never read: the reader is specs/173's, which refuses any half-B protein.
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
import lun_a_173 as A  # noqa: E402
from common_173 import split_of  # noqa: E402

BACKGROUND_RATE = 0.25
ROUTE_MIN_PAIRS, ROUTE_MIN_RATE, ROUTE_RATIO = 10, 0.20, 2.0
TOP_READOUTS = 5
GROUP_COVER, MAX_GROUPS = 0.5, 8
FAILURES = ("R2", "R3", "R4", "N1")
PASSES = ("R1", "R5", "N2")

# readout gene -> UniProt accession (canonical entries; checked against the catalog's leaf_refs in `trace`)
READOUT_ACC = {
    "CREB1": "P16220", "ATF1": "P18846", "STAT5A": "P42229", "SRC": "P12931", "PTK2": "Q05397", "MAP2K1": "Q02750",
    "MAP2K2": "P36507", "MAPKAPK2": "P49137", "RPS6KB1": "P23443", "MAP2K4": "P45985", "STAT1": "P42224",
    "TP53": "P04637", "RELA": "Q04206", "MAPK14": "Q16539", "MAPK11": "Q15759", "MAPK12": "P53778", "MAPK13": "O15264",
    "PRKAA1": "Q13131", "AKT1": "P31749", "AKT2": "P31751", "AKT3": "Q9Y243", "MAPK3": "P27361", "MAPK1": "P28482",
    "MARCKS": "P29966", "GSK3B": "P49841", "MAP2K3": "P46734", "MAP2K6": "P52564", "PDPK1": "O15530", "BTK": "Q06187",
    "ITK": "Q08881", "RPS6KA1": "Q15418", "RPS6KA2": "Q15349", "RPS6KA3": "P51812", "SMAD2": "Q15796",
    "SMAD3": "P84022", "CTNNB1": "P35222", "MAPK8": "P45983", "MAPK9": "P45984", "MAPK10": "P53779",
    "PLCG2": "P16885", "RPS6": "P62753", "MTOR": "P42345", "RAF1": "P04049", "EIF4EBP1": "Q13541",
}


def changed(c):
    return c in ("up", "down")


def case_class(reach, frozen, measured):
    f = frozen if changed(frozen) else "not"
    if reach != "up":
        return "N1" if changed(measured) else "N2"
    if changed(f):
        return "R1" if f == measured else "R2" if changed(measured) else "R3"
    return "R4" if changed(measured) else "R5"


def readout_types(cases):
    """{readout: background | route-specific | quiet} by the registered rule, with the rates used."""
    by = defaultdict(list)
    for c in cases:
        by[c["readout"]].append(c)
    out = {}
    for r, cs in by.items():
        on = [c for c in cs if c["reach"] == "up"]
        off = [c for c in cs if c["reach"] != "up"]
        r_on = sum(changed(c["measured"]) for c in on) / len(on) if on else None
        r_off = sum(changed(c["measured"]) for c in off) / len(off) if off else None
        if r_off is not None and r_off >= BACKGROUND_RATE:
            kind = "background"
        elif (len(on) >= ROUTE_MIN_PAIRS and r_on >= ROUTE_MIN_RATE
              and (r_off is None or r_on >= ROUTE_RATIO * r_off)):
            kind = "route-specific"
        else:
            kind = "quiet"
        out[r] = {"type": kind, "route_pairs": len(on), "no_route_pairs": len(off),
                  "rate_route": r_on, "rate_no_route": r_off}
    return out


def build_cases(a_file, pred_file):
    genes = A.marker_genes()
    sl = A.slopes(A.load_a(a_file, genes))
    calls = A.measured_calls(sl, genes, A.thresholds(sl))
    frozen, reach, per = defaultdict(list), defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if r["self"] == "1" or r["status"] != "scored" or split_of(r["gene"]) != "A":
            continue
        k = (r["gene"], r["readout"])
        frozen[k].append(r["predicted"])
        reach[k].append(r["reach_predicted"])
        per[k].append(f"{r['pathway_id']}:{r['predicted']}:{r['reach_predicted']}:{r['pred_ui']}")
    import lun_score as L
    cases = []
    for k in sorted(frozen):
        if k not in calls:
            continue
        f = L.combine_networks(frozen[k])
        rc = L.combine_networks(reach[k])
        cases.append({"gene": k[0], "readout": k[1], "measured": calls[k][0], "E": f"{calls[k][1]:.4g}",
                      "frozen": f, "reach": rc, "cls": case_class(rc, f, calls[k][0]), "networks": ";".join(per[k])})
    return cases


def summarise(cases, types):
    lib = {r["gene"]: r["library"] for r in A.read_tsv(A.MAPS / "pois.tsv")}
    cls = Counter(c["cls"] for c in cases)
    by_type = defaultdict(Counter)
    for c in cases:
        by_type[c["cls"]][types[c["readout"]]["type"]] += 1
    r1 = by_type["R1"]
    wrong = Counter(c["readout"] for c in cases if c["cls"] in ("R2", "R3"))
    top = wrong.most_common(TOP_READOUTS)
    per_readout = {r: dict(Counter(c["cls"] for c in cases if c["readout"] == r)) for r in sorted(types)}
    return {
        "cases": len(cases), "classes": dict(sorted(cls.items())),
        "class_by_readout_type": {k: dict(v) for k, v in sorted(by_type.items())},
        "readout_types": types, "per_readout": per_readout,
        "Q1_R1_on_background": {"R1": sum(r1.values()), "background": r1["background"],
                                "share": r1["background"] / sum(r1.values()) if r1 else None},
        "Q2_concentration": {"R2_R3": sum(wrong.values()), "top5": top,
                             "share_top5": sum(n for _, n in top) / sum(wrong.values()) if wrong else None},
        "class_by_library": {k: dict(Counter(c["cls"] for c in cases if lib.get(c["gene"]) == k)) for k in ("k", "p")},
    }


def groups(cases):
    """Registered groups: per failure class, readouts largest first until half the class is covered (at most eight);
    for passes, R1 on route-specific readouts the same way."""
    out = {}
    for k in FAILURES + ("R1",):
        sel = [c for c in cases if c["cls"] == k and (k != "R1" or c["rtype"] == "route-specific")]
        cnt = Counter(c["readout"] for c in sel)
        picked, covered = [], 0
        for r, n in cnt.most_common():
            if covered >= GROUP_COVER * len(sel) or len(picked) >= MAX_GROUPS:
                break
            picked.append((r, n))
            covered += n
        out[k] = {"cases": len(sel), "groups": picked, "covered": covered}
    return out


# ---------------------------------------------------------------- network tracing

class Net:
    def __init__(self, pdir):
        self.fwd = defaultdict(list)
        self.indeg = Counter()
        with open(pdir / "logic_network.csv") as f:
            for r in csv.DictReader(f):
                self.fwd[r["source_id"]].append((r["target_id"], r["pos_neg"], r["edge_type"]))
                if r["edge_type"] != "depletion":
                    self.indeg[r["target_id"]] += 1
        self.leaves, self.kind, self.stid = {}, {}, {}
        with open(pdir / "nodes.csv") as f:
            for r in csv.DictReader(f):
                self.leaves[r["uuid"]] = [x for x in r["member_leaves"].split("|") if x]
                self.kind[r["uuid"]] = r["node_kind"]
                self.stid[r["uuid"]] = r["diagram_entity_id"]
        self.ref, self.mods = {}, {}
        with open(pdir / "leaf_refs.csv") as f:
            for r in csv.DictReader(f):
                self.ref[r["leaf_stid"]] = r["reference_stid"]
                self.mods[r["leaf_stid"]] = r["modifications"]

    def nodes_with(self, acc, site=None, unmodified=False):
        want = f"uniprot:{acc}"
        hit = []
        for u, ls in self.leaves.items():
            if self.kind[u] == "reaction":
                continue
            for leaf in ls:
                if self.ref.get(leaf) != want:
                    continue
                m = self.mods.get(leaf, "")
                if site is not None and not re.search(rf"phospho[^;]* at {site}\b", m):
                    continue
                if unmodified and m:
                    continue
                hit.append(u)
                break
        return hit

    def path(self, sources, targets):
        """Shortest path (BFS) from any source to any target: [(node, sign_so_far, edge_type)]."""
        targets = set(targets)
        prev = {s: None for s in sources}
        q = deque(sources)
        while q:
            u = q.popleft()
            if u in targets:
                out, x = [], u
                while x is not None:
                    out.append(x)
                    x = prev[x][0] if prev[x] else None
                out.reverse()
                steps, sign = [], 1
                for a, b in zip(out, out[1:]):
                    e = next(e for e in self.fwd[a] if e[0] == b)
                    sign *= -1 if e[1] == "neg" else 1
                    steps.append((b, sign, e[2]))
                return out[0], steps
            for v, pn, et in self.fwd.get(u, ()):
                if v not in prev:
                    prev[v] = (u, pn, et)
                    q.append(v)
        return None, None

    def label(self, u, names):
        if self.kind[u] == "reaction":
            return f"[{self.stid[u]}]"
        gs = sorted({names.get(self.ref.get(leaf, "").split(":")[-1], self.ref.get(leaf, leaf).split(":")[-1])
                     for leaf in self.leaves[u]})
        return ":".join(gs[:6]) + ("…" if len(gs) > 6 else "")


def readout_targets(net, readout_row):
    if readout_row["kind"] == "total":
        acc = READOUT_ACC.get(readout_row["sites"])
        return net.nodes_with(acc, unmodified=True) if acc else []
    hit = []
    for s in readout_row["sites"].split(","):
        g, pos = s.split(":")
        if g in READOUT_ACC:
            hit += net.nodes_with(READOUT_ACC[g], site=int(pos))
    return list(dict.fromkeys(hit))


def trace(out_dir, build, everything=False):
    """`everything` (post hoc, not registered): every case with a route, not only the registered groups."""
    out_dir = Path(out_dir)
    cases = list(csv.DictReader(open(out_dir / "cases.tsv"), delimiter="\t"))
    pois = A.read_tsv(A.MAPS / "pois.tsv")
    acc = {r["gene"]: r["accession"] for r in pois}
    names = {v: k for k, v in READOUT_ACC.items()} | {r["accession"]: r["gene"] for r in pois}
    ro = {r["readout"]: r for r in A.read_tsv(A.MAPS / "readouts.tsv")}
    g = groups(cases)
    want = {(k, r) for k, v in g.items() for r, _ in v["groups"]}
    todo = [c for c in cases if (c["cls"], c["readout"]) in want
            and (c["cls"] != "R1" or c["rtype"] == "route-specific")]
    if everything:
        todo = [c for c in cases if c["reach"] == "up"]
    by_pw = defaultdict(list)
    for c in todo:
        for n in c["networks"].split(";"):
            pw, pred, rp, ui = n.split(":")
            if rp == "up" or c["cls"] == "N1":
                by_pw[pw].append((c, pred, ui))
    rows = []
    for pw in sorted(by_pw):
        net = Net(Path(build) / pw)
        tcache = {}
        for c, pred, ui in by_pw[pw]:
            r = c["readout"]
            if r not in tcache:
                tcache[r] = readout_targets(net, ro[r])
            src = net.nodes_with(acc[c["gene"]])
            roots = [u for u in src if net.indeg[u] == 0] or src
            start, steps = net.path(roots, tcache[r]) if tcache[r] else (None, None)
            rows.append({"cls": c["cls"], "readout": r, "gene": c["gene"], "measured": c["measured"],
                         "frozen": c["frozen"], "pathway": pw, "pred": pred, "pred_ui": ui,
                         "readout_nodes_matched": len(tcache[r]), "hops": "" if steps is None else len(steps),
                         "path_sign": "" if not steps else ("+" if steps[-1][1] > 0 else "-"),
                         "path": "" if steps is None else " ".join(
                             [net.label(start, names)] + [f"-{'|' if s < 0 else '>'}{et}-> {net.label(b, names)}"
                                                          for b, s, et in steps if net.kind[b] != "reaction"])})
    with open(out_dir / ("traces_all.tsv" if everything else "traces.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if not everything:
        (out_dir / "groups.json").write_text(json.dumps(g, indent=1))
        print(json.dumps(g, indent=1))
    print(f"traced rows: {len(rows)}; no matched readout node: {sum(not r['readout_nodes_matched'] for r in rows)}; "
          f"no path: {sum(r['hops'] == '' for r in rows)}")


def main():
    cmd = sys.argv[1]
    import dataset_roles
    dataset_roles.require_role("lun2019_a", "fit")
    if cmd == "cases":
        a_file, pred_file, out_dir = sys.argv[2], sys.argv[3], Path(sys.argv[4])
        cases = build_cases(a_file, pred_file)
        types = readout_types(cases)
        for c in cases:
            c["rtype"] = types[c["readout"]]["type"]
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / "cases.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(cases[0]), delimiter="\t")
            w.writeheader()
            w.writerows(cases)
        s = summarise(cases, types)
        (out_dir / "classes.json").write_text(json.dumps(s, indent=1))
        print(json.dumps({k: s[k] for k in ("cases", "classes", "class_by_readout_type", "Q1_R1_on_background",
                                            "Q2_concentration", "class_by_library")}, indent=1))
    elif cmd in ("trace", "trace-all"):
        trace(sys.argv[2], sys.argv[3], everything=cmd == "trace-all")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
