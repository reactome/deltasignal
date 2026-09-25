#!/usr/bin/env python3
"""Classify benchmark failures by the network structure between the pinned root
inputs and the readout.

The benchmark perturbs root inputs and scores terminal readouts; if the ends
are right, the middle is taken to be right. So a failure is a question about
the region in between: the nodes on some path from a pinned root to the
readout. For each scored case this computes, on that region:

  reach        no path / a path exists but none with the expected sign /
               a path with the expected sign exists
  loop         the largest strongly connected component the region passes
               through, and whether that loop is
                 negative   it contains an inhibitory edge (negative feedback)
                 positive   every edge in it is activating
  derived      whether the loop survives without derived edges (assembly,
               depletion, catalyst); if not it is welded by our own expansion
  self_inh     a self-contained inhibitor inside the region: an inhibitor that
               contains an input of the same reaction (specs/012, 022)
  assembly     an assembly edge inside the region (a complex built from the
               perturbed branch; matters for overexpression)

and the error type: false_change (expected NORMAL), missed (predicted NORMAL),
wrong_direction.

The pinned and readout uuids come from a case dump that records them
(`gene_uuids`, `output_uuids`). They depend only on the pinning protocol, not
on the solver, so `--uuids` can supply them for an arm whose dump predates the
columns.

Usage:
  failure_structure.py --cases ARM/curator_cases.tsv --catalog BUILD \
      [--uuids OTHER/curator_cases.tsv] [--out annotated.tsv]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cycle_structure import tarjan_sccs  # noqa: E402
from holdout_report import TUNING_PATHWAYS  # noqa: E402
from benchmark_vs_mpbiopath import self_contained_inhibitor_pairs  # noqa: E402

DERIVED = {"assembly", "depletion", "catalyst"}


class Network:
    def __init__(self, pathway_dir: Path):
        self.fwd = defaultdict(list)          # u -> [(v, is_neg, edge_type)]
        self.bwd = defaultdict(list)
        nodes = set()
        with open(pathway_dir / "logic_network.csv", newline="") as fh:
            for r in csv.DictReader(fh):
                u, v = r["source_id"], r["target_id"]
                neg = r.get("pos_neg") == "neg"
                et = r.get("edge_type", "")
                self.fwd[u].append((v, neg, et))
                self.bwd[v].append(u)
                nodes |= {u, v}
        adj = {u: [v for v, _, _ in es] for u, es in self.fwd.items()}
        self.comp = {}
        self.comp_size = {}
        self.comp_neg = {}
        self.comp_derived_only = {}
        for i, c in enumerate(tarjan_sccs(adj, nodes)):
            if len(c) < 2:
                continue
            cs = set(c)
            for n in c:
                self.comp[n] = i
            self.comp_size[i] = len(c)
            inner = [(u, v, neg, et) for u in c for v, neg, et in self.fwd[u] if v in cs]
            self.comp_neg[i] = any(neg for _, _, neg, _ in inner)
            # Does the component still cycle once derived edges are removed?
            core = defaultdict(list)
            for u, v, _, et in inner:
                if et not in DERIVED:
                    core[u].append(v)
            self.comp_derived_only[i] = not any(
                len(cc) > 1 for cc in tarjan_sccs(core, cs))
        self.self_inh = self_contained_inhibitor_pairs(pathway_dir)

    def region(self, sources, targets):
        fwd = self._bfs(sources, lambda n: (v for v, _, _ in self.fwd.get(n, ())))
        bwd = self._bfs(targets, lambda n: self.bwd.get(n, ()))
        return fwd & bwd

    @staticmethod
    def _bfs(start, nbrs):
        seen = set(start)
        q = deque(start)
        while q:
            n = q.popleft()
            for m in nbrs(n):
                if m not in seen:
                    seen.add(m)
                    q.append(m)
        return seen

    def parities(self, sources, targets, region):
        """Sign parities (0 = even number of inhibitions) with which any target
        is reachable from any source, walking only inside the region."""
        tset = set(targets)
        seen = {(s, 0) for s in sources if s in region}
        q = deque(seen)
        found = set()
        while q:
            n, p = q.popleft()
            if n in tset:
                found.add(p)
                if len(found) == 2:
                    break
            for v, neg, _ in self.fwd.get(n, ()):
                if v in region:
                    st = (v, p ^ int(neg))
                    if st not in seen:
                        seen.add(st)
                        q.append(st)
        return found


def classify(net: Network, row: dict) -> dict:
    src = [u for u in row["gene_uuids"].split("|") if u]
    tgt = [u for u in row["output_uuids"].split("|") if u]
    reg = net.region(src, tgt)
    out = {}
    pred, exp, d = row["predicted"], row["expected"], row["direction"]
    if exp == "1":
        out["error"] = "false_change"
    elif pred == "1":
        out["error"] = "missed"
    else:
        out["error"] = "wrong_direction"
    if not reg:
        out["reach"] = "no_path"
    else:
        par = net.parities(src, tgt, reg)
        if exp == "1":
            out["reach"] = "path"
        else:
            # KD (0) expecting DOWN, or OE (2) expecting UP, needs even parity.
            need = 0 if exp == d else 1
            out["reach"] = "signed_path" if need in par else "wrong_sign_only"
    comps = {net.comp[n] for n in reg if n in net.comp}
    if comps:
        big = max(comps, key=lambda c: net.comp_size[c])
        size = net.comp_size[big]
        out["loop"] = ("negative" if any(net.comp_neg[c] for c in comps) else "positive")
        out["loop_size"] = "small(<10)" if size < 10 else "medium(<100)" if size < 100 else "giant(>=100)"
        out["derived_only"] = "welded" if all(net.comp_derived_only[c] for c in comps) else "curated"
    else:
        out["loop"] = out["loop_size"] = out["derived_only"] = "none"
    out["self_inh"] = "yes" if any(s in reg and t in reg for s, t in net.self_inh) else "no"
    out["assembly"] = "yes" if any(
        et == "assembly" and u in reg and v in reg
        for u in reg for v, _, et in net.fwd.get(u, ())) else "no"
    return out


FEATURES = ("reach", "loop", "loop_size", "derived_only", "self_inh", "assembly")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--uuids", type=Path, help="a dump with gene_uuids/output_uuids for the same protocol")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.cases, newline=""), delimiter="\t"))
    src_rows = rows if a.uuids is None else list(csv.DictReader(open(a.uuids, newline=""), delimiter="\t"))
    key = lambda r: (r["pathway"], r["gene"], r["direction"], r["key_output"])
    uu = {key(r): (r.get("gene_uuids"), r.get("output_uuids")) for r in src_rows}
    if not any(g for g, _ in uu.values()):
        sys.exit("no gene_uuids in the uuid source; pass --uuids with a dump that records them")

    names = {}
    for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv"):
        if line.startswith("#") or line.startswith("id\t"):
            continue
        pid, name = line.rstrip("\n").split("\t")
        names[name] = pid

    nets: dict = {}
    scored = []
    for r in rows:
        if r["valid"] != "1":
            continue
        g, o = uu.get(key(r), (None, None))
        if not g or not o:
            continue
        r = dict(r, gene_uuids=g, output_uuids=o)
        pw = r["pathway"]
        if pw not in nets:
            nets[pw] = Network(a.catalog / names[pw])
        r["correct"] = r["predicted"] == r["expected"]
        r.update(classify(nets[pw], r) if not r["correct"] else
                 {k: v for k, v in classify(nets[pw], r).items() if k != "error"})
        scored.append(r)

    held = [r for r in scored if r["pathway"] not in TUNING_PATHWAYS]
    for label, sub in (("HELD-OUT", held), ("ALL", scored)):
        wrong = [r for r in sub if not r["correct"]]
        print(f"\n{label}: {len(sub):,} scored, {len(wrong):,} wrong ({len(wrong) / max(1, len(sub)):.1%})")
        print("  error types:", dict(Counter(r["error"] for r in wrong).most_common()))
        for f in FEATURES:
            c_all = Counter(r[f] for r in sub)
            c_bad = Counter(r[f] for r in wrong)
            cells = "  ".join(f"{v}: {c_bad[v]:,}/{c_all[v]:,} wrong ({c_bad[v] / c_all[v]:.1%})"
                              for v, _ in c_all.most_common())
            print(f"  {f:<13} {cells}")
        print("  wrong, by error type x reach:")
        for (e, rc), n in Counter((r["error"], r["reach"]) for r in wrong).most_common(8):
            print(f"    {e:<16} {rc:<16} {n:,}")

    if a.out:
        cols = list(rows[0].keys()) + [c for c in ("correct", "error") + FEATURES if c not in rows[0]]
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(scored)
        print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
