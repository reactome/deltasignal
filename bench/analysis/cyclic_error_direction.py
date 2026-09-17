#!/usr/bin/env python3
"""Split benchmark errors by direction, separately for cyclic and acyclic readouts.

specs/008 found that AND-joined cycles collapse and cannot be rescued by any
external input. If that mechanism drives the 17pp cyclic-readout accuracy gap,
cyclic readouts must show an EXCESS of false DOWN calls relative to acyclic
ones. This measures that, read-only, against an existing case dump — no solving.
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from cycle_structure import read_lng_network, tarjan_sccs  # noqa: E402
from benchmark_vs_mpbiopath import load_stid_to_uuids, neo4j_dbid_to_stid, \
    neo4j_set_member_stids  # noqa: E402

LAB = {"0": "DOWN", "1": "NORMAL", "2": "UP"}


def reach_index(pathway_dir: Path):
    """Reverse adjacency + node count, for the fraction of the pathway that can
    reach a given readout. A readout inside a giant SCC is reachable from most
    of the network; that is a confound for cyclicity, not a consequence of it."""
    edges = read_lng_network(pathway_dir)
    rev: dict[str, list[str]] = defaultdict(list)
    nodes: set[str] = set()
    for e in edges:
        rev[e["target"]].append(e["source"])
        nodes.add(e["source"]); nodes.add(e["target"])
    return rev, max(len(nodes), 1)


def ancestors(rev, seeds) -> int:
    seen, stack = set(seeds), list(seeds)
    while stack:
        for p in rev.get(stack.pop(), ()):
            if p not in seen:
                seen.add(p); stack.append(p)
    return len(seen)


def scc_members(pathway_dir: Path, min_size: int) -> set[str]:
    edges = read_lng_network(pathway_dir)
    adj: dict[str, list[str]] = defaultdict(list)
    nodes: set[str] = set()
    for e in edges:
        adj[e["source"]].append(e["target"])
        nodes.add(e["source"]); nodes.add(e["target"])
    out: set[str] = set()
    for comp in tarjan_sccs(adj, nodes):
        if len(comp) >= min_size:
            out.update(comp)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--pathway-list", type=Path,
                    default=Path.home() / "gitroot/mp-biopath-pathways/pathway_list.tsv",
                    help="Maps the curator pathway NAME recorded in a case dump "
                         "to its id, so directories are found by id.")
    ap.add_argument("--min-scc", type=int, default=2,
                    help="Smallest SCC counted as cyclic (default 2)")
    a = ap.parse_args()

    rows = list(csv.DictReader(a.cases.open(newline=""), delimiter="\t"))
    # Keyed by ID, not by name -- see pathway_dir_index().
    from _common import pathway_dir_index, name_to_id_map
    by_id = pathway_dir_index(a.catalog)
    name_to_id = name_to_id_map(a.pathway_list)
    dirs = type("D", (), {"get": staticmethod(
        lambda n: by_id.get(name_to_id.get(n) or ""))})()

    # Resolve every key_output once, via the same path the benchmark uses.
    kos = {r["key_output"] for r in rows if r.get("key_output")}
    dbid_to_stid = neo4j_dbid_to_stid(kos)

    cache: dict[str, tuple[set[str], dict[str, list[str]]]] = {}
    tally: dict[tuple[str, str, str], int] = defaultdict(int)
    matched: dict[tuple[int, str], list[int]] = defaultdict(lambda: [0, 0])
    per_pw: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    unresolved = 0

    for r in rows:
        pw = r["pathway"]
        if pw not in cache:
            d = dirs.get(pw)
            cache[pw] = ((scc_members(d, a.min_scc), load_stid_to_uuids(d))
                         + reach_index(d)) if d else (set(), {}, {}, 1)
        members, stid_to_uuids, rev, n_nodes = cache[pw]

        ko = r.get("key_output") or ""
        cands = []
        if dbid_to_stid.get(ko):
            cands.append(dbid_to_stid[ko])
        if f"R-HSA-{ko}" not in cands:
            cands.append(f"R-HSA-{ko}")
        uuids: list[str] = []
        for sid in cands:
            if stid_to_uuids.get(sid):
                uuids = stid_to_uuids[sid]
                break
        if not uuids:
            for sid, ms in neo4j_set_member_stids(cands).items():
                for m in ms:
                    uuids += stid_to_uuids.get(m, [])
        if not uuids:
            unresolved += 1
            continue

        kind = "cyclic" if any(u in members for u in uuids) else "acyclic"
        exp, pred = LAB.get(r["expected"], "?"), LAB.get(r["predicted"], "?")
        tally[(kind, exp, pred)] += 1
        if exp == "NORMAL":
            reach = ancestors(rev, uuids) / n_nodes
            bucket = min(int(reach * 5), 4)     # quintiles of pathway reach
            matched[(bucket, kind)][0] += 1
            per_pw[(pw, kind)][0] += 1
            if pred != "NORMAL":
                matched[(bucket, kind)][1] += 1
                per_pw[(pw, kind)][1] += 1

    print(f"readout uuids unresolved (skipped): {unresolved}\n")
    labs = ["DOWN", "NORMAL", "UP"]
    for kind in ("cyclic", "acyclic"):
        n = sum(v for (k, _, _), v in tally.items() if k == kind)
        ok = sum(v for (k, e, p), v in tally.items() if k == kind and e == p)
        if not n:
            continue
        print(f"=== {kind} readouts: {n} cases, accuracy {ok/n:.4f} ===")
        print(f"{'exp\\pred':<12}" + "".join(f"{l:>9}" for l in labs))
        for e in labs:
            print(f"{e:<12}" + "".join(f"{tally[(kind,e,p)]:>9}" for p in labs))
        err = n - ok
        fd = sum(v for (k, e, p), v in tally.items()
                 if k == kind and p == "DOWN" and e != "DOWN")
        fu = sum(v for (k, e, p), v in tally.items()
                 if k == kind and p == "UP" and e != "UP")
        fn = err - fd - fu
        if err:
            print(f"errors {err}: false DOWN {fd} ({fd/err:.1%})  "
                  f"false UP {fu} ({fu/err:.1%})  false NORMAL {fn} ({fn/err:.1%})\n")
    print("=== false-change rate on truly-NORMAL cases, by readout reach ===")
    print("(reach = fraction of the pathway that can reach the readout)")
    print(f"{'reach':<12}{'cyclic n':>10}{'cyclic FC':>11}"
          f"{'acyclic n':>11}{'acyclic FC':>12}")
    for b in range(5):
        cn, cf = matched[(b, "cyclic")]
        an, af = matched[(b, "acyclic")]
        cs = f"{cf/cn:.3f}" if cn >= 20 else "-"
        as_ = f"{af/an:.3f}" if an >= 20 else "-"
        print(f"{b*20:>3}-{b*20+20:<8}{cn:>10}{cs:>11}{an:>11}{as_:>12}")
    # WITHIN-PATHWAY control. A raw cyclic-vs-acyclic difference can be
    # "cyclic readouts live in harder pathways" -- the same between-pathway
    # confound that made readout in-degree look like a 4x effect.
    print("\n=== within-pathway control, truly-NORMAL cases ===")
    print("pathways with >=20 cyclic AND >=20 acyclic NORMAL cases")
    print(f"{'pathway':<52}{'cyc FC':>8}{'acyc FC':>9}{'diff':>8}")
    diffs = []
    for pw in sorted({k[0] for k in per_pw}):
        cn, cf = per_pw[(pw, "cyclic")]
        an, af = per_pw[(pw, "acyclic")]
        if cn < 20 or an < 20:
            continue
        d = cf / cn - af / an
        diffs.append(d)
        print(f"{pw[:50]:<52}{cf/cn:>8.3f}{af/an:>9.3f}{d:>+8.3f}")
    if diffs:
        diffs.sort()
        med = diffs[len(diffs) // 2] if len(diffs) % 2 else \
            (diffs[len(diffs)//2 - 1] + diffs[len(diffs)//2]) / 2
        print(f"\n{len(diffs)} pathways compared. median difference {med:+.3f}; "
              f"cyclic worse in {sum(1 for d in diffs if d > 0)}, "
              f"better in {sum(1 for d in diffs if d < 0)}")
        # COVERAGE. The control needs >=20 of BOTH kinds in a pathway, which
        # excludes pathways where cyclic readouts dominate -- exactly where a
        # real loop effect would hide. State what fraction of the cyclic
        # population the control actually saw, so the conclusion is not read
        # more widely than the evidence.
        compared = {pw for pw in {k[0] for k in per_pw}
                    if per_pw[(pw, "cyclic")][0] >= 20
                    and per_pw[(pw, "acyclic")][0] >= 20}
        seen = sum(per_pw[(pw, "cyclic")][0] for pw in compared)
        allc = sum(v[0] for k, v in per_pw.items() if k[1] == "cyclic")
        pw_with_cyclic = {k[0] for k, v in per_pw.items()
                          if k[1] == "cyclic" and v[0] > 0}
        print(f"COVERAGE: the control saw {seen} of {allc} cyclic NORMAL cases "
              f"({seen/allc:.1%}), in {len(compared)} of {len(pw_with_cyclic)} "
              f"pathways that have any.")
        excl = sorted(((per_pw[(pw,'cyclic')][0], pw) for pw in pw_with_cyclic
                       if pw not in compared), reverse=True)[:5]
        if excl:
            print("largest EXCLUDED pathways (cyclic NORMAL n, acyclic n):")
            for n, pw in excl:
                print(f"   {n:>5}  {per_pw[(pw,'acyclic')][0]:>5}  {pw[:52]}")
    else:
        print("\nNO pathway has enough of both -- cyclic and acyclic readouts do "
              "not co-occur, so the raw comparison is BETWEEN pathways only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
