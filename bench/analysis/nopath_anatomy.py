#!/usr/bin/env python3
"""Why is a `no_path` case unreachable? Four possibilities, told apart.

`no_path` says the logic network had no route from the perturbed gene to the
readout. It does not say why, and the answers need completely different work:

  1. UUID SILO      a route exists once positional variants of the same
                    stable id are collapsed. The biology connects; our node
                    identity split it. Fixable in the generator/solver.
  2. MPB CONNECTS   MP-BioPath's own network connects it and ours does not.
                    They had an edge we lack. Fixable upstream.
  3. REACTOME ONLY  neither logic network connects it, but Reactome does.
  4. NEITHER        nothing connects it; the curator's call rests on something
                    no directed graph encodes (co-regulation, feedback,
                    knowledge outside the pathway).

Only 4 is genuinely not ours to fix, and the split says how much is.
"""
from __future__ import annotations
import argparse, collections, csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from benchmark_vs_mpbiopath import (neo4j_gene_to_stids, neo4j_dbid_to_stid,  # noqa: E402
                                    load_stid_to_uuids)


def load_edges(pathway_dir: Path):
    f = pathway_dir / "logic_network.csv"
    if not f.exists():
        return []
    with f.open(newline="") as fh:
        return [(r["source_id"], r["target_id"]) for r in csv.DictReader(fh)]


def reaction_stids(pathway_dir: Path) -> set[str]:
    """Stable ids that name a REACTION, not an entity.

    These appear as nodes AND as edge_reaction_id -- 156 of 841 in WNT. They
    must NOT be collapsed by stable id: one curated reaction has many positional
    instances, and merging them manufactures routes between unrelated branches.
    Collapsing them inflated the silo figure on the first pass.
    """
    f = pathway_dir / "logic_network.csv"
    if not f.exists():
        return set()
    with f.open(newline="") as fh:
        return {r["edge_reaction_id"] for r in csv.DictReader(fh) if r.get("edge_reaction_id")}


def uuid_to_stid(pathway_dir: Path) -> dict[str, str]:
    f = pathway_dir / "stid_to_uuid_mapping.csv"
    if not f.exists():
        return {}
    with f.open(newline="") as fh:
        return {r["uuid"]: r["stable_id"] for r in csv.DictReader(fh)}


def adjacency(edges, collapse: dict[str, str] | None):
    adj = collections.defaultdict(set)
    for s, t in edges:
        if collapse is not None:
            s, t = collapse.get(s, s), collapse.get(t, t)
        adj[s].add(t)
    return adj


def hops(adj, sources: set[str], targets: set[str]) -> int | None:
    """Shortest number of hops from any source to any target, or None."""
    if not sources or not targets:
        return None
    if sources & targets:
        return 0
    seen = set(sources)
    frontier = list(sources)
    d = 0
    while frontier:
        d += 1
        nxt = []
        for node in frontier:
            for n in adj.get(node, ()):
                if n in targets:
                    return d
                if n not in seen:
                    seen.add(n); nxt.append(n)
        frontier = nxt
    return None


def reaches(adj, sources: set[str], targets: set[str]) -> bool:
    if not sources or not targets:
        return False
    if sources & targets:
        return True
    seen, q = set(sources), collections.deque(sources)
    while q:
        for nxt in adj.get(q.popleft(), ()):
            if nxt in targets:
                return True
            if nxt not in seen:
                seen.add(nxt); q.append(nxt)
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--ours", type=Path, required=True)
    ap.add_argument("--mpb", type=Path, required=True)
    ap.add_argument("--category", default="no_path")
    ap.add_argument("--pathway-list", type=Path,
                    default=Path.home() / "gitroot/mp-biopath-pathways/pathway_list.tsv",
                    help="Maps curator pathway NAME to id. Needed because the "
                         "curator files and the catalog dirs disagree on names "
                         "(a comma, a trailing underscore, a rename) -- matching "
                         "on name silently drops whole pathways.")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--separate", action="store_true",
                    help="Can the cases collapsing WOULD FIX be told apart from "
                         "the ones it would break? Reports the collapsed-path "
                         "length for each group, raw and controlled per pathway "
                         "-- a raw difference is usually a between-pathway "
                         "confound (it was for readout in-degree).")
    ap.add_argument("--exposure", action="store_true",
                    help="Score EVERY case, not just failures, and report what "
                         "collapsing variants would gain against what it puts "
                         "at risk. Four silo-bridge attempts have measured "
                         "negative; the ratio is why.")
    a = ap.parse_args()

    allrows = list(csv.DictReader(a.cases.open(newline=""), delimiter="\t"))
    if a.exposure or a.separate:
        # BOTH groups are needed: the RISK cases are currently CORRECT, so the
        # failures-only filter would silently leave the comparison one-sided.
        rows = [r for r in allrows if r["expected"] in ("0", "1", "2")]
        print(f"{len(rows)} scored cases, measuring silo exposure", flush=True)
    else:
        rows = [r for r in allrows
                if r["category"] == a.category and r["expected"] != r["predicted"]]
        print(f"{len(rows)} {a.category} failures to explain", flush=True)

    # Resolve by pathway ID, the way the benchmark does. Name matching dropped
    # 438 of 1,484 cases into a bogus "no network" bucket on the first run.
    name_to_id = {}
    with a.pathway_list.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            name_to_id[r["pathway_name"]] = r["pathway_id"]

    def index(root: Path):
        by_id = {}
        for d in root.iterdir():
            if d.is_dir() and "_R-HSA-" in d.name:
                by_id[d.name.rsplit("_R-HSA-", 1)[1]] = d
        return by_id

    ours_by_id, mpb_by_id = index(a.ours), index(a.mpb)

    def resolve(by_id, name):
        pid = name_to_id.get(name) or name_to_id.get(name.rstrip("_"))
        if pid is None:
            for k, v in name_to_id.items():
                if k.replace(",", "") == name.replace(",", "").rstrip("_"):
                    pid = v
                    break
        return by_id.get(pid) if pid else None

    ours_dirs = type("D", (), {"get": lambda self, n: resolve(ours_by_id, n)})()
    mpb_dirs = type("D", (), {"get": lambda self, n: resolve(mpb_by_id, n)})()

    genes = sorted({r["gene"] for r in rows})
    gene_stids = neo4j_gene_to_stids(genes)
    ko_stid = neo4j_dbid_to_stid({r["key_output"] for r in rows})

    cache: dict[str, tuple] = {}
    verdicts = collections.Counter()
    detail = []

    for r in rows:
        pw = r["pathway"]
        if pw not in cache:
            od = ours_dirs.get(pw)
            if od is None:
                cache[pw] = None
            else:
                u2s_all = uuid_to_stid(od)
                rxn = reaction_stids(od)
                # Entities collapse by stable id; reaction instances stay distinct.
                u2s = {u: sid for u, sid in u2s_all.items() if sid not in rxn}
                oe = load_edges(od)
                md = mpb_dirs.get(pw)
                cache[pw] = (
                    load_stid_to_uuids(od),          # stid -> [uuid] (ours)
                    adjacency(oe, None),             # uuid-level, ours
                    adjacency(oe, u2s),              # stid-level, ours
                    adjacency(load_edges(md), uuid_to_stid(md)) if md else None,
                    set(u2s_all.values()),
                )
        if cache[pw] is None:
            verdicts["no_network"] += 1
            continue
        s2u, adj_uuid, adj_stid, adj_mpb, present = cache[pw]

        g_stids = {s for g in [r["gene"]] for s in gene_stids.get(g, [])}
        k_stid = ko_stid.get(r["key_output"]) or f"R-HSA-{r['key_output']}"
        g_uuids = {u for s in g_stids for u in s2u.get(s, [])}
        k_uuids = set(s2u.get(k_stid, []))

        if not g_uuids or not k_uuids:
            verdicts["endpoint missing from our network"] += 1
            continue

        by_uuid = reaches(adj_uuid, g_uuids, k_uuids)
        by_stid = reaches(adj_stid, g_stids, {k_stid})
        by_mpb = reaches(adj_mpb, g_stids, {k_stid}) if adj_mpb is not None else None

        if a.separate:
            if by_uuid or not by_stid:
                continue
            d = hops(adj_stid, g_stids, {k_stid})
            grp = ("GAIN" if r["expected"] != "1" and r["predicted"] == "1"
                   else "RISK" if r["expected"] == "1" and r["predicted"] == "1"
                   else None)
            if grp and d is not None:
                detail.append({"pathway": pw, "group": grp, "hops": d})
            continue
        if a.exposure:
            # Only cases the silo currently BLOCKS can be changed by collapsing.
            if by_uuid or not by_stid:
                verdicts["unaffected by collapsing"] += 1
                continue
            if r["expected"] == "1":
                v = ("AT RISK: NORMAL, currently correct" if r["predicted"] == "1"
                     else "AT RISK: NORMAL, already wrong")
            else:
                v = ("GAIN: real change we currently miss" if r["predicted"] == "1"
                     else "change case, already predicted as a change")
            verdicts[v] += 1
            continue
        if by_uuid:
            v = "REACHABLE at uuid level (no_path label suspect)"
        elif by_stid:
            v = "1. UUID SILO — collapsing variants connects it"
        elif by_mpb:
            v = "2. MPB CONNECTS, we do not"
        elif adj_mpb is None:
            v = "no MP-BioPath network for this pathway"
        else:
            v = "4. NEITHER network connects it"
        verdicts[v] += 1
        detail.append({**{k: r[k] for k in ("pathway", "gene", "key_output", "expected")},
                       "verdict": v})

    if a.separate:
        import statistics
        by = collections.defaultdict(list)
        for d in detail:
            by[d["group"]].append(d["hops"])
        print()
        for g in ("GAIN", "RISK"):
            v = sorted(by[g])
            print(f"{g:<5} n={len(v):<6} median hops {statistics.median(v):>4.1f}  "
                  f"quartiles {v[len(v)//4]}/{v[3*len(v)//4]}  max {max(v)}")
        print("\nWITHIN PATHWAY (the control that killed readout in-degree):")
        per = collections.defaultdict(lambda: collections.defaultdict(list))
        for d in detail:
            per[d["pathway"]][d["group"]].append(d["hops"])
        diffs = []
        for pw, g in per.items():
            if len(g["GAIN"]) >= 5 and len(g["RISK"]) >= 5:
                diffs.append(statistics.median(g["GAIN"]) - statistics.median(g["RISK"]))
        if diffs:
            diffs.sort()
            m = diffs[len(diffs)//2]
            print(f"  {len(diffs)} pathways with >=5 of each. median difference "
                  f"in hops (GAIN - RISK): {m:+.1f}")
            print(f"  GAIN closer in {sum(1 for d in diffs if d < 0)}, "
                  f"further in {sum(1 for d in diffs if d > 0)}, "
                  f"equal in {sum(1 for d in diffs if d == 0)}")
        else:
            print("  no pathway has >=5 of each -- cannot control")
        # Would a hop cap actually help?
        print("\nIf we only reconnected paths within N hops:")
        for cap in (1, 2, 3, 4, 6, 8):
            g = sum(1 for d in detail if d["group"] == "GAIN" and d["hops"] <= cap)
            rk = sum(1 for d in detail if d["group"] == "RISK" and d["hops"] <= cap)
            print(f"  <= {cap:>2} hops: gain {g:>5}, risk {rk:>5}, ratio 1:{rk/g:.1f}" if g
                  else f"  <= {cap:>2} hops: gain     0")
        return 0

    print()
    total = sum(verdicts.values())
    for v, n in verdicts.most_common():
        print(f"{n:>6}  ({n/total:5.1%})  {v}")
    if a.out:
        with a.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(detail[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(detail)
        print(f"\nper-case detail: {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
