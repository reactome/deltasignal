#!/usr/bin/env python3
"""specs/052: score the zero-supply arm against its control, as registered.

    zero_supply_report.py BUILD CONTROL_DIR ARM_DIR

Per axis (curator tuning, curator held-out, experimental): cases fixed and
broken, net, exact McNemar, macro-F1 for both arms, pathways moved, distinct
readouts among discordant cases, and the knockout->UP subgroup.

**By-construction cases are reported apart.** A knockout case whose readout
nodes are ALL held at 0 by the rule reads DOWN because of the rule, not the
propagation. `held_by_rule` mirrors src/core/zero_supply.jl on the case's
pins (the dumped `gene_uuids`) so each case can be classified.
"""
import collections
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from holdout_report import TUNING_PATHWAYS, load, macro_f1  # noqa: E402
from network_files import open_network  # noqa: E402

CARRYING = {"input", "output", "assembly", "dissociation", "set_member", "variant_pool", "variant_split"}
REQUIRING = {"simple_entity", "simple_complex", "set_variant", "dissociation_sink"}


def held_by_rule(nodes, edges, pinned):
    """Python mirror of zero_supply_uuids: `nodes` uuid -> (kind, leaves set),
    `edges` [(src, dst, type)], `pinned` uuid -> UI value."""
    held = set()
    cands = set().union(*(nodes[u][1] for u, v in pinned.items() if v == 0.0 and u in nodes)) \
        if any(v == 0.0 for v in pinned.values()) else set()
    out = collections.defaultdict(list)
    for s, d, t in edges:
        if t in CARRYING:
            out[s].append((d, t))
    kind = lambda u: nodes[u][0] if u in nodes else ""  # noqa: E731
    for L in sorted(cands):
        carrier = {u for u, (k, ls) in nodes.items() if k != "reaction" and L in ls}
        for u in list(carrier):
            for d, t in out.get(u, ()):
                if t == "input" and kind(d) == "reaction":
                    carrier.add(d)
        has_in = {d for u in carrier for d, _ in out.get(u, ()) if d in carrier}
        zero = lambda u: pinned.get(u, -1.0) == 0.0  # noqa: E731
        frontier = [u for u in carrier if not zero(u) and
                    ((kind(u) != "reaction" and u not in has_in) or pinned.get(u, 0.0) > 0.0)]
        supplied = set(frontier)
        while frontier:
            u = frontier.pop()
            for d, _ in out.get(u, ()):
                if d in carrier and d not in supplied and not zero(d):
                    supplied.add(d)
                    frontier.append(d)
        held |= {u for u in carrier if u not in supplied and u not in pinned and kind(u) in REQUIRING}
    return held


def load_network(pdir):
    nodes = {}
    with open(pdir / "nodes.csv", newline="") as f:
        for r in csv.DictReader(f):
            nodes[r["uuid"]] = (r["node_kind"], set(filter(None, (r.get("member_leaves") or "").split("|"))))
    with open_network(pdir) as f:
        edges = [(r["source_id"], r["target_id"], r["edge_type"]) for r in csv.DictReader(f)]
    return nodes, edges


def mcnemar(b, c):
    n = b + c
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n)


def raw_rows(path):
    with open(path, newline="") as f:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(f, delimiter="\t") if r["expected"] in ("0", "1", "2")}


def by_construction(build, rows):
    """Keys of knockout cases whose every readout node the rule holds."""
    dirs = {p.name.split("_R-HSA")[0] if "_R-HSA" in p.name else p.name: p
            for p in build.iterdir() if p.is_dir() and (p / "nodes.csv").exists()}
    names = {}
    pt = build / "pathways.tsv"
    if pt.exists():
        for r in csv.reader(open(pt), delimiter="\t"):
            if len(r) >= 2:
                names[r[1]] = r[0]
    cache, out = {}, set()
    for k, r in rows.items():
        if r["direction"] != "0" or not r.get("output_uuids"):
            continue
        pid = names.get(r["pathway"])
        pdir = build / pid if pid and (build / pid).exists() else dirs.get(r["pathway"])
        if pdir is None:
            continue
        if pdir not in cache:
            cache[pdir] = load_network(pdir)
        nodes, edges = cache[pdir]
        pins = {u: 0.0 for u in r["gene_uuids"].split("|") if u}
        outs = [u for u in r["output_uuids"].split("|") if u]
        if pins and outs and set(outs) <= held_by_rule(nodes, edges, pins):
            out.add(k)
    return out


def compare(ctrl, arm, keys, bc):
    fixed = [k for k in keys if arm[k]["predicted"] == arm[k]["expected"] != ctrl[k]["predicted"]]
    broken = [k for k in keys if ctrl[k]["predicted"] == ctrl[k]["expected"] != arm[k]["predicted"]]
    ko_up = [k for k in keys if k[2] == "0" and ctrl[k]["predicted"] == "2"]
    disc = fixed + broken
    return {
        "n": len(keys),
        "control_correct": sum(ctrl[k]["predicted"] == ctrl[k]["expected"] for k in keys),
        "arm_correct": sum(arm[k]["predicted"] == arm[k]["expected"] for k in keys),
        "fixed": len(fixed), "broken": len(broken), "net": len(fixed) - len(broken),
        "mcnemar_p": mcnemar(len(fixed), len(broken)),
        "macro_f1": (macro_f1([ctrl[k] for k in keys]), macro_f1([arm[k] for k in keys])),
        "pathways_moved": len({k[0] for k in disc}),
        "distinct_readouts": len({(k[0], k[3]) for k in disc}),
        "by_construction": {"fixed": sum(k in bc for k in fixed), "broken": sum(k in bc for k in broken)},
        "net_without_by_construction": sum(k not in bc for k in fixed) - sum(k not in bc for k in broken),
        "ko_up_control": len(ko_up),
        "ko_up_wrong_fixed": sum(1 for k in ko_up if k in fixed),
        "ko_up_right_broken": sum(1 for k in ko_up if k in broken),
        "by_pathway": dict(collections.Counter(k[0] for k in fixed).most_common()),
        "by_pathway_broken": dict(collections.Counter(k[0] for k in broken).most_common()),
    }


def main(argv):
    build, ctrl_dir, arm_dir = map(Path, argv[1:4])
    res = {}
    for axis in ("curator", "experimental"):
        ctrl = raw_rows(ctrl_dir / f"{axis}_cases.tsv")
        arm = raw_rows(arm_dir / f"{axis}_cases.tsv")
        keys = sorted(set(ctrl) & set(arm))
        bc = by_construction(build, {k: arm[k] for k in keys})
        slices = ({"tuning": [k for k in keys if k[0] in TUNING_PATHWAYS],
                   "held-out": [k for k in keys if k[0] not in TUNING_PATHWAYS]}
                  if axis == "curator" else {"all": keys})
        res[axis] = {"unpaired": (len(ctrl), len(arm), len(keys))}
        for sl, ks in slices.items():
            m = compare(ctrl, arm, ks, bc)
            res[axis][sl] = m
            print(f"{axis:12s} {sl:8s} n={m['n']:6d} correct {m['control_correct']} -> {m['arm_correct']} "
                  f"(+{m['fixed']} / -{m['broken']}, net {m['net']:+d}, p {m['mcnemar_p']:.3g}); "
                  f"mF1 {m['macro_f1'][0]:.4f} -> {m['macro_f1'][1]:.4f}; pathways {m['pathways_moved']}, "
                  f"readouts {m['distinct_readouts']}; by construction +{m['by_construction']['fixed']} "
                  f"/ -{m['by_construction']['broken']} (net without {m['net_without_by_construction']:+d}); "
                  f"KO->UP {m['ko_up_control']}: wrong fixed {m['ko_up_wrong_fixed']}, "
                  f"right broken {m['ko_up_right_broken']}")
            print(f"    fixed by pathway  {m['by_pathway']}")
            print(f"    broken by pathway {m['by_pathway_broken']}")
    return res


if __name__ == "__main__":
    main(sys.argv)
