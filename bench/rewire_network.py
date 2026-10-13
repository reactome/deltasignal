"""specs/051: the rewired-network null, as its own catalog build.

Copies a complete build and rewires each pathway's reaction-level edges by
degree-preserving double-edge swaps inside an edge class (edge_type, pos_neg,
and_or). Entity-structure edges (boundary_edges.csv, variant_pool,
variant_split) are left alone. The pool and carrier-loop tables name specific
steps of the real wiring, so they are deleted, not emptied: a run that does
not switch those rules off refuses instead of solving with zero pools.

The design, the classes and the decision rule were registered in
specs/051-dataset-roles/spec.md before this script was written.

Usage: rewire_network.py SOURCE_BUILD_ID SEED
       -> ~/deltasignal-catalogs/builds/<SOURCE_BUILD_ID>_rewire<SEED>/
"""
import csv
import io
import json
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

REWIRED_TYPES = ("input", "output", "catalyst", "regulator", "depletion")
DROPPED_TABLES = ("pools.csv", "pool_transitions.csv", "pool_carriers.csv", "carrier_loops.csv")
NOT_COPIED = ("cache",)            # generator internals that describe the real wiring
SWAPS_PER_EDGE = 10
TRIES_PER_EDGE = 100
BUILDS = Path.home() / "deltasignal-catalogs" / "builds"


def edge_class(row):
    return (row["edge_type"], row["pos_neg"], row["and_or"])


def swap_class(pairs, rng, swaps_per_edge=SWAPS_PER_EDGE, tries_per_edge=TRIES_PER_EDGE):
    """Double-edge swaps over a list of (source, target) pairs, in place.

    (a->r1), (b->r2) becomes (a->r2), (b->r1). Every source keeps its
    out-degree and every target its in-degree. A swap is rejected if it makes
    a self-loop or repeats a pair already present. Returns (accepted, tried)."""
    n = len(pairs)
    if n < 2:
        return 0, 0
    present = Counter(pairs)
    want, cap = swaps_per_edge * n, tries_per_edge * n
    accepted = tried = 0
    while accepted < want and tried < cap:
        tried += 1
        i, j = rng.randrange(n), rng.randrange(n)
        (a, r1), (b, r2) = pairs[i], pairs[j]
        if a == b or r1 == r2 or a == r2 or b == r1:
            continue
        if present[(a, r2)] or present[(b, r1)]:
            continue
        present[(a, r1)] -= 1
        present[(b, r2)] -= 1
        present[(a, r2)] += 1
        present[(b, r1)] += 1
        pairs[i], pairs[j] = (a, r2), (b, r1)
        accepted += 1
    return accepted, tried


def rewire_rows(rows, seed, pathway):
    """Rewire the REWIRED_TYPES classes of one network. Returns (rows, stats).
    Row order and every column but target_id (and the then-stale
    edge_reaction_id) are kept."""
    out = [dict(r) for r in rows]
    by_class = {}
    for k, r in enumerate(out):
        if r["edge_type"] in REWIRED_TYPES:
            by_class.setdefault(edge_class(r), []).append(k)
    stats = {}
    for cls in sorted(by_class):
        idx = sorted(by_class[cls], key=lambda k: (out[k]["source_id"], out[k]["target_id"], k))
        pairs = [(out[k]["source_id"], out[k]["target_id"]) for k in idx]
        before = list(pairs)
        rng = random.Random(f"{seed}:{pathway}:{'|'.join(cls)}")
        accepted, tried = swap_class(pairs, rng)
        moved = 0
        for k, old, (_, tgt) in zip(idx, before, pairs):
            if tgt != old[1]:
                moved += 1
                out[k]["target_id"] = tgt
                if "edge_reaction_id" in out[k]:
                    out[k]["edge_reaction_id"] = ""
        stats["|".join(cls)] = {"edges": len(pairs), "swaps": accepted, "tried": tried, "moved": moved}
    return out, stats


def degrees(rows):
    """(out, in) degree per (node, class): what the rewiring must preserve."""
    deg = Counter()
    for r in rows:
        c = edge_class(r)
        deg[("out", r["source_id"], c)] += 1
        deg[("in", r["target_id"], c)] += 1
    return deg


def read_network(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return reader.fieldnames, list(reader)


def write_network(path, fields, rows):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    Path(path).write_text(buf.getvalue())


def rewire_pathway(src, dst, seed):
    dst.mkdir()
    for p in src.iterdir():
        if p.name in DROPPED_TABLES or p.name in NOT_COPIED:
            continue
        if p.is_dir():
            shutil.copytree(p, dst / p.name)
        else:
            shutil.copy2(p, dst / p.name)
    fields, rows = read_network(src / "logic_network.csv")
    new, stats = rewire_rows(rows, seed, src.name)
    if degrees(new) != degrees(rows):
        raise SystemExit(f"{src.name}: degrees changed; refusing to write a broken null")
    write_network(dst / "logic_network.csv", fields, new)
    return stats


def git_head(repo):
    sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", "bench"],
                           capture_output=True, text=True, check=True).stdout.strip()
    return sha, bool(dirty)


def main(argv):
    if len(argv) != 3 or not argv[2].isdigit():
        raise SystemExit(__doc__)
    source_id, seed = argv[1], int(argv[2])
    src = BUILDS / source_id
    meta = json.loads((src / "BUILD.json").read_text())
    if meta.get("status") != "complete":
        raise SystemExit(f"{source_id} is not a complete build")
    repo = Path(__file__).resolve().parent.parent
    sha, dirty = git_head(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit("bench/ is uncommitted; a null build must name the commit that made it")
    dst_id = f"{source_id}_rewire{seed}"
    dst = BUILDS / dst_id
    if dst.exists():
        raise SystemExit(f"{dst} exists; builds are never overwritten")
    dst.mkdir()
    started = time.time()
    per_pathway = {}
    for p in sorted(src.iterdir()):
        if p.is_dir() and (p / "logic_network.csv").exists():
            per_pathway[p.name] = rewire_pathway(p, dst / p.name, seed)
        elif p.is_file() and p.name not in ("BUILD.json",):
            shutil.copy2(p, dst / p.name)    # pathways.tsv etc.: decryptm_coverage reads it
    # A class that is (nearly) complete bipartite cannot be randomised; the null
    # then equals the real wiring there, which makes it harder, not easier, to beat.
    stuck = sorted((pw, cls, st["edges"], st["moved"]) for pw, stats in per_pathway.items()
                   for cls, st in stats.items() if st["edges"] >= 20 and st["moved"] < st["edges"] / 2)
    totals = Counter()
    for stats in per_pathway.values():
        for s in stats.values():
            totals.update(s)
    out = dict(meta)
    out.update({
        "build_id": dst_id,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "variant": f"rewire{seed}",
        "rewire": {
            "spec": "specs/051-dataset-roles (rewired-network null)",
            "source_build": source_id,
            "seed": seed,
            "classes": list(REWIRED_TYPES),
            "swaps_per_edge": SWAPS_PER_EDGE,
            "tries_per_edge": TRIES_PER_EDGE,
            "tables_deleted": list(DROPPED_TABLES),
            "script_commit": sha + ("-dirty" if dirty else ""),
            "edges": totals["edges"], "swaps": totals["swaps"], "moved": totals["moved"],
            "classes_mostly_unmoved": [{"pathway": pw, "class": c, "edges": e, "moved": m} for pw, c, e, m in stuck],
            "per_pathway": per_pathway,
            "seconds": round(time.time() - started, 1),
        },
        "status": "complete",
    })
    # status goes last, so a crash before this line leaves no BUILD.json and
    # run_arm refuses the directory
    (dst / "BUILD.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"{dst_id}: {totals['moved']} of {totals['edges']} reaction-level edges moved "
          f"({totals['swaps']} swaps) in {len(per_pathway)} pathways")


if __name__ == "__main__":
    main(sys.argv)
