#!/usr/bin/env python3
"""specs/166: build the comparison arms from the FullTextMining composite-h12b
bundles.

- with:      the bundles as delivered (curated rows + literature rows/nodes).
- without:   the literature rows of logic_network.csv removed (matched on
             source_id, target_id, edge_type, pos_neg from literature_origin.csv)
             and the literature nodes (context_origin.csv) removed from
             nodes.csv and stid_to_uuid_mapping.csv: the generator's bundle as
             built.
- shuf<N>:   the literature edges rewired at random: within each class
             (edge_type, pos_neg, and_or) the TARGETS are permuted among the
             class's edges, so every node keeps its in- and out-degree per class
             and every edge keeps its sign and type; self-loops, duplicates and
             collisions with curated pairs are repaired by swapping targets between two
             edges of the class (degree-preserving); a build with any left fails.
             Curated rows untouched. Seed N.
- only-<slice>: `without` plus only the literature edges in one slice (layer,
             sign_source, group, prior_confidence bin); descriptive ablations.

Each arm is a catalog build directory with a BUILD.json whose build_id is the
directory name (the API and the benches check it).

Usage: build_arms_166.py SRC_DIR BUILDS_ROOT N_SHUFFLES
"""
import csv
import json
import random
import shutil
from collections import Counter
import sys
from pathlib import Path

PREFIX = "20261118-lit-h12b"
SKIP = {"cache"}


def read_rows(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return r.fieldnames, list(r)


def write_rows(path, fields, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def lit_keys(pdir):
    _, lo = read_rows(pdir / "literature_origin.csv")
    return {(r["source_id"], r["target_id"], r["edge_type"], r["pos_neg"]) for r in lo}


def copy_pathway(src, dst):
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name in SKIP:
            continue
        if f.is_file():
            shutil.copy2(f, dst / f.name)


def build_json(arm_dir, note, pathways):
    doc = {"build_id": arm_dir.name, "status": "complete",
           "generator_commit": "d6d026bb43f4c92d34bd51af501e2df937c58889",
           "reactome": {"release": 98},
           "variant": note, "pathways_built": sorted(pathways),
           "pathways_requested": sorted(pathways), "pathways_missing": []}
    (arm_dir / "BUILD.json").write_text(json.dumps(doc, indent=1))


def make_without(src, dst):
    copy_pathway(src, dst)
    keys = lit_keys(src)
    fields, rows = read_rows(src / "logic_network.csv")
    kept = [r for r in rows if (r["source_id"], r["target_id"], r["edge_type"], r["pos_neg"]) not in keys]
    write_rows(dst / "logic_network.csv", fields, kept)
    _, co = read_rows(src / "context_origin.csv")
    added = {r["uuid"] for r in co}
    for name, col in (("nodes.csv", "uuid"), ("stid_to_uuid_mapping.csv", "uuid")):
        f2, r2 = read_rows(src / name)
        write_rows(dst / name, f2, [r for r in r2 if r[col] not in added])
    for name in ("context_origin.csv", "literature_origin.csv", "literature_features.csv",
                 "literature_evidence.csv", "literature_support.csv", "literature_report.json"):
        (dst / name).unlink(missing_ok=True)
    return len(rows) - len(kept)


def pc_bin(pc):
    pc = float(pc)
    return "pclt04" if pc < 0.4 else ("pc04to06" if pc < 0.6 else "pcge06")


SLICES = {
    "verified": lambda o, f: o["layer"] == "verified",
    "extracted": lambda o, f: o["layer"] == "extracted",
    "stated": lambda o, f: o["sign_source"] == "stated",
    "inferred": lambda o, f: o["sign_source"] == "inferred",
    "signna": lambda o, f: o["sign_source"] == "na",
    "binding": lambda o, f: f["group"] == "binding/other",
    "regulation": lambda o, f: f["group"] == "regulation",
    "expression": lambda o, f: f["group"] == "expression",
    "pclt04": lambda o, f: pc_bin(o["prior_confidence"]) == "pclt04",
    "pc04to06": lambda o, f: pc_bin(o["prior_confidence"]) == "pc04to06",
    "pcge06": lambda o, f: pc_bin(o["prior_confidence"]) == "pcge06",
}


def make_slice(src, dst, name):
    """`without` plus the literature edges of one slice. Literature nodes stay
    (an untouched one is inert); dropped rows are those outside the slice."""
    copy_pathway(src, dst)
    keep = SLICES[name]
    _, feats = read_rows(src / "literature_features.csv")
    fe = {r["edge_key"]: r for r in feats}
    _, lo = read_rows(src / "literature_origin.csv")
    drop = {(r["source_id"], r["target_id"], r["edge_type"], r["pos_neg"])
            for r in lo if not keep(r, fe[r["edge_key"]])}
    fields, rows = read_rows(src / "logic_network.csv")
    kept = [r for r in rows if (r["source_id"], r["target_id"], r["edge_type"], r["pos_neg"]) not in drop]
    write_rows(dst / "logic_network.csv", fields, kept)
    return len(lo) - len(drop)


def src_dir_name(d):
    return f"{d.parent.name}/{d.name}"


def make_shuffle(src, dst, seed):
    copy_pathway(src, dst)
    keys = lit_keys(src)
    fields, rows = read_rows(src / "logic_network.csv")
    lit_idx = [i for i, r in enumerate(rows)
               if (r["source_id"], r["target_id"], r["edge_type"], r["pos_neg"]) in keys]
    rng = random.Random(seed)
    by_class = {}
    for i in lit_idx:
        r = rows[i]
        by_class.setdefault((r["edge_type"], r["pos_neg"], r["and_or"]), []).append(i)
    lit_set = set(lit_idx)
    curated = {(r["source_id"], r["target_id"]) for i, r in enumerate(rows) if i not in lit_set}
    unresolved = 0
    for cls, idx in by_class.items():
        srcs = [rows[i]["source_id"] for i in idx]
        tgt = [rows[i]["target_id"] for i in idx]
        rng.shuffle(tgt)
        used = Counter(zip(srcs, tgt))

        def ok(k, t):
            return srcs[k] != t and (srcs[k], t) not in curated and used[(srcs[k], t)] == 0

        def bad(k):
            return (srcs[k] == tgt[k] or (srcs[k], tgt[k]) in curated
                    or used[(srcs[k], tgt[k])] > 1)

        # Swap repair: a target swap between two edges of one class keeps every
        # node's in- and out-degree in the class.
        for _ in range(20):
            todo = [k for k in range(len(idx)) if bad(k)]
            if not todo:
                break
            for k in todo:
                if not bad(k):
                    continue
                for _ in range(200):
                    j = rng.randrange(len(idx))
                    if j == k:
                        continue
                    used[(srcs[k], tgt[k])] -= 1
                    used[(srcs[j], tgt[j])] -= 1
                    if ok(k, tgt[j]) and ok(j, tgt[k]) and not (srcs[k] == srcs[j]):
                        tgt[k], tgt[j] = tgt[j], tgt[k]
                        used[(srcs[k], tgt[k])] += 1
                        used[(srcs[j], tgt[j])] += 1
                        break
                    used[(srcs[k], tgt[k])] += 1
                    used[(srcs[j], tgt[j])] += 1
        unresolved += sum(bad(k) for k in range(len(idx)))
        for i, t in zip(idx, tgt):
            rows[i] = dict(rows[i], target_id=t)
    write_rows(dst / "logic_network.csv", fields, rows)
    if unresolved:
        raise SystemExit(f"{src_dir_name(dst)}: {unresolved} shuffled edges could not be repaired")
    return len(lit_idx), unresolved


def main():
    src_root, builds, n_shuf = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    pathways = sorted(p.name for p in src_root.iterdir() if (p / "logic_network.csv").exists())
    arms = {"with": None, "without": None, **{f"shuf{s:02d}": s for s in range(1, n_shuf + 1)},
            **{f"only-{k}": k for k in SLICES}}
    for arm, seed in arms.items():
        arm_dir = builds / f"{PREFIX}-{arm}"
        if arm_dir.exists():
            shutil.rmtree(arm_dir)
        arm_dir.mkdir(parents=True)
        info = []
        for p in pathways:
            if arm == "with":
                copy_pathway(src_root / p, arm_dir / p)
            elif arm == "without":
                info.append(make_without(src_root / p, arm_dir / p))
            elif arm.startswith("only-"):
                info.append(make_slice(src_root / p, arm_dir / p, seed))
            else:
                info.append(make_shuffle(src_root / p, arm_dir / p, seed)[0])
        build_json(arm_dir, f"specs/166 arm {arm}", pathways)
        print(f"{arm_dir.name}: {len(pathways)} pathways"
              + (f"; literature rows handled {sum(info)}" if info else ""))


if __name__ == "__main__":
    main()
