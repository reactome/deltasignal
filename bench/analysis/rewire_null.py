#!/usr/bin/env python3
"""specs/051: score the rewired-network null against its matched control.

    rewire_null.py WIRE0_DIR SEED_DIR [SEED_DIR ...] [--canonical DIR] [--json OUT]

Each DIR is a run_arm.sh output holding curator_cases.tsv and
experimental_cases.tsv. Cases are paired across every arm (the same pathway,
gene, direction and readout), so an arm that lost a case to a solve failure
cannot shift the comparison. The decision rule was registered in
specs/051-dataset-roles/spec.md before any arm ran:

  the wiring carries signal on an axis if wire0's macro-F1 exceeds every
  seed's AND the 95% interval of (wire0 - seed mean) correctness, bootstrapped
  over genes, excludes 0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from holdout_report import TUNING_PATHWAYS, load, macro_f1  # noqa: E402
from protocol_stats import cluster_bootstrap_ci  # noqa: E402


def paired(arms: list[dict]) -> list:
    """Keys present in every arm, in a stable order."""
    keys = set(arms[0])
    for a in arms[1:]:
        keys &= set(a)
    return sorted(keys)


def slices(keys, axis):
    if axis == "curator":
        return {"tuning": [k for k in keys if k[0] in TUNING_PATHWAYS],
                "held-out": [k for k in keys if k[0] not in TUNING_PATHWAYS]}
    return {"all": list(keys)}


def compare(wire0: dict, seeds: list[dict], keys: list, canonical: dict | None = None) -> dict:
    rows = []
    for k in keys:
        exp = wire0[k]["expected"]
        rows.append({"gene": k[1], "expected": exp,
                     "real": wire0[k]["predicted"],
                     "seeds": [s[k]["predicted"] for s in seeds]})
    if not rows:
        return {"n": 0}

    def acc(rs, pick):
        return sum(pick(r) == r["expected"] for r in rs) / len(rs)

    def mf1(rs, pick):
        return macro_f1([{"predicted": pick(r), "expected": r["expected"]} for r in rs])

    def diff(rs):
        return sum((r["real"] == r["expected"]) -
                   sum(p == r["expected"] for p in r["seeds"]) / len(r["seeds"]) for r in rs) / len(rs)

    def mf1_diff(rs):
        seed_m = sum(mf1(rs, lambda r, i=i: r["seeds"][i]) for i in range(len(seeds))) / len(seeds)
        return mf1(rs, lambda r: r["real"]) - seed_m

    seed_mf1 = [mf1(rows, lambda r, i=i: r["seeds"][i]) for i in range(len(seeds))]
    seed_acc = [acc(rows, lambda r, i=i: r["seeds"][i]) for i in range(len(seeds))]
    out = {
        "n": len(rows),
        "genes": len({r["gene"] for r in rows}),
        "wire0": {"accuracy": acc(rows, lambda r: r["real"]), "macro_f1": mf1(rows, lambda r: r["real"])},
        "seeds": {"accuracy": seed_acc, "macro_f1": seed_mf1},
        "seed_mean": {"accuracy": sum(seed_acc) / len(seed_acc), "macro_f1": sum(seed_mf1) / len(seed_mf1)},
        "correct_diff_ci95_by_gene": list(cluster_bootstrap_ci(rows, diff, lambda r: r["gene"])),
        "macro_f1_diff_ci95_by_gene": list(cluster_bootstrap_ci(rows, mf1_diff, lambda r: r["gene"], n_boot=500)),
    }
    if canonical is not None:
        cr = [canonical[k] for k in keys if k in canonical]
        out["canonical"] = {"n": len(cr), "accuracy": sum(r["predicted"] == r["expected"] for r in cr) / len(cr),
                            "macro_f1": macro_f1(cr)}
    lo = out["correct_diff_ci95_by_gene"][1]
    out["carries_signal"] = out["wire0"]["macro_f1"] > max(seed_mf1) and lo > 0
    return out


def run(wire0_dir: Path, seed_dirs: list[Path], canonical_dir: Path | None) -> dict:
    result = {"wire0": str(wire0_dir), "seeds": [str(d) for d in seed_dirs],
              "canonical": str(canonical_dir) if canonical_dir else None}
    for axis in ("curator", "experimental"):
        name = f"{axis}_cases.tsv"
        w0 = load(wire0_dir / name)
        ss = [load(d / name) for d in seed_dirs]
        can = load(canonical_dir / name) if canonical_dir else None
        keys = paired([w0, *ss])
        result[axis] = {"unpaired": {"wire0": len(w0), "seeds": [len(s) for s in ss], "paired": len(keys)}}
        for sl, ks in slices(keys, axis).items():
            result[axis][sl] = compare(w0, ss, ks, can)
    return result


def show(res: dict) -> None:
    for axis in ("curator", "experimental"):
        print(f"\n{axis}: cases {res[axis]['unpaired']}")
        for sl, m in res[axis].items():
            if sl == "unpaired" or not m.get("n"):
                continue
            lo, hi = m["correct_diff_ci95_by_gene"][1:]
            flo, fhi = m["macro_f1_diff_ci95_by_gene"][1:]
            seeds = " ".join(f"{x:.3f}" for x in m["seeds"]["macro_f1"])
            can = f"  canonical mF1 {m['canonical']['macro_f1']:.4f}" if "canonical" in m else ""
            print(f"  {sl:9s} n={m['n']:6d} genes={m['genes']:4d}  wire0 acc {m['wire0']['accuracy']:.4f} "
                  f"mF1 {m['wire0']['macro_f1']:.4f} | seeds acc {m['seed_mean']['accuracy']:.4f} "
                  f"mF1 [{seeds}]{can}")
            print(f"            wire0 - seeds: correct [{lo:+.4f}, {hi:+.4f}], mF1 [{flo:+.4f}, {fhi:+.4f}]"
                  f"  -> carries signal: {m['carries_signal']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wire0", type=Path)
    ap.add_argument("seeds", type=Path, nargs="+")
    ap.add_argument("--canonical", type=Path)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args(argv)
    res = run(a.wire0, a.seeds, a.canonical)
    show(res)
    if a.json:
        a.json.write_text(json.dumps(res, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
