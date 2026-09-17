#!/usr/bin/env python3
"""Report accuracy split into the tuning pathways and the held-out rest.

The MP-BioPath paper tuned on ten pathways and reported on the others. This
project drifted off that protocol: the ten-pathway set kept reversing decisions
that held at scale, so decisions moved to the wide curator set -- which is the
set we then report accuracy on. That fixed the misleading and silently turned
the test set into the training set.

This tool restores the split at REPORTING time, over an existing case dump, so
no benchmark has to be re-run. Use it for any number that leaves the building.

Two things it will show you:

  * the headline all-pathways figure is dragged down by the tuning ten, which
    are hard for everyone (TP53, WNT, PIP3, cell cycle);
  * whether a config decision made on the wide set also wins on pathways it
    never saw -- pass --compare with a second dump to check that directly.
    A decision that only wins on the tuning half is overfitting.
"""
from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

# The paper's tuning set. Names as they appear in the curator files.
TUNING_PATHWAYS = {
    "Cell_Cycle_Checkpoints",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA_",
    "Mitotic_G2-G2_M_phases",
    "Mitotic_Prophase",
    "Mitotic_G1-G1_S_phases",
    "PIP3_activates_AKT_signaling",
    "RAF_MAP_kinase_cascade",
    "Signaling_by_ERBB2",
    "Signaling_by_WNT",
    "S_Phase",
    "Transcriptional_Regulation_by_TP53",
}

LABELS = ("0", "1", "2")


def load(path: Path) -> dict:
    with path.open(newline="") as fh:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(fh, delimiter="\t")
                if r["expected"] in LABELS}


def macro_f1(rows) -> float:
    f1s = []
    for lab in LABELS:
        tp = sum(1 for r in rows if r["predicted"] == lab and r["expected"] == lab)
        fp = sum(1 for r in rows if r["predicted"] == lab and r["expected"] != lab)
        fn = sum(1 for r in rows if r["predicted"] != lab and r["expected"] == lab)
        p = tp / (tp + fp) if tp + fp else 0.0
        rc = tp / (tp + fn) if tp + fn else 0.0
        f1s.append(2 * p * rc / (p + rc) if p + rc else 0.0)
    return sum(f1s) / len(f1s)


def report(name: str, cases: dict) -> None:
    groups = {"TUNING (the paper's ten)": [], "HELD-OUT (report this)": []}
    for k, r in cases.items():
        key = "TUNING (the paper's ten)" if k[0] in TUNING_PATHWAYS else "HELD-OUT (report this)"
        groups[key].append(r)
    print(f"\n=== {name} ===")
    print(f"{'split':<26}{'pathways':>9}{'cases':>8}{'accuracy':>10}{'macro-F1':>10}")
    for g, rows in groups.items():
        if not rows:
            continue
        pw = len({r["pathway"] for r in rows})
        ok = sum(1 for r in rows if r["predicted"] == r["expected"])
        print(f"{g:<26}{pw:>9}{len(rows):>8}{ok/len(rows):>10.4f}{macro_f1(rows):>10.4f}")
    allrows = [r for rs in groups.values() for r in rs]
    ok = sum(1 for r in allrows if r["predicted"] == r["expected"])
    print(f"{'(all pathways)':<26}{len({r['pathway'] for r in allrows}):>9}"
          f"{len(allrows):>8}{ok/len(allrows):>10.4f}{macro_f1(allrows):>10.4f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--compare", type=Path,
                    help="A second dump. Reports whether the change wins on the "
                         "HELD-OUT half as well as the tuning half -- a change "
                         "that only wins on tuning is overfitting.")
    a = ap.parse_args()

    base = load(a.cases)
    report(a.cases.name, base)

    if a.compare:
        arm = load(a.compare)
        report(a.compare.name, arm)
        shared = sorted(set(base) & set(arm))
        print(f"\n=== does it generalise?  {a.compare.name} vs {a.cases.name} ===")
        print(f"{'split':<26}{'cases':>8}{'baseline':>10}{'arm':>10}{'net':>8}")
        for label, keys in (
            ("TUNING (the paper's ten)", [k for k in shared if k[0] in TUNING_PATHWAYS]),
            ("HELD-OUT (report this)", [k for k in shared if k[0] not in TUNING_PATHWAYS]),
        ):
            if not keys:
                continue
            b = sum(1 for k in keys if base[k]["predicted"] == base[k]["expected"])
            m = sum(1 for k in keys if arm[k]["predicted"] == arm[k]["expected"])
            print(f"{label:<26}{len(keys):>8}{b/len(keys):>10.4f}{m/len(keys):>10.4f}{m-b:>+8d}")
        print("\nA decision that wins on TUNING but not HELD-OUT is overfitting.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
