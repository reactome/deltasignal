#!/usr/bin/env python3
"""Where DeltaSignal loses to MP-BioPath on the experimental axis, case by case.

MP-BioPath's per-case calls are published in
mp-biopath-pathways/analysis_results/MP_BioPathReactomePathwayAccuracy-Tests-updated.tsv
(`mp_biopath_state` against `experimental_value`, 849 cases, 643 correct).
This joins them with an arm's experimental_cases.tsv and splits the cases into
both right / only MP-BioPath / only DeltaSignal / neither, by pathway,
perturbation, our failure category and (truth, our call).

A case that does not join, or whose experimental value disagrees with our
expected label, is an error, not a skip: a silent mis-join would make the gap
look like whatever the join happened to keep.

Usage: experimental_gap.py ARM_RESULTS_DIR [--mpb TSV] [--out cases.tsv]
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

DEFAULT_MPB = (Path.home() / "gitroot" / "mp-biopath-pathways" / "analysis_results"
               / "MP_BioPathReactomePathwayAccuracy-Tests-updated.tsv")
LABEL = {"0": "DOWN", "1": "NORMAL", "2": "UP"}


def load_mpb(path: Path) -> dict[tuple, tuple[str, str]]:
    """(pathway, gene, direction, key_output) -> (experimental, mp-biopath state)."""
    out: dict[tuple, tuple[str, str]] = {}
    with path.open(newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["experimental_value"] in ("NA", ""):
                continue
            gene, direction = r["scenario"].rsplit("_", 1)
            key = (r["pathway_name"], gene, direction, r["keyoutput_id"])
            if key in out:
                raise SystemExit(f"duplicate MP-BioPath row {key}")
            out[key] = (r["experimental_value"], r["mp_biopath_state"])
    return out


def split(mpb: dict, ours_rows: list[dict]) -> list[dict]:
    rows = []
    for o in ours_rows:
        key = (o["pathway"], o["gene"], o["direction"], o["key_output"])
        if key not in mpb:
            raise SystemExit(f"case not in MP-BioPath's table: {key}")
        truth, them = mpb[key]
        if truth != o["expected"]:
            raise SystemExit(f"experimental value {truth} != our expected {o['expected']} for {key}")
        us = o["predicted"] if o["valid"] == "1" else "invalid"
        rows.append({"pathway": key[0], "gene": key[1], "direction": key[2], "key_output": key[3],
                     "truth": truth, "mpbiopath": them, "deltasignal": us,
                     "quadrant": ("MPB" if them == truth else "-") + "/" + ("DS" if us == truth else "-"),
                     "category": o.get("category", "")})
    missing = set(mpb) - {(r["pathway"], r["gene"], r["direction"], r["key_output"]) for r in rows}
    if missing:
        raise SystemExit(f"{len(missing)} MP-BioPath cases have no row in the arm, e.g. {sorted(missing)[0]}")
    return rows


def report(rows: list[dict]) -> str:
    out = []
    q = collections.Counter(r["quadrant"] for r in rows)
    mpb_right = q["MPB/DS"] + q["MPB/-"]; ds_right = q["MPB/DS"] + q["-/DS"]
    out.append(f"cases {len(rows)}: MP-BioPath {mpb_right}, DeltaSignal {ds_right}, gap {mpb_right - ds_right:+d}")
    out.append(f"  both right {q['MPB/DS']}, only MP-BioPath {q['MPB/-']}, "
               f"only DeltaSignal {q['-/DS']}, neither {q['-/-']}; either {len(rows) - q['-/-']}")
    pw = collections.defaultdict(collections.Counter)
    for r in rows:
        pw[r["pathway"]][r["quadrant"]] += 1
    out.append("per pathway: only-MPB | only-DS | net")
    for p, c in sorted(pw.items(), key=lambda x: -(x[1]["MPB/-"] - x[1]["-/DS"])):
        out.append(f"  {p[:40]:40s} {c['MPB/-']:4d} | {c['-/DS']:4d} | {c['MPB/-'] - c['-/DS']:+4d}  (n={sum(c.values())})")
    lost = [r for r in rows if r["quadrant"] == "MPB/-"]
    out.append("only-MPB by our category: " + str(collections.Counter(r["category"] for r in lost).most_common()))
    out.append("only-MPB by (truth, our call): " + str(collections.Counter(
        (LABEL.get(r["truth"]), LABEL.get(r["deltasignal"], r["deltasignal"])) for r in lost).most_common()))
    pert = collections.Counter((r["pathway"][:24], r["gene"], "KO" if r["direction"] == "0" else "UP") for r in lost)
    out.append("only-MPB by perturbation: " + str(pert.most_common(15)))
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arm", type=Path)
    ap.add_argument("--mpb", type=Path, default=DEFAULT_MPB)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    with (a.arm / "experimental_cases.tsv").open(newline="") as handle:
        ours = list(csv.DictReader(handle, delimiter="\t"))
    rows = split(load_mpb(a.mpb), ours)
    print(report(rows))
    if a.out:
        with a.out.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
