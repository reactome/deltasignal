#!/usr/bin/env python3
"""specs/169: measured transcript changes for split-A perturbed genes only (split B is never read).

Per perturbed gene G: among its non-control rows with knockdown fold_expr <= 0.5 and num_cells_filtered >= 25, the
row with the most cells (none: G is not engaged and is dropped). Control level of each measured gene = the mean of
X over the core_control rows. Measured fold of T under G = X[row, T] / control[T] (control 0: dropped).
Labels: `measured_50` at korkut's convention (|log2 fold| >= 0.585, i.e. 1.5x) and `measured_15` at the model's
cutoffs (< 0.85 down, >= 1.15 up).

Usage: truth_169.py H5AD_RAW PREDICTIONS_TSV OUT_TSV     (only (G, T) pairs present in the predictions are written)
"""
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent))
from common_169 import split_of, gene_of, h5_col, obs_index  # noqa: E402

KD_MAX, CELLS_MIN = 0.5, 25


def label(fold, lo, hi):
    return "down" if fold < lo else ("up" if fold >= hi else "not")


def choose_rows(genes, kd, cells, ids, ctrl):
    """{gene: row index} for engaged genes; the most cells wins."""
    best = {}
    for i, (oid, c) in enumerate(zip(ids, ctrl)):
        g = gene_of(oid)
        if c or g not in genes or not (kd[i] <= KD_MAX and cells[i] >= CELLS_MIN):
            continue
        if g not in best or cells[i] > cells[best[g]]:
            best[g] = i
    return best


def main():
    import h5py
    import numpy as np
    import dataset_roles
    dataset_roles.require_role("replogle2022_a", "score")
    h5, preds, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    pairs = {(r["gene"], r["target"]) for r in csv.DictReader(open(preds), delimiter="\t")}
    genes = {g for g, _ in pairs}
    if any(split_of(g) != "A" for g in genes):
        raise SystemExit("a split-B gene reached the truth reader: refusing (split B is sealed)")
    with h5py.File(h5, "r") as f:
        ids, ctrl = obs_index(f), h5_col(f, "obs", "core_control").astype(bool)
        kd, cells = h5_col(f, "obs", "fold_expr"), h5_col(f, "obs", "num_cells_filtered")
        col = {g: j for j, g in enumerate(h5_col(f, "var", "gene_name"))}
        rows = choose_rows(genes, kd, cells, ids, ctrl)
        X = f["X"]
        control = np.asarray(X[np.where(ctrl)[0], :]).mean(axis=0)
        out_rows = []
        for g in sorted(rows):
            x = np.asarray(X[rows[g], :])
            for t in sorted(t for gg, t in pairs if gg == g):
                j = col.get(t)
                if j is None or control[j] <= 0:
                    continue
                fold = float(x[j] / control[j])
                out_rows.append({"gene": g, "target": t, "measured_fold": f"{fold:.6g}",
                                 "measured_50": label(fold, 2 ** -0.585, 2 ** 0.585),
                                 "measured_15": label(fold, 0.85, 1.15),
                                 "knockdown": f"{kd[rows[g]]:.4g}", "cells": int(cells[rows[g]])})
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]) if out_rows else ["gene"], delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"{len(rows)} engaged genes of {len(genes)}; {len(out_rows)} (gene, target) truths -> {out}")


if __name__ == "__main__":
    main()
