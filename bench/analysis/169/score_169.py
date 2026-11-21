#!/usr/bin/env python3
"""specs/169: the specs/167 numbers on the genetic panel, per truth threshold, bootstrap over perturbed genes.

Usage: score_169.py PREDICTIONS_TSV TRUTH_TSV OUT_JSON
"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "167"))
import base_rate as br  # noqa: E402
from fold_rules import record_fold  # noqa: E402


def cases(preds, truth, col, matched=None):
    """{(pathway, gene, target): (call, truth, fold)}; fold 1 at record precision = unreached.
    matched: re-call the predicted fold at the truth's own threshold (descriptive)."""
    out = {}
    for r in preds:
        t = truth.get((r["gene"], r["target"]))
        if t is None:
            continue
        fold = record_fold(float(r["fold"]))
        call = r["call"]
        if matched is not None:
            lo, hi = matched
            call = "down" if fold < lo else ("up" if fold >= hi else "not")
        out[(r["pathway"], r["gene"], r["target"])] = (call, t[col], fold)
    return out


def main():
    preds = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    truth = {(r["gene"], r["target"]): r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
    unit = lambda k: k[1]  # noqa: E731
    res = {
        "truth_1.5x": br.stats(cases(preds, truth, "measured_50"), unit),
        "truth_15pct": br.stats(cases(preds, truth, "measured_15"), unit),
        "truth_1.5x_model_called_at_1.5x": br.stats(cases(preds, truth, "measured_50", (2 ** -0.585, 2 ** 0.585)), unit),
    }
    b = res["truth_1.5x"]["base_rate_reached"]
    res["reading"] = ("none: nothing reached" if b is None else
                      "genetic perturbations carry through curated routes (gene axes' precision is not only selection)"
                      if b >= 0.5 else
                      "comparable to the drug panels: the gene axes' precision is readout selection"
                      if b <= 0.25 else "intermediate: neither reading holds")
    Path(sys.argv[3]).write_text(json.dumps(res, indent=1))
    for k, r in res.items():
        print(k, r if isinstance(r, str) else {x: r[x] for x in ("cases", "units", "reached", "base_rate_reached",
              "base_rate_unreached", "call_rate", "precision_change", "lift", "lift_ci", "declined_change_rate")})


if __name__ == "__main__":
    main()
