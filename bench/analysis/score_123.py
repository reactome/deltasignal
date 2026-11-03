#!/usr/bin/env python3
"""specs/123 M1: CPPA training per-drug correlation (fold arm, specs/121's measure) for several drug-action arms, and
each arm minus the control with a 95% paired bootstrap over drugs (2,000 resamples, seed 123).

Each arm is NAME=CASES_GLOB (cppa_bench case logs). The first arm is the control. With --labels TABLE, the same is
reported for the drugs whose CPPA label is known (not `unknown`). Training portion only (cppa_data_121.load).

Usage: score_123.py CPPA_UNPACKED SPLIT_TSV READOUTS_TSV [--labels TABLE] NAME=GLOB NAME=GLOB [...]
"""
import csv
import glob
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_data_121  # noqa: E402
import score_121  # noqa: E402


def paired(drug_corr, control, arm, drugs, seed=123, n_boot=2000):
    """(mean arm - control, 95% interval) over `drugs`, resampling drugs."""
    diff = {d: drug_corr[arm][d] - drug_corr[control][d] for d in drugs}
    rng = random.Random(seed)
    boots = sorted(sum(diff[rng.choice(drugs)] for _ in drugs) / len(drugs) for _ in range(n_boot))
    return sum(diff.values()) / len(drugs), [boots[int(0.025 * n_boot)], boots[int(0.975 * n_boot) - 1]]


def report(drug_corr, names, drugs):
    out = {"drugs": len(drugs), "mean": {a: sum(drug_corr[a][d] for d in drugs) / len(drugs) for a in names}}
    for a in names[1:]:
        m, ci = paired(drug_corr, names[0], a, drugs)
        out[f"{a}-{names[0]}"] = {"mean": m, "ci95": ci, "better": sum(drug_corr[a][d] > drug_corr[names[0]][d]
                                                                          for d in drugs),
                                  "worse": sum(drug_corr[a][d] < drug_corr[names[0]][d] for d in drugs)}
    return out


def main():
    args = sys.argv[1:]
    if len(args) < 5:
        raise SystemExit(__doc__)
    unpacked, split_tsv, readouts_tsv, *rest = args
    labels = None
    if rest and rest[0] == "--labels":
        with open(rest[1]) as f:
            labels = {r["drug"]: r["label"] for r in csv.DictReader((x for x in f if not x.startswith("#")),
                                                                    delimiter="\t") if r["axis"] == "cppa"}
        rest = rest[2:]
    arms = [a.split("=", 1) for a in rest]
    readouts = {r["readout"] for r in csv.DictReader((x for x in open(readouts_tsv) if not x.startswith("#")),
                                                     delimiter="\t") if r["kind"] != "exclude"}
    samples = cppa_data_121.load(unpacked, split_tsv, "train")
    drug_corr = {}
    for name, pattern in arms:
        files = sorted(glob.glob(pattern))
        if not files:
            raise SystemExit(f"{name}: no case log matches {pattern}")
        preds = score_121.predictions(files)["fold"]
        drug_corr[name] = score_121.per_drug(samples, {"fold": preds}, readouts)["fold"]
    names = [n for n, _ in arms]
    drugs = sorted(set.intersection(*(set(drug_corr[a]) for a in names)))
    out = {"portion": "train", "all": report(drug_corr, names, drugs)}
    if labels is not None:
        known = [d for d in drugs if labels.get(d, "unknown") != "unknown"]
        out["known_label"] = report(drug_corr, names, known)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
