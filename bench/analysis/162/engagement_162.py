#!/usr/bin/env python3
"""specs/162: per-(drug, cell line) engagement on CPPA, from the drug's OWN
specs/139 markers only.

e = min(1, 2^m), with m the mean marker delta-p (log2) over every marker
antibody present in the pair's single-compound samples (all times and doses).
Non-stimulus drugs with at least one marker only. Reads CPPA through
cppa_data_121.load for every portion (all CPPA is FIT after specs/121's single
scoring, so final=True only lifts the development guard). It uses only marker
antibodies; no other value is read here. Output stays OUTSIDE the repository
(CPPA terms).

Usage: engagement_162.py CPPA_UNPACKED SPLIT_TSV MARKERS_TSV DRUG_TARGETS_TSV OUT_TSV
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import cppa_data_121 as cd  # noqa: E402
import dataset_roles  # noqa: E402


def markers(markers_tsv, drug_targets_tsv):
    stim = {r["drug"] for r in csv.DictReader(open(drug_targets_tsv), delimiter="\t")
            if r["role"] == "stimulus"}
    out = {}
    for r in csv.DictReader(open(markers_tsv), delimiter="\t"):
        ms = [m for m in r["markers"].split(",") if m.strip()]
        if ms and r["drug"] not in stim:
            out[r["drug"]] = set(ms)
    return out


def engagement(unpacked, split_tsv, mk):
    acc = defaultdict(list)
    for portion in cd.PORTIONS:
        for s in cd.load(unpacked, split_tsv, portion, final=True):
            ms = mk.get(s["drug"])
            if not ms:
                continue
            acc[(s["drug"], s["cell_line"])].extend(v for a, v in s["values"].items() if a in ms)
    return {k: min(1.0, 2.0 ** (sum(v) / len(v))) for k, v in acc.items() if v}


def main():
    unpacked, split_tsv, markers_tsv, drug_targets_tsv, out_tsv = sys.argv[1:6]
    for ds in ("cppa_train", "cppa_test"):       # specs/051: engagement is fitted on every portion
        dataset_roles.require_role(ds, "fit")
    mk = markers(markers_tsv, drug_targets_tsv)
    eng = engagement(unpacked, split_tsv, mk)
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["drug", "cell_line", "engagement"])
        for (d, c), e in sorted(eng.items()):
            w.writerow([d, c, f"{e:.4f}"])
    print(f"{len(eng)} (drug, cell line) pairs over {len({d for d, _ in eng})} drugs -> {out_tsv}")


if __name__ == "__main__":
    main()
