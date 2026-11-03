#!/usr/bin/env python3
"""specs/121: the CPPA train / test split, fixed before any CPPA value is read.

Reads only the CPPA delta-p treatment summaries (drug, cell line) and the drug names of the datasets already used.
Single-compound samples only. Seen drugs (in korkut2015, decryptM, HPN-DREAM or LINCS P100) and stimuli (ligands,
serum, vehicle) are training only. Of the rest, 40% of drugs are test drugs; 25% of cell lines are test cell lines;
seed 121. Writes the split (kept outside the repository, under the CPPA terms) and prints its SHA-256.

Usage: split_121.py CPPA_UNPACKED KORKUT_DRUGS_TSV DECRYPTM_CASES_TSV P100_INST_GZ OUT_TSV
"""
import csv
import glob
import gzip
import hashlib
import random
import re
import sys

SEED, TEST_DRUGS, TEST_CELLS = 121, 0.40, 0.25
# HPN-DREAM's four inhibitor conditions (Hill et al. 2016)
HPN = ("GSK690693", "GSK1120212", "PD173074", "TEMSIROLIMUS")
STIMULI = {"PBS", "SERUM", "SERUMGLUCOSESTARVATION", "EGF", "FGF1", "FGF2", "HGF", "IGF1", "INSULIN", "NRG1"}
# one name per compound, so a synonym cannot make a seen drug look unseen
SYNONYMS = {"GLEEVEC": "IMATINIB", "GSK2118436": "DABRAFENIB", "GSK1120212": "TRAMETINIB", "PD325901": "PD0325901",
            "BAFILOMYCIN": "BAFILOMYCINA1"}


def norm(name):
    n = re.sub(r"[^A-Z0-9]", "", name.upper())
    return SYNONYMS.get(n, n)


def cppa_singles(unpacked):
    """[(drug, cell line)] for samples with exactly one compound."""
    out = []
    for f in sorted(glob.glob(f"{unpacked}/delta-p_TreatmentSummary/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t"):
            names = [(r.get(f"compound_name_{i}") or "").strip() for i in range(1, 5)]
            names = [x for x in names if x and x.upper() != "NA"]
            if len(names) == 1:
                out.append((norm(names[0]), norm(r["cell_line_name"])))
    return out


def seen_drugs(korkut_drugs, decryptm_cases, p100_inst):
    seen = {norm(h) for h in HPN}
    seen |= {norm(r["drug"]) for r in csv.DictReader((l for l in open(korkut_drugs) if not l.startswith("#")),
                                                       delimiter="\t")}
    seen |= {norm(r["drug"]) for r in csv.DictReader(open(decryptm_cases), delimiter="\t")}
    with gzip.open(p100_inst, "rt") as f:
        seen |= {norm(r["pert_iname"]) for r in csv.DictReader(f, delimiter="\t") if r["pert_type"] == "trt_cp"}
    return seen


def split(pairs, seen, seed=SEED):
    drugs = sorted({d for d, _ in pairs})
    cells = sorted({c for _, c in pairs})
    rng = random.Random(seed)
    free = [d for d in drugs if d not in seen and d not in STIMULI]
    test_drugs = set(rng.sample(free, round(TEST_DRUGS * len(free))))
    test_cells = set(rng.sample(cells, round(TEST_CELLS * len(cells))))
    rows = []
    for d in drugs:
        kind = "stimulus" if d in STIMULI else "seen" if d in seen else "test" if d in test_drugs else "train"
        rows.append(("drug", d, kind))
    rows += [("cell_line", c, "test" if c in test_cells else "train") for c in cells]
    return rows


def main():
    unpacked, korkut, decryptm, p100, out = sys.argv[1:6]
    rows = split(cppa_singles(unpacked), seen_drugs(korkut, decryptm, p100))
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["unit", "name", "role"])
        w.writerows(rows)
    counts = {}
    for unit, _, role in rows:
        counts[f"{unit}:{role}"] = counts.get(f"{unit}:{role}", 0) + 1
    print(counts)
    print("sha256", hashlib.sha256(open(out, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
