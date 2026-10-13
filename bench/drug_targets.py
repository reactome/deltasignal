"""specs/051 R4: drug -> target genes from measured affinity.

A protein is a target of a drug if the Kuster kinobeads deposit in ChEMBL 37
(Klaeger 2017) has a measured `Kd apparent` (relation '=') at or below the
cutoff. Only single-protein targets count. A drug the deposit does not cover,
or that has no target at the cutoff, keeps its hand-written targets and is
marked `hand`. The rule was registered in specs/051 research R4 before any
value was read. Amendment A1 (after the map was seen, before any score): a
hand target whose gene the deposit never measured is kept beside the
measured targets, marked `hand_unmeasured`.

Usage: drug_targets.py DATA_DIR [CUTOFF_NM] > drug_targets.tsv
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIMARY_CUTOFF_NM = 1000.0
SENSITIVITY_NM = (100.0, 10000.0)


def load_molecules(path=HERE / "drug_molecules.tsv"):
    with open(path, newline="") as f:
        return {r["drug"]: r["chembl_molecule"] for r in csv.DictReader(f, delimiter="\t")}


def load_targets(path, acc_gene_path):
    """target_chembl_id -> (accession, gene) for SINGLE PROTEIN human targets; others -> None.

    The gene is Reactome's own name for the accession (`reactome97_acc_gene.tsv`,
    isoforms excluded), so it matches what the bench resolves; ChEMBL's synonym
    lists carry old aliases. An accession Reactome lacks maps to ("", None)."""
    with open(acc_gene_path, newline="") as f:
        acc_gene = {r["acc"]: r["gene"] for r in csv.DictReader(f, delimiter="\t")}
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ok = r["target_type"] == "SINGLE PROTEIN" and r["organism"] == "Homo sapiens"
            out[r["target_chembl_id"]] = (r["accession"], acc_gene.get(r["accession"])) if ok else None
    return out


def build(activities, targets, molecules, hand, cutoff_nm=PRIMARY_CUTOFF_NM):
    """Rows (drug, gene, uniprot, kd_nm, source) and a tally. `activities` are
    the deposit's rows as dicts; `hand` is the old {drug: [genes]} table."""
    by_mol = defaultdict(dict)
    measured = set()            # genes with ANY record, for any drug, either relation
    tally = defaultdict(int)
    for a in activities:
        if a["standard_type"] != "Kd" or a["standard_units"] != "nM":
            tally["other_type"] += 1
            continue
        t = targets.get(a["target_chembl_id"])
        if t is None:
            tally["not_single_protein"] += 1
            continue
        acc, gene = t
        if gene is None:
            tally["not_in_reactome"] += 1
            continue
        measured.add(gene)
        if a["standard_relation"] != "=":
            continue
        kd = float(a["standard_value"])
        prev = by_mol[a["molecule_chembl_id"]].get(gene)
        if prev is None or kd < prev[1]:
            by_mol[a["molecule_chembl_id"]][gene] = (acc, kd)
    rows = []
    for drug in sorted(molecules):
        mol = molecules[drug]
        hits = sorted((kd, g, acc) for g, (acc, kd) in by_mol.get(mol, {}).items() if kd <= cutoff_nm)
        if mol and hits:
            rows += [(drug, g, acc, kd, "kinobeads") for kd, g, acc in hits]
            # Amendment A1: a hand target the deposit never measured, for ANY
            # drug, was absent from the profiled lysates; that is not evidence
            # it does not bind. One measured for this drug is decided by the
            # measurement.
            got = {g for _, g, _ in hits}
            unmeasured = [g for g in hand.get(drug, ()) if g not in measured and g not in got]
            rows += [(drug, g, "", "", "hand_unmeasured") for g in unmeasured]
            tally["hand_unmeasured_kept"] += len(unmeasured)
            tally["drugs_kinobeads"] += 1
        else:
            rows += [(drug, g, "", "", "hand") for g in hand.get(drug, ())]
            tally["drugs_hand"] += 1
    return rows, dict(tally)


def as_map(rows):
    m = defaultdict(list)
    for drug, gene, *_ in rows:
        m[drug].append(gene)
    return dict(m)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main(argv):
    if len(argv) not in (2, 3):
        raise SystemExit(__doc__)
    data = Path(argv[1])
    cutoff = float(argv[2]) if len(argv) == 3 else PRIMARY_CUTOFF_NM
    sys.path.insert(0, str(HERE / "analysis"))
    os.environ["DRUG_TARGETS"] = "hand"     # the hand table only; the map is what we are writing
    from decryptm_coverage import HAND_DRUG_TARGETS
    rows, tally = build(read_tsv(data / "activities.tsv"), load_targets(data / "targets.tsv", data / "reactome97_acc_gene.tsv"),
                        load_molecules(), HAND_DRUG_TARGETS, cutoff)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["drug", "gene", "uniprot", "kd_app_nM", "source"])
    for drug, gene, acc, kd, src in rows:
        w.writerow([drug, gene, acc, f"{kd:g}" if kd != "" else "", src])
    print(f"cutoff {cutoff:g} nM: {tally}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
