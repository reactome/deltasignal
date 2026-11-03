#!/usr/bin/env python3
"""specs/121 phase 1: CPPA drug -> target proteins, by specs/113's rule (bench/p100_targets.py).

A named compound is found in ChEMBL 37 by name, then synonym (CPPA gives no structure), and taken to its parent.
Targets: kinobeads at Kd_app <= 1 uM, else ChEMBL drug_mechanism, else ChEMBL binding (pChEMBL >= 6), as in specs/113.
An anonymised compound ("MEK inhibitor-01") takes the genes of the class it names (CLASSES below), source "class".
Stimuli (ligands) take their own gene, source "ligand". Reads only the CPPA treatment summaries, no measurement.

CPPA's own spelling is looked up first; a few misspellings are corrected (SPELLING) before the ChEMBL search.

Usage: cppa_targets_121.py SPLIT_TSV CPPA_UNPACKED KLAEGER_DIR CACHE_DIR > map.tsv   (the map stays outside the repository)
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import p100_targets as pt  # noqa: E402
import split_121  # noqa: E402

# the target class an anonymised CPPA compound names -> HGNC genes (family members a class inhibitor is made for)
CLASSES = {
    "AKT": ["AKT1", "AKT2", "AKT3"], "AKT12": ["AKT1", "AKT2"], "ATR": ["ATR"], "BET": ["BRD2", "BRD3", "BRD4"],
    "BRAF": ["BRAF"], "CHK1": ["CHEK1"], "EGFR": ["EGFR"], "FAK": ["PTK2"], "FAS": ["FASN"], "GSK3": ["GSK3A", "GSK3B"],
    "HDAC": ["HDAC1", "HDAC2", "HDAC3"], "KDM5": ["KDM5A", "KDM5B", "KDM5C", "KDM5D"], "MEK": ["MAP2K1", "MAP2K2"],
    "MET": ["MET"], "MTOR": ["MTOR"], "NFKB": ["IKBKB"], "PARP": ["PARP1", "PARP2"], "PDK1": ["PDPK1"],
    "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"], "WEE1": ["WEE1"], "IGFR": ["IGF1R"],
}
# CPPA spellings that ChEMBL does not know, and the compound meant
SPELLING = {"METROTREXATE": "METHOTREXATE", "AMITRIPTYLENE": "AMITRIPTYLINE", "PODOPHYLLOTOCIN": "PODOPHYLLOTOXIN",
            "AUY": "LUMINESPIB", "ZVAD": "Z-VAD-FMK", "3-MA": "3-METHYLADENINE", "2-DG": "2-DEOXY-D-GLUCOSE",
            "SPIROMYCIN": "SPIRAMYCIN", "MEVINOLIN": "LOVASTATIN"}
LIGANDS = {"EGF": ["EGF"], "FGF1": ["FGF1"], "FGF2": ["FGF2"], "HGF": ["HGF"], "IGF1": ["IGF1"], "INSULIN": ["INS"],
           "NRG1": ["NRG1"]}


def anonymised_class(name):
    """'MEKINHIBITOR06' -> 'MEK'; 'AKT12INHIBITOR01' -> 'AKT12'; else None."""
    m = re.match(r"^([A-Z0-9]+?)INHIBITOR\d+$", name)
    return m.group(1) if m else None


def spellings(name):
    """CPPA's spelling, then code-number forms ChEMBL lists (PD0325901 -> PD-0325901, PLX4720 -> PLX-4720)."""
    out = [name]
    hyphen = re.sub(r"^([A-Za-z]+)-?(\d)", r"\1-\2", name)
    for v in (hyphen, hyphen.replace("-", " "), name.replace("-", "")):
        if v not in out:
            out.append(v)
    return out


def originals(unpacked):
    """{normalised name: CPPA's spelling} from the treatment summaries."""
    import glob
    out = {}
    for f in sorted(glob.glob(f"{unpacked}/delta-p_TreatmentSummary/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t"):
            for i in range(1, 5):
                x = (r.get(f"compound_name_{i}") or "").strip()
                if x and x.upper() != "NA":
                    out.setdefault(split_121.norm(x), x)
    return out


def build(split_tsv, unpacked, klaeger_dir, cache):
    Path(cache).mkdir(parents=True, exist_ok=True)
    orig = originals(unpacked)
    kb, kb_mols = pt.kinobeads(klaeger_dir, 1000.0)
    drugs = [r for r in csv.DictReader(open(split_tsv), delimiter="\t") if r["unit"] == "drug"]
    rows, log = [], []
    for d in drugs:
        name, role = d["name"], d["role"]
        cls = "IGFR" if name == "IGFR" else anonymised_class(name)     # CPPA's "IGFR" is an IGF1R inhibitor
        if name in LIGANDS:
            rows += [(name, role, "", g, "", "ligand") for g in LIGANDS[name]]
            log.append((name, role, "ligand", len(LIGANDS[name])))
            continue
        if cls in CLASSES:
            rows += [(name, role, "", g, "", "class") for g in CLASSES[cls]]
            log.append((name, role, "class", len(CLASSES[cls])))
            continue
        if role == "stimulus":
            log.append((name, role, "none", 0))
            continue
        mol, accs, src = None, {}, "none"
        spelled = orig.get(name, name)
        for variant in spellings(SPELLING.get(spelled.upper(), spelled)):
            for mol, how in pt.find_molecules(cache, "", variant):
                accs, src = pt.targets_of(cache, mol, kb, kb_mols)
                if accs:
                    break
            if accs:
                break
        rows += [(name, role, mol, a, accs[a], src) for a in sorted(accs)]
        log.append((name, role, src, len(accs)))
    return rows, log


def main(argv):
    if len(argv) != 5:
        raise SystemExit(__doc__)
    rows, log = build(*argv[1:5])
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["drug", "role", "chembl_parent", "target", "kd_app_nM", "source"])
    w.writerows(rows)
    for r in log:
        print("\t".join(map(str, r)), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
