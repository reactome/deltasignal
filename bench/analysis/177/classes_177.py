#!/usr/bin/env python3
"""specs/177: the two class tables the tyrosine-kinase rule reads. Identifiers only; no Lun value is read.

  tyrosine_kinases.tsv  every Lun protein of interest (both halves) with Reactome's own UniProt keywords for its
                        reference entity: tyrosine-specific = 'Tyrosine-protein kinase' without
                        'Serine/threonine-protein kinase' (dual-specificity kinases such as the MAP2Ks are left out)
  readout_sites.tsv     each phospho readout's residues at its registered sites, from the canonical sequences

Usage: classes_177.py FASTA OUT_DIR     (Neo4j: NEO4J_URL, the release the predictions were built from)
"""
import csv
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1] / "173"))
sys.path.insert(0, str(HERE.parents[1] / "174"))
sys.path.insert(0, str(HERE.parents[1] / "176"))
sys.path.insert(0, str(HERE.parents[2]))
import lun_a_173 as A  # noqa: E402
from axis_176 import readout_residues  # noqa: E402

TYR, ST = "Tyrosine-protein kinase", "Serine/threonine-protein kinase"
KEYWORDS = """MATCH (re:ReferenceGeneProduct) WHERE re.identifier IN $accs AND re.variantIdentifier IS NULL
RETURN re.identifier AS acc, re.keyword AS kw"""


def main():
    from py2neo import Graph
    fasta, out = sys.argv[1], Path(sys.argv[2])
    g = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"))
    release = g.run("MATCH (d:DBInfo) RETURN d.releaseNumber AS r").evaluate()
    pois = A.read_tsv(A.MAPS / "pois.tsv")
    kw = {r["acc"]: set(r["kw"] or ()) for r in g.run(KEYWORDS, accs=sorted({r["accession"] for r in pois if r["accession"]}))}
    with open(out / "tyrosine_kinases.tsv", "w", newline="") as f:
        f.write(f"# specs/177: UniProt keywords as Reactome release {release} stores them on the reference entity. No Lun value read.\n")
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "accession", "in_reactome", "tyrosine_kinase_keyword", "ser_thr_kinase_keyword", "tyrosine_specific"])
        seen = set()
        for r in pois:
            if r["gene"] in seen:
                continue
            seen.add(r["gene"])
            k = kw.get(r["accession"])
            w.writerow([r["gene"], r["accession"], int(k is not None), int(TYR in (k or ())), int(ST in (k or ())),
                        int(TYR in (k or ()) and ST not in (k or ()))])
    with open(out / "readout_sites.tsv", "w", newline="") as f:
        f.write("# specs/177: residues at each phospho readout's registered sites (bench/analysis/127/readouts.tsv), canonical\n"
                "# UniProt sequences; '?' = the coordinate is not S/T/Y in the canonical sequence (isoform numbering).\n")
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["readout", "residues", "has_tyrosine"])
        for name, res in sorted(readout_residues(fasta).items()):
            w.writerow([name, res, int("Y" in res)])


if __name__ == "__main__":
    main()
