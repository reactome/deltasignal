"""specs/114 part 1: coverage_114.py counts on a one-pathway fixture."""
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "analysis" / "coverage_114.py"


def run(tmp_path, dataset, readout_file):
    b = tmp_path / "build"
    p = b / "R-HSA-9"
    p.mkdir(parents=True)
    # leaf PE is P1 phosphorylated at T10; leaf Q is protein P2 unmodified
    (p / "leaf_refs.csv").write_text("leaf_stid,reference_stid,modifications\n"
                                     "R-HSA-PE,uniprot:P1,\nR-HSA-Q,UniProt:P2-2,\nR-ALL-1,chebi:1,\n")
    (p / "nodes.csv").write_text("uuid,node_kind,diagram_entity_id,member_leaves\n"
                                 "u1,complex,R-HSA-C,R-HSA-PE|R-ALL-1\nu2,entity,R-HSA-Q,R-HSA-Q\n")
    ph = tmp_path / "ph.csv"
    ph.write_text('h\n"R-HSA-PE", "P1", "ERK", 10, "O-phospho-L-threonine"\n'
                  '"R-HSA-PY", "P1", "ERK", 12, "O4\'-phospho-L-tyrosine"\n'
                  '"R-HSA-QX", "P2", "MEK", 5, "glycyl phospho-5\'-adenosine"\n')
    genes = tmp_path / "genes.csv"
    genes.write_text('g, r.identifier\n"ERK", "P1"\n"MEK", "P2"\n"NOPE", "P9"\n')
    out = subprocess.run([sys.executable, str(SCRIPT), str(ph), str(genes), dataset, str(readout_file), str(b)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert out[0] == "dataset\tbuild\treadouts\tprotein_in\tsite_in_reactome\tsite_in"
    return [int(x) for x in out[1].split("\t")[2:]]


def test_decryptm(tmp_path):
    f = tmp_path / "all_sites.tsv"
    # P1 T10 twice (one site), P1 Y12 curated but not in the network, P2 S5 a non-phospho mod, P3 absent
    f.write_text("uniprot\tposition\tresidue\tregulation\n"
                 "P1\t10\tT\tdown\nP1\t10\tT\tup\nP1\t12\tY\tnot\nP2-1\t5\tS\tnot\nP3\t7\tS\tnot\n")
    assert run(tmp_path, "decryptm", f) == [4, 3, 2, 1]


def test_korkut(tmp_path):
    f = tmp_path / "readouts.tsv"
    f.write_text("# comment\nreadout\tkind\tsites\tnote\n"
                 "pERK\tphospho\tERK:12,ERK:10\t\npMEK\tphospho\tMEK:5\t\n"
                 "MEK\ttotal\tMEK\t\nX\ttotal\tNOPE\t\nY\texclude\t\tnot curated\n")
    assert run(tmp_path, "korkut", f) == [4, 3, 1, 1]


def test_p100(tmp_path):
    f = tmp_path / "analyte.txt.gz"
    with gzip.open(f, "wt") as g:
        g.write("pr_analyte_id\tpr_p100_phosphosite\tpr_uniprot_id\n"
                "a1\tS3_T10\tP1\na2\tY12\tP1\na3\tS5\tP2\na4\tS1\tP9\n")
    assert run(tmp_path, "p100", f) == [4, 3, 2, 1]
