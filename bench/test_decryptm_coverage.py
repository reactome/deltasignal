"""specs/021 pilot: decryptm_coverage.py end to end on a one-pathway fixture."""
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "analysis" / "decryptm_coverage.py"
sys.path.insert(0, str(SCRIPT.parent))
import decryptm_coverage as dc  # noqa: E402


def test_drug_of():
    assert dc.drug_of("ddPTM_KYSE-520_SHP099_60min_R1") == ("SHP099", "KYSE-520")
    assert dc.drug_of("ddPTM_A549_Selumetinib+MK2206_1to2_60min_R1")[0] is None   # combination
    assert dc.drug_of("bad") == (None, None)


def test_coverage_counts_and_conflicts(tmp_path):
    b = tmp_path / "build"
    p = b / "R-HSA-9"
    p.mkdir(parents=True)
    (b / "pathways.tsv").write_text("id\tpathway_name\nR-HSA-9\tSome_heldout_pathway\n")
    (p / "stid_to_uuid_mapping.csv").write_text("uuid,stable_id\nu1,R-HSA-C\n")
    # complex C contains the phospho form PE (two sites) and the target PTPN11 form T
    (p / "containment.csv").write_text("stable_id,contains_stable_id,reactome_release\n"
                                       "R-HSA-C,R-HSA-PE,97\nR-HSA-C,R-HSA-T,97\n")
    ph = tmp_path / "ph.csv"
    ph.write_text('h\n"R-HSA-PE", "P1", "ERK", 10, "O-phospho-L-threonine"\n'
                  '"R-HSA-PE", "P1", "ERK", 12, "O4\'-phospho-L-tyrosine"\n'
                  '"R-HSA-PE", "P1", "ERK", 20, "glycyl phospho-5\'-adenosine"\n')
    tg = tmp_path / "tg.csv"
    tg.write_text('h\n"R-HSA-T", "PTPN11"\n')
    s = tmp_path / "s.tsv"
    s.write_text("uniprot\tposition\texperiment\tregulation\n"
                 "P1\t10\tddPTM_X_SHP099_60min_R1\tdown\n"
                 "P1\t12\tddPTM_X_SHP099_60min_R1\tdown\n"
                 "P1\t12\tddPTM_Y_SHP099_60min_R1\tup\n"
                 "P1\t10\tddPTM_X_Gefitinib_60min_R1\tup\n")      # EGFR not in this network
    out = tmp_path / "out.tsv"
    r = subprocess.run([sys.executable, str(SCRIPT), str(b), str(ph), str(tg), str(out), str(s)],
                       capture_output=True, text=True, check=True).stdout
    assert "held-out: (gene, site) pairs possible 2, measured 2 (down 1, up 0, unchanged 0, conflicting 1)" in r
    rows = out.read_text().splitlines()[1:]
    assert sorted(x.split("\t")[6] for x in rows) == ["conflict", "down"]
