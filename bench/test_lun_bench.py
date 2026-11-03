"""specs/127: the Lun 2019 design tables as the prediction bench reads them."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lun_bench as lb  # noqa: E402
from korkut_bench import read_tsv  # noqa: E402


def test_proteins_are_genes_with_an_accession_once_each():
    p = lb.pois()
    genes = [g for g, _ in p]
    assert len(genes) == len(set(genes))
    rows = read_tsv(lb.MAPS / "pois.tsv")
    assert len(rows) == 649
    assert set(genes) == {r["gene"] for r in rows if r["accession"]}
    assert all(a and "-" not in a for _, a in p)


def test_readout_panel_matches_specs_127():
    rows = read_tsv(lb.MAPS / "readouts.tsv")
    scored = [r for r in rows if r["kind"] in ("phospho", "total")]
    assert len(rows) == 36 and len(scored) == 32
    assert {r["readout"] for r in rows if r["kind"] == "exclude"} == {"GAPDH", "GFP", "p-HH3", "cleaved PARP"}
    for r in scored:
        for s in r["sites"].split(","):
            if r["kind"] == "phospho":
                g, c = s.split(":")
                assert g and int(c) > 0
    by = {r["readout"]: r for r in rows}
    assert lb.readout_genes(by["p-MKK3/6"]) == {"MAP2K3", "MAP2K6"}
    assert "MAP2K3:218" in by["p-MKK3/6"]["sites"]
    assert lb.readout_genes(by["p-STAT5"]) == {"STAT5A"}
