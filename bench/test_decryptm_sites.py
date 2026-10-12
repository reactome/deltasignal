"""specs/021 pilot: decryptM peptide sites are placed on reviewed canonical sequences."""
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import decryptm_sites as ds  # noqa: E402


def test_phospho_offsets_skip_other_modifications():
    assert ds.phospho_offsets("_(ac)AAS(ph)PT(ph)K_") == ("AASPTK", [2, 4])
    assert ds.phospho_offsets("M(ox)PEPT(ph)IDE") == ("MPEPTIDE", [4])
    assert ds.phospho_offsets("PEPTIDE") == ("PEPTIDE", [])


def test_read_fasta_keeps_only_reviewed(tmp_path):
    f = tmp_path / "p.fasta"
    f.write_text(">sp|P27361|MK03_HUMAN x\nAAAATEYVAA\n>tr|B3KR49|B3KR49_HUMAN y\nTEYVAA\n")
    assert ds.read_fasta(str(f)) == {"P27361": "AAAATEYVAA"}
    assert set(ds.read_fasta(str(f), reviewed_only=False)) == {"P27361", "B3KR49"}


def test_placed_on_reviewed_entry_not_the_trembl_leading_protein():
    # decryptM leads MAPK3's activation-loop peptide with a TrEMBL fragment; the
    # site must land on the canonical Swiss-Prot position
    seqs = {"P27361": "AAAATEYVAA"}
    rows = [{"Modified sequence": "_T(ph)EY(ph)VAA_", "Leading proteins": "B3KR49",
             "Proteins": "B3KR49;P27361-2;P27361", "Experiment": "e", "Regulation": "down"}]
    stats = collections.Counter()
    out = list(ds.sites(rows, seqs, stats))
    assert [(r["uniprot"], r["position"], r["residue"]) for r in out] == [("P27361", 5, "T"), ("P27361", 7, "Y")]
    assert out[0]["regulation"] == "down" and out[0]["n_sites_on_peptide"] == 2


def test_ambiguous_or_unplaceable_peptides_are_skipped_and_counted():
    seqs = {"P1": "TEYTEY", "P2": "AAAA"}
    stats = collections.Counter()
    rows = [{"Modified sequence": "T(ph)EY", "Leading proteins": "P1", "Proteins": "P1", "Experiment": "e"},
            {"Modified sequence": "T(ph)EY", "Leading proteins": "P2", "Proteins": "P2", "Experiment": "e"}]
    assert list(ds.sites(rows, seqs, stats)) == []
    assert stats["no reviewed sequence places it uniquely"] == 2
    # a blank regulation means not regulated
    rows = [{"Modified sequence": "A(ph)A", "Leading proteins": "P3", "Proteins": "P3", "Experiment": "e",
             "Regulation": ""}]
    assert list(ds.sites(rows, {"P3": "CAAC"}, stats))[0]["regulation"] == "not"
