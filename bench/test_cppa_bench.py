import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_bench as c  # noqa: E402


def test_drug_genes_translates_accessions_and_drops_stimuli(tmp_path):
    (tmp_path / "t.tsv").write_text("drug\trole\tchembl_parent\ttarget\tkd_app_nM\tsource\n"
                                    "A\ttest\tCHEMBL1\tP1\t5\tkinobeads\nA\ttest\tCHEMBL1\tP9\t5\tkinobeads\n"
                                    "MEKINHIBITOR01\ttrain\t\tMAP2K1\t\tclass\nEGF\tstimulus\t\tEGF\t\tligand\n")
    (tmp_path / "g.tsv").write_text("acc\tgene\nP1\tBRAF\n")
    assert c.drug_genes(tmp_path / "t.tsv", tmp_path / "g.tsv") == {"A": ("test", ["BRAF"]),
                                                                       "MEKINHIBITOR01": ("train", ["MAP2K1"])}
