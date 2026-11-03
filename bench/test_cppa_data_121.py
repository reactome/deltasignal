import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import cppa_data_121 as d  # noqa: E402


def test_portions():
    assert d.portion_of("train", "train") == "train"
    assert d.portion_of("seen", "train") == "train" and d.portion_of("stimulus", "train") == "train"
    assert d.portion_of("test", "train") == "A"
    assert d.portion_of("seen", "test") == "B"
    assert d.portion_of("test", "test") == "C"


def test_a_test_portion_needs_the_final_flag(tmp_path):
    for p in ("A", "B", "C"):
        with pytest.raises(SystemExit):
            d.load(tmp_path, tmp_path / "split.tsv", p)


def test_load_returns_only_the_portion_asked_for(tmp_path):
    ts = tmp_path / "delta-p_TreatmentSummary"
    rp = tmp_path / "delta-p_RPPA"
    ts.mkdir(), rp.mkdir()
    cols = "UID\tset\tresource\tcell_line_name\ttime\t2D_3D\tstimuli\tcompound_name_1\tdosage_1\tcompound_name_2\tdosage_2"
    (ts / "x_set7.tsv").write_text(cols + "\nu1\tset7\t\tC1\t24hr\t\t\tDRUGA\t1\t\t\nu2\tset7\t\tC1\t24hr\t\t\tDRUGB\t1\t\t\n")
    (rp / "x_delta-p_RPPA_set7.tsv").write_text("UID\tAkt\nu1\t0.5\nu2\t-9\n")
    (tmp_path / "split.tsv").write_text("unit\tname\trole\ndrug\tDRUGA\ttrain\ndrug\tDRUGB\ttest\ncell_line\tC1\ttrain\n")
    tr = d.load(tmp_path, tmp_path / "split.tsv", "train")
    assert [(x["uid"], x["values"]) for x in tr] == [("u1", {"Akt": 0.5})]
    a = d.load(tmp_path, tmp_path / "split.tsv", "A", final=True)
    assert [x["uid"] for x in a] == ["u2"]
