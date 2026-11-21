"""specs/169: split, labels, row choice and case assembly for the genetic-panel comparison."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "169"))
import common_169 as c  # noqa: E402
import truth_169 as t  # noqa: E402
import score_169 as s  # noqa: E402


def test_split_is_fixed_by_symbol_and_roughly_even():
    genes = [f"G{i}" for i in range(2000)]
    a = sum(c.split_of(g) == "A" for g in genes)
    assert 900 < a < 1100 and c.split_of("TP53") == c.split_of("TP53")


def test_gene_of_parses_screen_ids():
    assert c.gene_of("10000_ZBTB46_P1P2_ENSG00000130584") == "ZBTB46"
    assert c.gene_of("non-targeting") is None


def test_labels_at_both_thresholds():
    assert t.label(0.6, 2 ** -0.585, 2 ** 0.585) == "down" and t.label(0.7, 2 ** -0.585, 2 ** 0.585) == "not"
    assert t.label(0.84, 0.85, 1.15) == "down" and t.label(1.15, 0.85, 1.15) == "up"


def test_choose_rows_requires_engagement_and_takes_most_cells():
    ids = ["1_A_P1_E1", "2_A_P2_E1", "3_B_P1_E2", "4_C_P1_E3", "ctrl"]
    kd = [0.3, 0.2, 0.8, 0.1, 1.0]
    cells = [40, 90, 100, 10, 500]
    ctrl = [False, False, False, False, True]
    assert t.choose_rows({"A", "B", "C"}, kd, cells, ids, ctrl) == {"A": 1}   # B not knocked down, C too few cells


def test_cases_read_reach_at_record_precision_and_matched_calls():
    preds = [{"pathway": "P", "gene": "A", "target": "X", "fold": "1.0000000000000218", "call": "not"},
             {"pathway": "P", "gene": "A", "target": "Y", "fold": "0.7", "call": "down"},
             {"pathway": "P", "gene": "A", "target": "Z", "fold": "0.5", "call": "down"}]
    truth = {("A", "X"): {"measured_50": "not"}, ("A", "Y"): {"measured_50": "down"}}
    cs = s.cases(preds, truth, "measured_50")
    assert cs[("P", "A", "X")] == ("not", "not", 1.0) and ("P", "A", "Z") not in cs
    m = s.cases(preds, truth, "measured_50", (2 ** -0.585, 2 ** 0.585))
    assert m[("P", "A", "Y")][0] == "not"            # 0.7 is not a 1.5x change


def test_h5_col_reads_both_categorical_layouts(tmp_path):
    h5py = pytest.importorskip("h5py")
    import numpy as np
    p = tmp_path / "x.h5ad"
    with h5py.File(p, "w") as f:
        obs = f.create_group("obs")
        g = obs.create_group("new_cat")
        g["categories"] = np.array([b"a", b"b"])
        g["codes"] = np.array([1, 0, 1])
        var = f.create_group("var")
        var["gene_name"] = np.array([1, 0])
        var.create_group("__categories")["gene_name"] = np.array([b"TP53", b"MDM2"])
    with h5py.File(p, "r") as f:
        assert list(c.h5_col(f, "obs", "new_cat")) == ["b", "a", "b"]
        assert list(c.h5_col(f, "var", "gene_name")) == ["MDM2", "TP53"]
