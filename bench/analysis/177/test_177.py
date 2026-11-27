"""specs/177: the rule, the half reader's refusal and the class tables, on synthetic inputs and the shipped tables."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_177 as S  # noqa: E402
from common_173 import split_of  # noqa: E402

HEADER = "experiment,marker,perturbation,row_col,timepoint,origin,target,bin,counts,fit_var,median,overall_var,var_ratio\n"


def gene_in(half):
    return next(g for g in (f"G{i}" for i in range(1000)) if split_of(g) == half)


def test_rule_adds_up_only_where_the_model_is_silent():
    tk, tyr = {"SRC"}, {"p-STAT5"}
    assert S.apply_rule("not", "SRC", "p-STAT5", tk, tyr) == "up"
    assert S.apply_rule("down", "SRC", "p-STAT5", tk, tyr) == "down"      # the model's own call stands
    assert S.apply_rule("not", "SRC", "p-AKT", tk, tyr) == "not"
    assert S.apply_rule("not", "AKT1", "p-STAT5", tk, tyr) == "not"
    assert S.apply_rule("up", "AKT1", "p-AKT", tk, tyr) == "up"


def test_half_reader_refuses_the_other_half(tmp_path):
    a, b = gene_in("A"), gene_in("B")
    f = tmp_path / "x.csv"
    f.write_text(HEADER + "E1,FLAG-GFP1,EGF,A1,0,GFP,p-AKT,0,50,0,1.0,0,0\nE1,MB,EGF,A2,0,GFP,p-AKT,0,50,0,1.0,0,0\n")
    assert len(S.load_half(f, {"MB": b}, "B")) == 2
    with pytest.raises(SystemExit):
        S.load_half(f, {"MB": b}, "A")
    with pytest.raises(SystemExit):
        S.load_half(f, {"MA": a}, "B")                                    # MB is unmapped here


def test_sign_test_and_cases():
    assert S.sign_test(5, 0) == pytest.approx(1 / 32)
    m = {("SRC", "p-STAT5"): (0.5, 0.2), ("SRC", "p-AKT"): (-0.5, 0.2), ("AKT1", "p-STAT5"): (0.1, 0.2)}
    pred = {("SRC", "p-STAT5"): ("not", "not"), ("SRC", "p-AKT"): ("up", "up"), ("AKT1", "p-STAT5"): ("conflict", "up")}
    cs = {(c["gene"], c["readout"]): c for c in S.build_cases(m, pred, {"SRC"}, {"p-STAT5"})}
    assert cs[("SRC", "p-STAT5")]["candidate"] == "up" and cs[("SRC", "p-STAT5")]["model"] == "not"
    assert cs[("SRC", "p-AKT")]["candidate"] == "up" and cs[("SRC", "p-AKT")]["measured"] == "down"
    assert cs[("AKT1", "p-STAT5")]["model"] == "not"                      # a conflict is no call
    assert S.build_cases(m, pred, {"SRC"}, {"p-STAT5"}, mult=3)[1]["measured"] == "not"   # 0.5 < 3 x 0.2


def test_shipped_tables():
    tk, tyr = S.read_classes()
    assert tyr == {"p-BTK/ITK", "p-ERK1/2", "p-FAK", "p-JNK", "p-PLCg2", "p-SRC", "p-STAT3", "p-STAT5", "p-p38"}
    assert {"SRC", "EGFR", "JAK3", "FGFR1"} <= tk and len(tk) == 60
    assert not tk & {"MAP2K1", "MAP2K6", "AKT1", "BRAF", "EPHB6"}          # dual-specificity and non-kinases are out
