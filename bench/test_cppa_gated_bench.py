import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_gated_bench as b  # noqa: E402


def test_model_for_exact_then_longest_parent():
    models = {"A375": "ACH-1", "A37": "ACH-0", "MCF7": "ACH-2"}
    assert b.model_for("MCF7", models) == ("ACH-2", "exact")
    assert b.model_for("A375BR", models) == ("ACH-1", "parent")
    assert b.model_for("BAF3", models) == (None, "none")


def test_depmap_models_reads_names_only(tmp_path):
    f = tmp_path / "m.csv"
    f.write_text("ModelID,PatientID,CellLineName,StrippedCellLineName,Age,Sex\nACH-1,P1,A-375,A375,50,X\n")
    assert b.depmap_models(f) == {"A375": "ACH-1"}


def test_expression_rows(tmp_path):
    f = tmp_path / "e.csv"
    f.write_text(",AKT1 (207),EGFR (1956)\nACH-1,5.0,0.2\nACH-2,1,1\n")
    assert b.expression_rows(f, {"ACH-1"}, {"EGFR"}) == {"ACH-1": {"EGFR": 0.2}}
