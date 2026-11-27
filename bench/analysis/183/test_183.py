import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pakt_183 as P  # noqa: E402

CALLS = {("K1", "p-AKT"): "down", ("K1", "p-GSK3b"): "down", ("K1", "p-p38"): "up",
         ("K2", "p-AKT"): "not", ("K2", "p-GSK3b"): "not", ("K2", "p-p38"): "not",
         ("K3", "p-AKT"): "up", ("K3", "p-GSK3b"): "not",
         ("K4", "p-GSK3b"): "down"}                      # not measured at p-AKT: in neither group


def test_co_movement_groups_and_rates():
    c = P.co_movement(P.by_protein(CALLS))
    assert (c["proteins_lowering"], c["proteins_rest"]) == (1, 2)
    g = c["readouts"]["p-GSK3b"]
    assert g["lowers"] == {"n": 1, "up": 0.0, "down": 1.0} and g["rest"] == {"n": 2, "up": 0.0, "down": 0.0}
    assert c["readouts"]["p-p38"]["rest"]["n"] == 1


def test_activity_counts_other_readouts_only():
    a = P.activity(P.by_protein(CALLS))
    assert a["1-3"] == {"proteins": 1, "lower": 1.0, "raise": 0.0}     # K1: two other readouts changed
    assert a["0-0"] == {"proteins": 2, "lower": 0.0, "raise": 0.5}     # K2 and K3 change nothing else
    assert [P.band(n) for n in (0, 3, 4, 9, 10, 30)] == ["0-0", "1-3", "4-9", "4-9", "10+", "10+"]
