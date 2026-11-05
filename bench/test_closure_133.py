"""specs/133: closure_133.py's choice of source pathway is deterministic (review 2026-11-05)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "analysis"))
import closure_133 as cl  # noqa: E402


def test_pick_diagram_takes_least_depth_then_smallest_stid_and_reports_ties():
    rows = [{"s": "R-S1", "q_diagram": "R-HSA-9", "name": "nine", "L": 1},
            {"s": "R-S1", "q_diagram": "R-HSA-2", "name": "two", "L": 1},
            {"s": "R-S1", "q_diagram": "R-HSA-1", "name": "one", "L": 2},
            {"s": "R-S2", "q_diagram": "R-HSA-5", "name": "five", "L": 3},
            {"s": "R-S2", "q_diagram": "R-HSA-3", "name": "three", "L": 1}]
    for order in (rows, rows[::-1]):              # row order (Neo4j's) never decides
        best, names, ties = cl.pick_diagram(order)
        assert best == {"R-S1": "R-HSA-2", "R-S2": "R-HSA-3"}
        assert ties == {"R-S1": ["R-HSA-2", "R-HSA-9"]}
        assert names["R-HSA-9"] == "nine"
