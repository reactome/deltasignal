import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import split_121 as s  # noqa: E402


def test_norm_merges_synonyms_and_spelling():
    assert s.norm("Gleevec") == s.norm("imatinib")
    assert s.norm("PD-0325901") == s.norm("PD325901")
    assert s.norm("I-BET726") != s.norm("I-BET762")


def test_seen_and_stimuli_never_test_and_split_is_seeded():
    pairs = [(f"D{i}", f"C{j}") for i in range(20) for j in range(8)] + [("EGF", "C0"), ("SEEN", "C1")]
    rows = s.split(pairs, {"SEEN"})
    role = {(u, n): r for u, n, r in rows}
    assert role[("drug", "SEEN")] == "seen" and role[("drug", "EGF")] == "stimulus"
    assert sum(r == "test" for (u, _), r in role.items() if u == "drug") == 8        # 40% of the 20 free drugs
    assert sum(r == "test" for (u, _), r in role.items() if u == "cell_line") == 2   # 25% of 8 cell lines
    assert rows == s.split(list(reversed(pairs)), {"SEEN"})                          # input order does not matter
