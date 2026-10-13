"""specs/051: the rewired-null comparison, against hand-computed values."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import rewire_null as rn  # noqa: E402

HEAD = "pathway\tgene\tdirection\tkey_output\tpredicted\texpected\n"


def _write(d, axis, rows):
    d.mkdir(exist_ok=True)
    (d / f"{axis}_cases.tsv").write_text(HEAD + "".join("\t".join(r) + "\n" for r in rows))


def test_pairing_slices_and_rule(tmp_path):
    held, tun = "Some_Pathway", "Signaling_by_WNT"
    # wire0 right on 4 held-out cases over 2 genes; seeds right on 1 of them
    real = [(held, g, "0", str(i), "0", "0") for g, i in (("A", 1), ("A", 2), ("B", 3), ("B", 4))]
    real += [(tun, "C", "0", "9", "2", "2")]
    s1 = [r[:4] + ("1", r[5]) for r in real]
    s1[0] = real[0]
    s2 = [r[:4] + ("1", r[5]) for r in real][:-1]      # loses the tuning case: dropped from pairing
    s2[0] = real[0]
    for d, rows in (("w0", real), ("s1", s1), ("s2", s2)):
        _write(tmp_path / d, "curator", rows)
        _write(tmp_path / d, "experimental", rows[:2])
    res = rn.run(tmp_path / "w0", [tmp_path / "s1", tmp_path / "s2"], None)
    cur = res["curator"]
    assert cur["unpaired"] == {"wire0": 5, "seeds": [5, 4], "paired": 4}
    assert "n" not in cur["tuning"] or cur["tuning"]["n"] == 0
    h = cur["held-out"]
    assert (h["n"], h["genes"]) == (4, 2)
    assert h["wire0"]["accuracy"] == 1.0 and h["seed_mean"]["accuracy"] == 0.25
    # diff per case: 0, 1, 1, 1 -> 0.75
    assert abs(h["correct_diff_ci95_by_gene"][0] - 0.75) < 1e-12
    # gene A diffs (0,1) mean .5, gene B (1,1) mean 1: bootstrap draws lie in [.5, 1]
    assert 0.5 <= h["correct_diff_ci95_by_gene"][1] <= h["correct_diff_ci95_by_gene"][2] <= 1.0
    assert h["carries_signal"] is True
    # experimental: one gene only -> interval collapses to the point
    e = res["experimental"]["all"]
    assert e["n"] == 2 and e["correct_diff_ci95_by_gene"][1] == e["correct_diff_ci95_by_gene"][0] == 0.5


def test_rule_fails_when_a_seed_matches(tmp_path):
    rows = [("P", g, "0", "1", "0", "0") for g in "ABCD"]
    for d in ("w0", "s1"):
        _write(tmp_path / d, "curator", rows)
        _write(tmp_path / d, "experimental", rows)
    res = rn.run(tmp_path / "w0", [tmp_path / "s1"], tmp_path / "w0")
    m = res["experimental"]["all"]
    assert m["carries_signal"] is False and m["correct_diff_ci95_by_gene"] == [0.0, 0.0, 0.0]
    assert m["canonical"]["accuracy"] == 1.0
