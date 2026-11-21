"""Issue #278: the korkut parity gate's comparison, and the readout mode the ensemble sends per kind."""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "analysis"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


par = _load("korkut_parity_157", HERE / "analysis" / "157" / "korkut_parity_157.py")
ens = _load("ensemble_157", HERE / "analysis" / "157" / "ensemble_157.py")


def test_noise_below_print_precision_is_equal_and_a_real_difference_is_not():
    kinds = {"A": "phospho", "T": "total"}
    record = {("P", "SR", "A"): (0.538, "down"), ("P", "SR", "T"): (1.0, "not"),
              ("P", "AK", "A"): (0.105, "down")}
    mean = (0.00538 * 3) / 3 / 0.01                  # computed, as the engine's mean is
    replay = {("P", "SR", "A"): (mean, "down"), ("P", "SR", "T"): (1.0000000000000218, "not"),
              ("P", "AK", "A"): (0.006, "down")}     # the forms-mode value issue #278 found
    c, diffs = par.compare(replay, record, kinds)
    assert c["phospho_total"] == 2 and c["phospho_equal"] == 1
    assert c["total_total"] == 1 and c["total_equal"] == 1
    assert [d[2] for d in diffs] == ["fold"]


def test_missing_and_call_differences_fail():
    record = {("P", "SR", "A"): (0.5, "down"), ("P", "SR", "B"): (2.0, "up")}
    replay = {("P", "SR", "A"): (0.5, "not")}
    _, diffs = par.compare(replay, record, {"A": "phospho", "B": "phospho"})
    assert sorted(d[2] for d in diffs) == ["call", "missing from replay"]


def test_ensemble_sends_phospho_as_sites_and_totals_as_abundance():
    assert ens.korkut_query_mode("phospho", "sites") == "sites"
    assert ens.korkut_query_mode("total", "sites") == "abundance"
    assert ens.korkut_query_mode("phospho", "forms") == "forms"   # the 2026-11-15 runs
    assert ens.korkut_query_mode("total", "forms") == "abundance"


mrg = _load("merge_korkut_278", HERE / "analysis" / "157" / "merge_korkut_278.py")


def _write(path, rows):
    import gzip
    with gzip.open(path, "wt") as f:
        f.write("\t".join(mrg.HEADER) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_merge_replaces_korkut_rows_and_keeps_the_rest(tmp_path, monkeypatch):
    old, new, out = tmp_path / "old", tmp_path / "new", tmp_path / "out"
    old.mkdir(); new.mkdir()
    cur = ["curator", "P", "G", "0", "1", "0.5", "down", "1"]
    _write(old / "calls_canonical.tsv.gz", [cur, ["korkut", "P", "SR", "", "A", "1.0", "not", "1"]])
    _write(new / "calls_canonical.tsv.gz", [["korkut", "P", "SR", "", "A", "0.538", "down", "1"]])
    monkeypatch.setattr(sys, "argv", ["m", str(old), str(new), str(out)])
    mrg.main()
    got = (out / "calls_canonical.tsv").read_text().splitlines()
    assert got[1].split("\t") == cur
    assert got[2].split("\t")[5:7] == ["0.538", "down"]


def test_merge_refuses_a_different_korkut_population(tmp_path, monkeypatch):
    import pytest
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir(); new.mkdir()
    _write(old / "calls_canonical.tsv.gz", [["korkut", "P", "SR", "", "A", "1.0", "not", "1"]])
    _write(new / "calls_canonical.tsv.gz", [["korkut", "P", "SR", "", "B", "1.0", "not", "1"]])
    monkeypatch.setattr(sys, "argv", ["m", str(old), str(new), str(tmp_path / "out")])
    with pytest.raises(SystemExit):
        mrg.main()
