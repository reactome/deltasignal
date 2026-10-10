"""A pathway whose parse or solve FAILS must fail the benchmark run (code
review 2026-10-02). It used to print "[skip] solve_failed", drop the pathway's
cases from every metric and exit 0, so catalog.sh and run_arm.sh recorded a
partial run as complete. Runtime test: main() itself is driven with one pathway
scoring and one failing, so the guard is exercised, not just present."""
import sys
from collections import Counter
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / "analysis"))

import benchmark_vs_mpbiopath as bench  # noqa: E402


def ok_result(name):
    return {"status": "ok", "name": name, "id": name, "total": 2, "correct": 2,
            "accuracy": 1.0, "valid_total": 2, "valid_correct": 2, "valid_accuracy": 1.0,
            "confusion": {("UP", "UP"): 1, ("DOWN", "DOWN"): 1},
            "failure_categories": {}, "n_perturbations": 1, "n_key_outputs": 2, "case_log": []}


def run_main(tmp_path, monkeypatch, outcomes):
    plist = tmp_path / "pathways.tsv"
    plist.write_text("pathway_id\tpathway_name\n" + "".join(f"{n}\t{n}\n" for n in outcomes))
    monkeypatch.setattr(bench, "apply_catalog_ids", lambda p, _: p)
    monkeypatch.setattr(bench, "network_edge_count", lambda pid: 10)
    monkeypatch.setattr(bench, "run_pathway",
                        lambda pid, pname, **kw: outcomes[pname](pname))
    monkeypatch.setattr(bench, "SOLVE_TALLY", Counter())
    monkeypatch.setattr(sys, "argv", ["bench", "--pathway-list", str(plist),
                                      "--report", str(tmp_path / "report.tsv")])
    bench.main()


@pytest.mark.parametrize("status", sorted(bench.FATAL_STATUSES))
def test_a_failed_pathway_fails_the_run(tmp_path, monkeypatch, status):
    outcomes = {"PathA": ok_result,
                "PathB": lambda n: {"status": status, "name": n, "error": "HTTP 500"}}
    with pytest.raises(SystemExit) as e:
        run_main(tmp_path, monkeypatch, outcomes)
    assert e.value.code not in (0, None)
    assert "failed to parse or solve" in str(e.value.code)


def test_a_pathway_without_data_is_still_a_plain_skip(tmp_path, monkeypatch):
    outcomes = {"PathA": ok_result,
                "PathB": lambda n: {"status": "no_network", "name": n}}
    run_main(tmp_path, monkeypatch, outcomes)    # returns normally: not a failure


def _pathogen_outcome(listed, held):
    def outcome(name):
        bench.SOLVE_TALLY["pathogen_rule=inert"] += 1
        bench.SOLVE_TALLY["pathogens_listed"] += listed
        bench.SOLVE_TALLY["pathogens_held"] += held
        return ok_result(name)
    return outcome


def test_default_inert_on_pathways_listing_no_pathogen_is_not_refused(tmp_path, monkeypatch, capsys):
    # specs/048: inert is the default and 83 of 92 pathways list no pathogen, so a
    # subset run that holds nothing is the expected case, not "the arm is the control".
    run_main(tmp_path, monkeypatch, {"PathA": _pathogen_outcome(0, 0)})
    assert "no pathway in this run lists a pathogen" in capsys.readouterr().out


def test_inert_refused_when_listed_pathogens_match_no_node(tmp_path, monkeypatch):
    with pytest.raises(SystemExit) as e:
        run_main(tmp_path, monkeypatch, {"PathA": _pathogen_outcome(5, 0)})
    assert "matched no node" in str(e.value.code)


def test_inert_holding_listed_pathogens_runs(tmp_path, monkeypatch):
    run_main(tmp_path, monkeypatch, {"PathA": _pathogen_outcome(5, 3)})
