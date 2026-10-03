"""open_network reads logic_network.csv and boundary_edges.csv as one network
(specs/044), the way the solver does."""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from network_files import open_network  # noqa: E402

HDR = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n"


def rows(d):
    return [(r["source_id"], r["target_id"], r["edge_type"]) for r in csv.DictReader(open_network(d))]


def test_pre_split_bundle_reads_as_before(tmp_path):
    (tmp_path / "logic_network.csv").write_text(HDR + "x,r,pos,and,input,1\n")
    assert rows(tmp_path) == [("x", "r", "input")]


def test_boundary_rows_follow_the_curated_ones(tmp_path):
    (tmp_path / "logic_network.csv").write_text(HDR + "x,r,pos,and,input,1\n")
    (tmp_path / "boundary_edges.csv").write_text(HDR + "a,x,pos,and,assembly,1\n")
    assert rows(tmp_path) == [("x", "r", "input"), ("a", "x", "assembly")]
    # line iteration and csv.reader see one header too
    lines = list(open_network(tmp_path))
    assert len(lines) == 3 and lines[0] == HDR


def test_no_trailing_newline_and_empty_boundary(tmp_path):
    (tmp_path / "logic_network.csv").write_text(HDR + "x,r,pos,and,input,1")
    (tmp_path / "boundary_edges.csv").write_text(HDR + "a,x,pos,and,assembly,1\n")
    assert rows(tmp_path) == [("x", "r", "input"), ("a", "x", "assembly")]
    (tmp_path / "boundary_edges.csv").write_text(HDR)
    assert rows(tmp_path) == [("x", "r", "input")]


def test_mismatched_header_is_an_error(tmp_path):
    (tmp_path / "logic_network.csv").write_text(HDR + "x,r,pos,and,input,1\n")
    (tmp_path / "boundary_edges.csv").write_text("source_id,target_id\na,x\n")
    with pytest.raises(ValueError):
        open_network(tmp_path)


def test_header_compared_as_columns(tmp_path):
    (tmp_path / "logic_network.csv").write_text(HDR + "x,r,pos,and,input,1\n")
    (tmp_path / "boundary_edges.csv").write_text(
        "\ufeff" + ",".join(f'"{c}"' for c in HDR.strip().split(",")) + "\na,x,pos,and,assembly,1\n")
    assert rows(tmp_path) == [("x", "r", "input"), ("a", "x", "assembly")]
