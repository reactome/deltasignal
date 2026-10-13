"""specs/050 R1: the competition census counts the motif it claims to.

Fixture: S feeds two reactions. R1 (partner G, the perturbed gene) makes an
inert complex; R2 (partner Q) makes the readout O. A G knockout should be
counted as 'competition agrees' when O is expected UP, and dropped when the
GO filter excludes R1 or when R1's branch acts on something.
"""
import csv, os, subprocess, sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "analysis" / "competition_census.py"


def write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def build(tmp, acts=False):
    b = tmp / "build"
    d = b / "R-HSA-1"
    d.mkdir(parents=True)
    (b / "pathways.tsv").write_text("id\tpathway_name\nR-HSA-1\tP\n")
    edges = [("S", "r1", "pos", "and", "input", 1, "R-1"), ("G", "r1", "pos", "and", "input", 1, "R-1"),
             ("r1", "C", "pos", "", "output", 1, "R-1"),
             ("S", "r2", "pos", "and", "input", 1, "R-2"), ("Q", "r2", "pos", "and", "input", 1, "R-2"),
             ("r2", "O", "pos", "", "output", 1, "R-2")]
    if acts:
        edges.append(("C", "r3", "pos", "and", "catalyst", 1, "R-3"))
    write(d / "logic_network.csv", ["source_id", "target_id", "pos_neg", "and_or", "edge_type",
                                    "stoichiometry", "edge_reaction_id"], edges)
    write(d / "stid_to_uuid_mapping.csv", ["uuid", "stable_id"],
          [(u, f"R-HSA-{u}") for u in ("S", "G", "Q", "C", "O", "r1", "r2", "r3")])
    cases = tmp / "cases.tsv"
    cols = ["pathway", "gene", "direction", "key_output", "predicted", "expected", "valid",
            "gene_uuids", "output_uuids"]
    with open(cases, "w") as f:
        f.write("\t".join(cols) + "\n")
        f.write("\t".join(["P", "GENE", "0", "o", "1", "2", "1", "G", "O"]) + "\n")
    return b, cases


def run(b, cases, **env):
    e = {**os.environ, **env}
    return subprocess.run([sys.executable, str(SCRIPT), str(b), str(cases)], env=e,
                          capture_output=True, text=True, check=True).stdout


def test_knockout_of_a_competitor_partner_is_counted(tmp_path):
    out = run(*build(tmp_path))
    assert "1  now wrong, competition agrees" in out


def test_go_filter_excludes_unannotated_reactions(tmp_path):
    b, cases = build(tmp_path)
    neg = tmp_path / "negreg.csv"
    neg.write_text('r.stId, kind, name\n"R-2", "reaction", "negative regulation of x"\n')
    assert "competition agrees" not in run(b, cases, NEGREG=str(neg))
    neg.write_text('r.stId, kind, name\n"R-1", "reaction", "negative regulation of x"\n')
    assert "1  now wrong, competition agrees" in run(b, cases, NEGREG=str(neg))


def test_inert_sink_filter_drops_an_acting_branch(tmp_path):
    b, cases = build(tmp_path, acts=True)
    assert "1  now wrong, competition agrees" in run(b, cases)   # counted without the filter
    assert "competition agrees" not in run(b, cases, SINK="1")
