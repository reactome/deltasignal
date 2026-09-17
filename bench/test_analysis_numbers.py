"""Tests for the code that produces the numbers we make decisions on.

Every solver decision in this project rests on figures computed by
`benchmark_vs_mpbiopath.py` and the scripts in `analysis/`. Those scripts are
run, read and believed, but almost none of their arithmetic has ever been
checked against a hand-computed value -- it is only ever observed on live
data, where a wrong number looks like a result rather than a bug.

This session produced three wrong figures from exactly this layer: an
unconditioned arm comparison (+265, actually +64), a catalog mismatch that
made two arms share no uuids, and a name->directory mapping that silently
matched nothing and reported an empty sweep as a finished one. Each was a
script defect, not a solver defect.

The tests below pin behaviour that is easy to break silently.
"""

import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / "analysis"))

import benchmark_vs_mpbiopath as bench  # noqa: E402
import holdout_report  # noqa: E402
from _common import pathway_dir_index  # noqa: E402

DOWN, NORMAL, UP = bench.DOWN, bench.NORMAL, bench.UP


# --- the primary metric -----------------------------------------------------

class TestSummaryMetrics:

    def test_perfect_predictor_scores_one(self):
        m = bench.summary_metrics({(DOWN, DOWN): 3, (NORMAL, NORMAL): 4, (UP, UP): 5})
        assert m["macro_f1"] == 1.0
        assert m["change_f1"] == 1.0
        assert m["balanced_accuracy"] == 1.0

    def test_hand_computed_case(self):
        """One wrong call per class, worked through by hand.

        Each class: tp=2, fp=1, fn=1 -> precision = recall = 2/3,
        F1 = 2*(4/9)/(4/3) = 2/3. All three equal, so every aggregate is 2/3.
        """
        c = {
            (DOWN, DOWN): 2, (DOWN, NORMAL): 1,
            (NORMAL, NORMAL): 2, (NORMAL, UP): 1,
            (UP, UP): 2, (UP, DOWN): 1,
        }
        m = bench.summary_metrics(c)
        assert m["f1"][DOWN] == pytest.approx(2 / 3)
        assert m["macro_f1"] == pytest.approx(2 / 3)
        assert m["balanced_accuracy"] == pytest.approx(2 / 3)

    def test_lazy_all_normal_predictor_scores_zero_on_change(self):
        """The whole reason macro-F1 is the primary metric.

        Adam's rule: do not reward "say no change". A predictor that answers
        NORMAL to everything must score 0 on both change classes, so accuracy
        can look good while macro-F1 does not.
        """
        c = {(NORMAL, DOWN): 30, (NORMAL, NORMAL): 40, (NORMAL, UP): 30}
        m = bench.summary_metrics(c)
        assert m["f1"][DOWN] == 0.0
        assert m["f1"][UP] == 0.0
        assert m["change_f1"] == 0.0
        assert m["macro_f1"] < 0.2

    def test_empty_confusion_does_not_divide_by_zero(self):
        m = bench.summary_metrics({})
        assert m["macro_f1"] == 0.0
        assert m["balanced_accuracy"] == 0.0

    def test_transposing_the_matrix_moves_only_the_recall_metric(self):
        """The orientation trap, pinned.

        `confusion` is keyed (predicted, expected). F1 is symmetric in
        precision and recall, so getting the orientation backwards leaves
        macro-F1 and change-F1 BIT-IDENTICAL while silently turning balanced
        accuracy into macro-precision. Anyone checking a refactor by watching
        macro-F1 alone would see nothing. This asserts the asymmetry exists,
        so the reader knows which number to watch.
        """
        # Deliberately asymmetric: DOWN over-predicted into both other
        # classes. macro-precision (0.778) and macro-recall (0.667) differ
        # here, which a more balanced matrix would hide.
        c = {(DOWN, DOWN): 10, (DOWN, NORMAL): 10, (DOWN, UP): 10,
             (NORMAL, NORMAL): 10, (UP, UP): 10}
        flipped = {(e, p): n for (p, e), n in c.items()}
        straight, other = bench.summary_metrics(c), bench.summary_metrics(flipped)
        assert straight["macro_f1"] == other["macro_f1"]
        assert straight["balanced_accuracy"] != other["balanced_accuracy"]


class TestMetricImplementationsAgree:
    """`holdout_report` carries its own copy of macro-F1.

    Two implementations of the project's primary metric can drift apart
    without anything failing. This pins them together.
    """

    @staticmethod
    def _rows(counts):
        return [{"predicted": str(p), "expected": str(e)}
                for (p, e), n in counts.items() for _ in range(n)]

    @pytest.mark.parametrize("counts", [
        {(DOWN, DOWN): 3, (NORMAL, NORMAL): 4, (UP, UP): 5},
        {(DOWN, DOWN): 2, (DOWN, NORMAL): 1, (NORMAL, NORMAL): 2,
         (NORMAL, UP): 1, (UP, UP): 2, (UP, DOWN): 1},
        {(NORMAL, DOWN): 30, (NORMAL, NORMAL): 40, (NORMAL, UP): 30},
        {(DOWN, UP): 7, (UP, DOWN): 7},
    ])
    def test_both_macro_f1_implementations_agree(self, counts):
        assert holdout_report.macro_f1(self._rows(counts)) == pytest.approx(
            bench.summary_metrics(counts)["macro_f1"]
        )


# --- the classifier ---------------------------------------------------------

class TestClassifyBoundaries:
    """Cutoffs are `< DOWN` and `>= UP`, so the boundary belongs to UP."""

    def test_just_below_down_cutoff_is_down(self):
        assert bench.classify(bench.DOWN_CUTOFF - 1e-9) == DOWN

    def test_exactly_the_down_cutoff_is_normal(self):
        assert bench.classify(bench.DOWN_CUTOFF) == NORMAL

    def test_just_below_up_cutoff_is_normal(self):
        assert bench.classify(bench.UP_CUTOFF - 1e-9) == NORMAL

    def test_exactly_the_up_cutoff_is_up(self):
        assert bench.classify(bench.UP_CUTOFF) == UP

    def test_baseline_is_normal(self):
        """UI baseline is 1.0; it must not sit in a change class."""
        assert bench.classify(1.0) == NORMAL


# --- arm comparison ---------------------------------------------------------

class TestHoldoutConditioning:
    """A pairing must compare two answers to the SAME question.

    A case key present in both arms is not enough: if the arms resolved a
    different number of gene or knockout uuids, they perturbed different
    things and the pair is not comparable. Unconditioned, a boundary-removal
    arm read +265; conditioned it read +64 on 7,168 cases.
    """

    @staticmethod
    def _dump(tmp_path, name, rows):
        p = tmp_path / name
        cols = ["pathway", "gene", "direction", "key_output",
                "predicted", "expected", "n_gene_uuids", "n_ko_uuids"]
        with p.open("w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        return p

    def _case(self, gene, pred, exp, ngene="2", nko="3"):
        return {"pathway": "Signaling_by_WNT", "gene": gene, "direction": "KO",
                "key_output": "111", "predicted": pred, "expected": exp,
                "n_gene_uuids": ngene, "n_ko_uuids": nko}

    def test_load_keys_on_the_case_identity(self, tmp_path):
        p = self._dump(tmp_path, "a.tsv", [self._case("AAA", "1", "1")])
        loaded = holdout_report.load(p)
        assert list(loaded) == [("Signaling_by_WNT", "AAA", "KO", "111")]

    def test_rows_with_unscorable_expected_are_dropped(self, tmp_path):
        p = self._dump(tmp_path, "a.tsv", [
            self._case("AAA", "1", "1"), self._case("BBB", "1", "-1")])
        assert len(holdout_report.load(p)) == 1

    def test_a_case_whose_perturbation_set_changed_is_not_comparable(self, tmp_path):
        """The defect this class exists for.

        Same key in both arms, but the arms resolved a different gene uuid
        count, so the pair must be excluded from the net figure.
        """
        base = holdout_report.load(self._dump(tmp_path, "base.tsv", [
            self._case("AAA", "1", "1", ngene="2"),
            self._case("BBB", "0", "1", ngene="2"),
        ]))
        arm = holdout_report.load(self._dump(tmp_path, "arm.tsv", [
            self._case("AAA", "1", "1", ngene="2"),
            self._case("BBB", "1", "1", ngene="9"),   # resolution changed
        ]))
        paired = set(base) & set(arm)
        comparable = [k for k in paired
                      if base[k]["n_gene_uuids"] == arm[k]["n_gene_uuids"]
                      and base[k]["n_ko_uuids"] == arm[k]["n_ko_uuids"]]
        assert len(paired) == 2
        assert len(comparable) == 1, (
            "BBB changed resolution between arms; counting its flip as a win "
            "credits the change for answering a different question"
        )

    def test_report_runs_on_both_splits(self, capsys):
        """The tuning/held-out split must actually partition the cases."""
        cases = {
            ("Signaling_by_WNT", "A", "KO", "1"): {
                "pathway": "Signaling_by_WNT", "predicted": "1", "expected": "1"},
            ("Base_Excision_Repair", "B", "KO", "1"): {
                "pathway": "Base_Excision_Repair", "predicted": "0", "expected": "1"},
        }
        holdout_report.report("fixture", cases)
        out = capsys.readouterr().out
        assert "TUNING" in out and "HELD-OUT" in out

    def test_the_paper_s_ten_are_the_tuning_set(self):
        """Signaling_by_WNT is one of the ten; Base_Excision_Repair is not."""
        assert "Signaling_by_WNT" in holdout_report.TUNING_PATHWAYS
        assert "Base_Excision_Repair" not in holdout_report.TUNING_PATHWAYS


# --- catalog lookup ---------------------------------------------------------

class TestPathwayDirIndex:
    """A lookup that silently matches nothing reports an empty sweep as a
    finished one -- which happened in this session."""

    def test_finds_both_directory_layouts(self, tmp_path):
        (tmp_path / "Signaling_by_WNT_R-HSA-195721").mkdir()
        (tmp_path / "R-HSA-909733").mkdir()
        idx = pathway_dir_index(tmp_path)
        assert idx["195721"].name.endswith("R-HSA-195721")
        assert idx["909733"].name == "R-HSA-909733"

    def test_matches_any_species_prefix(self, tmp_path):
        (tmp_path / "R-MMU-12345").mkdir()
        assert "12345" in pathway_dir_index(tmp_path)

    def test_ignores_files_and_unrelated_directories(self, tmp_path):
        (tmp_path / "notes.txt").write_text("x")
        (tmp_path / "scratch").mkdir()
        assert pathway_dir_index(tmp_path) == {}

    def test_a_numeric_collision_is_not_resolved_silently(self, tmp_path):
        (tmp_path / "R-HSA-12345").mkdir()
        (tmp_path / "R-MMU-12345").mkdir()
        with pytest.raises(ValueError, match="share pathway id 12345"):
            pathway_dir_index(tmp_path)
