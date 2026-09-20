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


class TestConcentrationColumns:
    """The compare report must say how CONCENTRATED a gain is, not just how big.

    Two real cases forced this. LNG #89 reported "+14 held-out, p = 0.0013" that
    was 18 cases on ONE readout. The loop-elasticity arm then reported "+19
    held-out over 18 distinct readouts" -- and 13 of the 19 were ONE gene (DOK1)
    read at 13 places. Readout count said distributed; the cause was shared.
    """

    COLS = ["pathway", "gene", "direction", "key_output",
            "predicted", "expected", "n_gene_uuids", "n_ko_uuids"]

    def _dump(self, tmp_path, name, rows):
        p = tmp_path / name
        with p.open("w", newline="") as fh:
            fh.write("\t".join(self.COLS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in self.COLS) + "\n")
        return p

    def _row(self, pathway, gene, readout, pred, exp):
        return {"pathway": pathway, "gene": gene, "direction": "2", "key_output": readout,
                "predicted": pred, "expected": exp, "n_gene_uuids": "1", "n_ko_uuids": "1"}

    def _run(self, tmp_path, base_rows, arm_rows, capsys):
        b = self._dump(tmp_path, "base.tsv", base_rows)
        a = self._dump(tmp_path, "arm.tsv", arm_rows)
        holdout_report.main(["--cases", str(b), "--compare", str(a)])
        return capsys.readouterr().out

    def test_one_gene_many_readouts_is_flagged_as_one_perturbation(self, tmp_path, capsys):
        pw = "Base_Excision_Repair"                       # held-out
        # 10 discordant cases: 8 from gene DOK1 at 8 readouts, 2 singletons.
        base = [self._row(pw, "DOK1", f"r{i}", "1", "2") for i in range(8)]
        base += [self._row(pw, "G1", "s1", "1", "2"), self._row(pw, "G2", "s2", "1", "2")]
        arm = [dict(r, predicted="2") for r in base]      # arm fixes all 10
        out = self._run(tmp_path, base, arm, capsys)
        assert "genes" in out                              # the column exists
        held = [l for l in out.splitlines() if l.startswith("HELD-OUT")][-1]
        cols = held.split()
        assert cols[-3] == "1/1"                           # pathways moved / scored
        assert cols[-2] == "10"                            # 10 distinct readouts...
        assert cols[-1] == "3"                             # ...but only 3 genes
        assert "one perturbation (DOK1" in out             # and the dominant one is named
        assert "80%" in out

    def test_a_genuinely_distributed_gain_is_not_flagged(self, tmp_path, capsys):
        # Two held-out pathways, ten genes, ten readouts: nothing is concentrated.
        # (A one-pathway fixture would rightly trip the single-pathway warning.)
        pws = ["Base_Excision_Repair", "Signaling_by_EGFR"]
        base = [self._row(pws[i % 2], f"G{i}", f"r{i}", "1", "2") for i in range(10)]
        arm = [dict(r, predicted="2") for r in base]
        out = self._run(tmp_path, base, arm, capsys)
        held = [l for l in out.splitlines() if l.startswith("HELD-OUT")][-1]
        assert held.split()[-1] == "10"                    # 10 genes
        assert held.split()[-3] == "2/2"                   # both pathways moved
        assert "one perturbation" not in out
        assert "NOT independent" not in out

    def test_single_readout_cluster_invalidates_mcnemar(self, tmp_path, capsys):
        """The LNG #89 shape: many genes, ONE readout."""
        pw = "Base_Excision_Repair"
        base = [self._row(pw, f"G{i}", "981545", "1", "2") for i in range(9)]
        arm = [dict(r, predicted="2") for r in base]
        out = self._run(tmp_path, base, arm, capsys)
        assert "McNemar does not apply" in out
        assert "single-readout" in out

    def test_coverage_fix_with_no_comparable_pairs_is_named(self, tmp_path, capsys):
        """Arms resolved different perturbation sets => nothing is like-for-like."""
        pw = "Base_Excision_Repair"
        # One case whose readout became resolvable (dropped as not like-for-like)
        # plus one unchanged comparable case, so the split is populated but has
        # NO discordant pair -- the shape LNG #89 produced.
        base = [self._row(pw, "G1", "r1", "1", "2"), self._row(pw, "G2", "r2", "1", "1")]
        arm = [dict(base[0], predicted="2", n_ko_uuids="5"), dict(base[1])]
        out = self._run(tmp_path, base, arm, capsys)
        assert "NO comparable case moved" in out
        assert "COVERAGE fix" in out


class TestCuratorOracleGraph:
    """The oracle's graph functions, on fixtures shaped like the real defects.

    specs/016's numbers (43% / 54% of curator routes severed at uuid level;
    392 -> 43 with a two-hop composition bridge) come from these functions, so
    they are pinned here without Neo4j.
    """

    def setup_method(self):
        import curator_oracle as co  # noqa: WPS433
        self.co = co

    # (reaction, class, input, output, catalyst, regulator, regulation_class)
    ROWS = [
        ("rx1", "Reaction", "U", "A", "", "", ""),
        ("rx1", "Reaction", "", "", "CAT", "", ""),
        ("rx2", "Reaction", "A", "B", "", "", ""),
        ("rx2", "Reaction", "", "", "", "INH", "NegativeRegulation"),
        ("rx3", "Reaction", "B", "T", "", "ACT", "PositiveRegulation"),
    ]

    def test_roles_and_signs_are_read_from_the_rows(self):
        adj, ronly, ents, rxns, role = self.co.build_reaction_graph(self.ROWS)
        assert rxns == {"rx1", "rx2", "rx3"}
        assert role["U"] == {"input"} and role["A"] == {"output", "input"}
        assert role["CAT"] == {"catalyst"} and role["INH"] == {"neg_regulator"} and role["ACT"] == {"pos_regulator"}
        assert ("rx2", -1) in adj["INH"] and ("rx3", +1) in adj["ACT"]
        # the reaction-only copy is independent of later composition additions
        self.co.add_composition(adj, [("CPLX", "T")], [], strict=True)
        assert "T" in adj and ("CPLX", +1) in adj["T"]
        assert "T" not in ronly

    def test_signed_reach_multiplies_signs_along_the_path(self):
        adj, *_ = self.co.build_reaction_graph(self.ROWS)
        r = self.co.signed_reach(adj, "INH")
        assert r["rx2"] == -1 and r["B"] == -1 and r["T"] == -1     # one inhibition flips everything downstream
        r = self.co.signed_reach(adj, "U")
        assert r["T"] == +1

    def test_roots_exclude_cofactors_and_non_protein_ids(self):
        rows = [("rx", "Reaction", "R-HSA-1", "R-HSA-2", "", "", ""),
                ("rx", "Reaction", "R-ALL-atp", "", "", "", ""),
                ("rx", "Reaction", "R-HSA-cof", "", "", "", "")]
        _adj, ronly, ents, _r, _role = self.co.build_reaction_graph(rows)
        roots, terms = self.co.roots_and_terminals(ronly, ents, cofactors={"R-HSA-cof"})
        assert roots == ["R-HSA-1"]          # ATP (R-ALL) and the cofactor are not perturbable roots
        assert terms == ["R-HSA-2"]

    def test_strict_composition_is_assembly_direction_only(self):
        adj = {}
        self.co.add_composition(adj, [("BIG", "small")], [("SET", "m1")], strict=True)
        assert ("BIG", +1) in adj["small"]                 # component -> container
        assert "BIG" not in adj                            # NOT container -> component (the broadcast route)
        assert ("SET", +1) in adj["m1"] and ("m1", +1) in adj["SET"]   # a set is its members, both ways
        adj2 = {}
        self.co.add_composition(adj2, [("BIG", "small")], [], strict=False)
        assert ("small", +1) in adj2["BIG"]

    def test_resolve_ids_falls_back_to_expanded_leaves(self):
        present = {"m1", "m2", "X"}
        leaves = {"SET": {"m1", "m2", "m3"}}
        assert self.co.resolve_ids("X", present, leaves) == {"X"}
        assert self.co.resolve_ids("SET", present, leaves) == {"m1", "m2"}   # m3 is not a node
        assert self.co.resolve_ids("GONE", present, leaves) == set()

    def test_first_break_locates_the_cdk4_shaped_split(self):
        # stable-id path U -> A -> B; A has two shards: the signal lands on a1 (a
        # dissociation sink), the onward edge leaves from a2.
        path = ["U", "A", "B"]
        edge_uuids = {("U", "A"): [("u", "a1")], ("A", "B"): [("a2", "b")]}
        step, frontier = self.co.first_break(path, {"u"}, edge_uuids)
        assert step == 2 and frontier == {"a1"}
        # and no break when the shards line up
        edge_uuids2 = {("U", "A"): [("u", "a2")], ("A", "B"): [("a2", "b")]}
        assert self.co.first_break(path, {"u"}, edge_uuids2) is None

    def test_simulated_composition_never_uses_a_sink_as_source(self):
        uadj = {"x_live": {"rx"}}
        s2u = {"X": ["x_live", "x_sink"], "Y": ["y1", "y2"]}
        out, added, fan = self.co.simulate_composition(uadj, [("X", "Y")], s2u, sink_uuids={"x_sink"})
        assert added == 2 and fan == [2]
        assert out["x_live"] == {"rx", "y1", "y2"}
        assert "x_sink" not in out
        # the original adjacency is untouched
        assert uadj == {"x_live": {"rx"}}

    def test_shortest_path_and_reachable_agree(self):
        adj = {"a": {"b"}, "b": {"c"}, "c": {"d"}, "z": set()}
        assert self.co.reachable(adj, {"a"}, {"d"}) and not self.co.reachable(adj, {"z"}, {"d"})
        assert self.co.shortest_path(adj, {"a"}, {"d"}) == ["a", "b", "c", "d"]
        assert self.co.shortest_path(adj, {"z"}, {"d"}) is None

    def test_drop_carriers_removes_small_molecule_links_both_ways(self):
        # A -> rx -> GTP -> rx2 -> B : the only link is through GTP, so B must NOT
        # be reachable from A once small molecules are excluded as carriers.
        adj = {"R-HSA-A": {("R-HSA-rx", 1)}, "R-HSA-rx": {("R-ALL-gtp", 1)},
               "R-ALL-gtp": {("R-HSA-rx2", 1)}, "R-HSA-rx2": {("R-HSA-B", 1)}}
        assert self.co.signed_reach(adj, "R-HSA-A").get("R-HSA-B") == 1
        pruned = self.co.drop_carriers(adj, lambda n: n.startswith("R-HSA-"))
        assert "R-HSA-B" not in self.co.signed_reach(pruned, "R-HSA-A")
        assert "R-ALL-gtp" not in pruned and all("R-ALL-gtp" not in {t for t, _ in ts} for ts in pruned.values())
        # a declared cofactor is treated the same way
        pruned2 = self.co.drop_carriers(adj, lambda n: n.startswith("R-HSA-") and n != "R-HSA-rx2")
        assert "R-HSA-B" not in self.co.signed_reach(pruned2, "R-HSA-A")

