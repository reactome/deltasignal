#!/usr/bin/env python3

import unittest

from benchmark_mpbiopath_cases import (
    DOWN,
    NORMAL,
    UP,
    metric_summary_for_field,
    reachable_path_signs,
    structural_prediction,
)


class StructuralBaselineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.graph = {
            "A": [("B", 1), ("C", 1)],
            "B": [("C", -1)],
            "C": [("A", 1)],
        }

    def test_all_reachable_signs_preserve_both_polarities(self) -> None:
        signs = reachable_path_signs(
            self.graph,
            ["A"],
            ["C"],
            shortest_only=False,
        )
        self.assertEqual(signs, {-1, 1})
        self.assertEqual(structural_prediction(signs, UP), NORMAL)

    def test_shortest_path_uses_only_minimum_distance(self) -> None:
        signs = reachable_path_signs(
            self.graph,
            ["A"],
            ["C"],
            shortest_only=True,
        )
        self.assertEqual(signs, {1})
        self.assertEqual(structural_prediction(signs, UP), UP)
        self.assertEqual(structural_prediction(signs, DOWN), DOWN)

    def test_negative_path_flips_perturbation_direction(self) -> None:
        signs = reachable_path_signs(
            self.graph,
            ["B"],
            ["C"],
            shortest_only=True,
        )
        self.assertEqual(structural_prediction(signs, UP), DOWN)
        self.assertEqual(structural_prediction(signs, DOWN), UP)

    def test_no_path_is_conservatively_unchanged(self) -> None:
        signs = reachable_path_signs(
            self.graph,
            ["missing"],
            ["C"],
            shortest_only=False,
        )
        self.assertEqual(signs, set())
        self.assertEqual(structural_prediction(signs, UP), NORMAL)


class BaselineMetricTests(unittest.TestCase):
    def test_unavailable_structural_prediction_reduces_coverage(self) -> None:
        rows = [
            {
                "expected": UP,
                "baseline": UP,
                "mapping_status": "mapped",
            },
            {
                "expected": DOWN,
                "baseline": None,
                "mapping_status": "key_output_absent",
            },
        ]
        summary = metric_summary_for_field(rows, "baseline")
        self.assertEqual(summary["eligible_cases"], 2)
        self.assertEqual(summary["scored_cases"], 1)
        self.assertEqual(summary["correct"], 1)
        self.assertEqual(summary["coverage"], 0.5)

    def test_convergence_reports_cases_and_unique_solves_separately(self) -> None:
        rows = [
            {
                "pathway_id": "R-HSA-1",
                "gene": "TP53",
                "direction": UP,
                "expected": UP,
                "prediction": UP,
                "mapping_status": "scored",
                "converged": False,
            },
            {
                "pathway_id": "R-HSA-1",
                "gene": "TP53",
                "direction": UP,
                "expected": DOWN,
                "prediction": DOWN,
                "mapping_status": "scored",
                "converged": False,
            },
            {
                "pathway_id": "R-HSA-1",
                "gene": "ATM",
                "direction": DOWN,
                "expected": DOWN,
                "prediction": DOWN,
                "mapping_status": "scored",
                "converged": True,
            },
        ]
        summary = metric_summary_for_field(rows, "prediction")
        # Baseline projections are intentionally treated as converged.
        self.assertEqual(summary["unique_nonconverged_solves"], 0)

        from benchmark_mpbiopath_cases import metric_summary

        summary = metric_summary(rows)
        self.assertEqual(summary["nonconverged_cases"], 2)
        self.assertEqual(summary["unique_scored_solves"], 2)
        self.assertEqual(summary["unique_nonconverged_solves"], 1)
        self.assertEqual(summary["converged_scored_cases"], 1)


if __name__ == "__main__":
    unittest.main()
