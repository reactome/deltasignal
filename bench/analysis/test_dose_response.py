"""Tests for dose_response.py. M1 is a pass/fail gate, so its definitions are
pinned: what counts as moving, as monotone, as railed, and the slope."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dose_response import (  # noqa: E402
    away, identity, input_values, load_ladder, monotone, moves, railed,
    summarise, transfer_slope, fname, LADDER,
)


def test_moves_needs_a_real_departure_from_baseline():
    assert not moves([1.0, 1.0, 1.0, 1.0])
    assert not moves([1.0, 1.0000001, 1.0, 1.0])            # float noise
    assert moves([1.0, 1.0, 1.0, 2.0])


def test_monotone_accepts_either_direction_and_flat_stretches():
    assert monotone([1.2, 2.0, 5.0, 40.0])                   # activated readout
    assert monotone([0.9, 0.5, 0.1, 0.0])                    # inhibited readout
    assert monotone([2.0, 2.0, 5.0, 5.0])                    # plateaus are fine
    assert monotone([100.0, 100.0, 100.0, 100.0])            # railed throughout


def test_monotone_rejects_a_reversal():
    assert not monotone([1.5, 3.0, 2.0, 4.0])                # up, down, up
    assert not monotone([2.0, 5.0, 1.0, 0.5])                # up then crosses below


def test_railed_at_zero_and_at_the_cap_only():
    assert railed(0.0) and railed(100.0) and railed(99.995)
    assert not railed(50.0) and not railed(1.0) and not railed(0.01)


def test_transfer_slope_recovers_a_power_law():
    ins = [2.0, 5.0, 20.0, 80.0]
    assert abs(transfer_slope(ins, ins) - 1.0) < 1e-9        # identity: slope 1
    halved = [i ** 0.5 for i in ins]                          # sqrt damping
    assert abs(transfer_slope(ins, halved) - 0.5) < 1e-9


def test_slope_is_direction_agnostic():
    # A knockdown ladder that passes through unchanged is also slope 1.
    ins = input_values("0")[:3]                               # 0.5, 0.2, 0.05
    assert abs(transfer_slope(ins, ins) - 1.0) < 1e-9


def test_summary_counts_a_reversal_and_a_switch():
    cases = {
        ("P", "G1", "2", "k1"): [1.5, 2.0, 5.0, 20.0],        # graded, monotone
        ("P", "G2", "2", "k2"): [100.0, 100.0, 100.0, 100.0],  # railed at mildest
        ("P", "G3", "2", "k3"): [1.5, 3.0, 2.0, 4.0],          # reversal
        ("P", "G4", "2", "k4"): [1.0, 1.0, 1.0, 1.0],          # never moves
    }
    s = summarise(cases)
    assert s["moving"] == 3
    assert s["monotone"] == 2
    assert s["railed_at_mildest"] == 1
    assert len(s["reversals"]) == 1


# --- Pins added after review of PR #69 ---

def test_identity_readouts_are_recognised_and_excluded_from_m1_denominator():
    ins = input_values("2")                                   # 2, 5, 20, 80
    assert identity(ins, list(ins))
    assert not identity(ins, [2.0, 5.0, 20.0, 60.0])
    s = summarise({("P", "G", "2", "self"): list(ins),
                   ("P", "G", "2", "k"): [1.5, 3.0, 2.0, 4.0]})
    assert s["identity"] == 1 and s["nontrivial"] == 1
    assert s["nt_monotone"] == 0                              # identity cannot prop it up


def test_away_requires_growing_distance_not_just_one_direction():
    assert away([0.9, 0.5, 0.1, 0.0])
    assert away([1.5, 3.0, 3.0, 9.0])
    assert monotone([0.5, 0.8, 1.5, 3.0]) and not away([0.5, 0.8, 1.5, 3.0])   # crosses baseline
    assert monotone([3.0, 2.0, 1.5, 1.2]) and not away([3.0, 2.0, 1.5, 1.2])   # falls back to baseline


def test_print_unit_tolerance_absorbs_a_single_rounding_step():
    outs = [0.5, 0.400001, 0.4, 0.400001]                     # one 1e-6 print unit back
    assert not monotone(outs)
    assert monotone(outs, 1e-6)


def test_knockdown_slope_ignores_the_zero_step():
    # Output passes the finite steps through exactly, then reads 1e-3 at input 0.
    # With 0 floored to 1e-6 that one point would drag the slope far below 1.
    ins = input_values("0")
    outs = ins[:3] + [1e-3]
    assert abs(transfer_slope(ins, outs) - 1.0) < 1e-9


def test_load_ladder_joins_on_the_case_key_and_warns_on_drops(tmp_path, capsys):
    cols = "pathway\tgene\tdirection\tkey_output\tvalid\tpred_ui\n"
    for n, (d, u) in enumerate(LADDER):
        rows = cols + f"P\tG\t2\tk\t1\t{u}\n"
        if n == 0:
            rows += "P\tG\t2\tonly_here\t1\t3\n"
        (tmp_path / fname(d, u)).write_text(rows)
    cases = load_ladder(str(tmp_path))
    assert list(cases) == [("P", "G", "2", "k")]
    assert cases[("P", "G", "2", "k")] == [u for _, u in LADDER]
    assert "1 cases are not valid in every run" in capsys.readouterr().err


def test_identity_tolerance_is_relative_so_both_ends_of_the_ladder_match():
    kd = input_values("0")                                    # 0.5 .. 0.05, 0
    assert identity(kd, [0.5, 0.2, 0.05, 0.0])
    assert not identity(kd, [0.5, 0.2, 0.0501, 0.0])          # 0.2% off at 0.05: not a copy
    oe = input_values("2")
    assert identity(oe, [2.0, 5.0, 20.0, 80.004])            # 5e-5 relative: a copy


def test_away_print_tolerance_reaches_the_distance_check():
    outs = [0.5, 0.4, 0.400001, 0.3]                          # one rounding step back toward 1
    assert not away(outs)
    assert away(outs, 1e-6, 1e-5)


def test_m2_counts_non_identity_readouts_only():
    ins = input_values("2")
    s = summarise({("P", "G", "2", "self"): list(ins),
                   ("P", "G", "2", "k"): [100.0, 100.0, 100.0, 100.0]})
    assert s["railed_at_mildest"] == 1 and s["nontrivial"] == 1
