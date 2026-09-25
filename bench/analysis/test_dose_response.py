"""Tests for dose_response.py. M1 is a pass/fail gate, so its definitions are
pinned: what counts as moving, as monotone, as railed, and the slope."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dose_response import (  # noqa: E402
    input_values, monotone, moves, railed, summarise, transfer_slope,
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
