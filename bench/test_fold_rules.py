"""2026-11-20 review: the record-precision fold rule shared by the fitters and the engine."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "review_1120"))
import fold_rules as f  # noqa: E402
import untouched_parity as u  # noqa: E402


def test_noise_reads_untouched_and_real_moves_do_not():
    noisy = sum([0.1 * 0.1] * 3) / 3 / 0.01          # computed, as an entity mean is
    assert noisy != 1.0 and f.untouched(noisy)
    assert not f.untouched(1.00001) and not f.untouched(0.99999)


def test_record_fold_matches_the_printed_precision():
    assert f.record_fold(0.0100000000011) == 0.01
    assert f.record_fold(98.99999999999) == 99.0
    assert f.record_fold(0.0) == 0.0


def test_compare_counts_both_rules():
    k = ("korkut", "P", "d", "", "r")
    replay = {k: 1.0000000000000218, ("korkut", "P", "d", "", "s"): 0.010000000011}
    rec = {k: 1.0, ("korkut", "P", "d", "", "s"): 0.01}
    c = u.compare(replay, rec)
    assert c["joined"] == 2 and c["record_untouched"] == 1
    assert c["untouched_mismatch_exact"] == 1 and c["untouched_mismatch_record"] == 0
    assert c["railed_mismatch_exact"] == 1 and c["railed_mismatch_record"] == 0
