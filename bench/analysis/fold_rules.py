"""The record's fold precision, shared by every fitter and scorer that reads replayed folds.

The benches print a predicted fold with 6 significant figures (korkut, decryptM: `:.6g`) or 6 decimals
(curator, experimental: `:.6f`), so in the recorded case files an untouched readout reads exactly 1.
A fold replayed straight from a solve carries float noise (1.0000000000000218, an entity mean over many
copies). Every label the record's tables were fitted on is therefore read off the fold rounded to 6
significant figures. No recorded fold lies within 1e-3 of 1 except exactly 1 (2026-11-20 review), so
either print precision gives the same classes. The engine applies the same rule (`record_fold` in
src/api/readouts.jl).
"""


def record_fold(fold):
    """The fold at the record's printed precision (6 significant figures)."""
    return float(f"{fold:.6g}")


def untouched(fold):
    """A readout the perturbation does not reach: its fold reads 1 at the record's precision."""
    return record_fold(fold) == 1.0
