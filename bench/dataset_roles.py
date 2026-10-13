"""specs/051: the dataset register (bench/datasets.tsv) and its guard.

Every dataset has one role: FIT (fit and trace anything), SELECT (choose
between registered candidates; never fitted to directly) or TEST (scored once
per frozen release; never used to choose). A fitting or selection script calls
`require_role` before it reads a dataset, and the call refuses a role that does
not allow the purpose. Scoring is allowed on any role and returns the role so
it is written next to the number.
"""
import csv
from pathlib import Path

REGISTER = Path(__file__).resolve().parent / "datasets.tsv"
ROLES = {"FIT", "SELECT", "TEST"}
ALLOWED = {
    "fit": {"FIT"},
    "select": {"FIT", "SELECT"},
    "score": {"FIT", "SELECT", "TEST"},
}


def load(path=REGISTER):
    rows = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["role"] not in ROLES:
            raise ValueError(f"{path}: {r['id']} has role {r['role']!r}, not one of {sorted(ROLES)}")
        if r["id"] in rows:
            raise ValueError(f"{path}: {r['id']} is registered twice")
        rows[r["id"]] = r
    return rows


def require_role(dataset, purpose, path=REGISTER):
    """The dataset's role, or SystemExit if `purpose` is not allowed for it."""
    if purpose not in ALLOWED:
        raise ValueError(f"purpose {purpose!r} must be one of {sorted(ALLOWED)}")
    reg = load(path)
    if dataset not in reg:
        raise SystemExit(f"{dataset!r} is not in {path}: register it (specs/051) before using it")
    role = reg[dataset]["role"]
    if role not in ALLOWED[purpose]:
        raise SystemExit(f"{dataset!r} is registered {role}; {purpose!r} is allowed only on "
                         f"{sorted(ALLOWED[purpose])} (specs/051)")
    return role
