#!/usr/bin/env python3
"""Fail if CLAUDE.md's assertion table disagrees with what the suites printed.

Why this exists. Three counts in that table were wrong from the day they were
written and survived months of edits: `test_and_curves.jl` was listed at 17
when it had 71, `test_cli_observations.jl` at 20 when it had 39, and
`test_cycle_handling.jl` at 39+2 when it had 34+2. The table even claimed the
figures were CI-verified. Nothing checked, so nothing caught it.

The counts are already printed by every run. This turns them into an assertion.

Input is one file per suite, named `test_<suite>.out`, each holding that
suite's output. The suite is identified by the FILE NAME, never by guessing from
surrounding text: an earlier version inferred it from nearby lines and got it
wrong in every concatenated layout tried.

Each `Test Summary:` header is read as named columns (Pass, Fail, Error, Broken,
Total, Time) and zipped with the row directly below it. Reading columns by name
rather than position matters: a `Pass Fail Total` row has the same shape as
`Pass Broken Total`, and positional parsing reads the Fail count as Broken --
which on `test_cycle_handling` (34 + 2) would have matched the table exactly and
passed a failing suite.

Exits non-zero on any disagreement, any suite reporting Fail or Error, any
unparseable summary, more than one top-level summary in one file, a suite in
one place and not the other, a duplicate table row, or a missing input path.

Usage:
  python scripts/check_doc_counts.py test-output/          # a directory
  python scripts/check_doc_counts.py a/test_x.out b/test_y.out
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CLAUDE_MD = REPO / "CLAUDE.md"
SUITE_LIST = REPO / "test" / "asserting_suites.txt"

# | `test/test_and_curves.jl` | 71 | |
# | `test/test_cycle_handling.jl` | 34 | + 2 `@test_broken` |
TABLE_ROW = re.compile(
    r"^\|\s*`test/(?P<file>test_\w+)\.jl`\s*\|\s*(?P<count>[^|]*?)\s*\|"
    r"(?P<note>[^|]*)\|\s*$"
)
BROKEN_NOTE = re.compile(r"\+\s*(\d+)\s*`@test_broken`")
# Julia prints e.g. "17.7s" or, past a minute, "1m01.1s".
DURATION = re.compile(r"^(?:\d+m)?\d+(?:\.\d+)?s$")


def read_expected(problems: list[str]) -> dict[str, tuple[int, int]]:
    """file stem -> (assertions, broken) as CLAUDE.md claims."""
    out: dict[str, tuple[int, int]] = {}
    for n, line in enumerate(CLAUDE_MD.read_text().splitlines(), 1):
        m = TABLE_ROW.match(line)
        if not m:
            continue
        name, count = m.group("file"), m.group("count")
        if not count.isdigit():
            problems.append(f"CLAUDE.md:{n}: {name} has a non-numeric count {count!r}")
            continue
        if name in out:
            problems.append(f"CLAUDE.md:{n}: {name} has more than one row in the table")
            continue
        broken = BROKEN_NOTE.search(m.group("note"))
        out[name] = (int(count), int(broken.group(1)) if broken else 0)
    return out


def parse_summary(header: str, row: str) -> dict[str, int] | str:
    """Zip a Test Summary header's column names with the row's numbers.

    Returns the column -> value mapping, or a string describing why it could
    not be parsed.
    """
    if "|" not in header or "|" not in row:
        return "summary header or row has no '|' separator"
    cols = header.split("|", 1)[1].split()
    vals = row.split("|", 1)[1].split()
    if not cols or cols[-1] != "Time":
        return f"unexpected summary header columns {cols}"
    if not vals or not DURATION.match(vals[-1]):
        return f"summary row does not end in a duration: {row.strip()!r}"
    names, nums = cols[:-1], vals[:-1]
    if len(names) != len(nums) or not all(v.isdigit() for v in nums):
        return f"summary row does not match its header: {cols} vs {vals}"
    return {k: int(v) for k, v in zip(names, nums)}


def read_actual(paths: list[Path], problems: list[str]) -> dict[str, tuple[int, int]]:
    """file stem -> (pass, broken), read from each suite's own output file."""
    out: dict[str, tuple[int, int]] = {}
    for p in paths:
        name = p.stem
        if not name.startswith("test_"):
            problems.append(
                f"{p.name}: input files must be named test_<suite>.out so the suite "
                "is identified by name, not guessed from its contents"
            )
            continue
        lines = p.read_text(errors="replace").splitlines()
        found = 0
        for i, line in enumerate(lines):
            if not line.lstrip().startswith("Test Summary:"):
                continue
            found += 1
            if found > 1:
                problems.append(
                    f"{p.name}: more than one top-level Test Summary; nest the "
                    "file's testsets in one outer testset so every assertion is counted"
                )
                break
            row = lines[i + 1] if i + 1 < len(lines) else ""
            parsed = parse_summary(line, row)
            if isinstance(parsed, str):
                problems.append(f"{p.name}: {parsed}")
                continue
            bad = {k: v for k, v in parsed.items() if k in ("Fail", "Error") and v}
            if bad:
                problems.append(
                    f"{p.name}: the suite did not pass ({', '.join(f'{k} {v}' for k, v in bad.items())}); "
                    "a Pass count understates a file that did not fully pass"
                )
                continue
            out[name] = (parsed.get("Pass", 0), parsed.get("Broken", 0))
        if found == 0:
            problems.append(f"{p.name}: no Test Summary found")
    return out


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    problems: list[str] = []
    paths: list[Path] = []
    for a in sys.argv[1:]:
        p = Path(a)
        if not p.exists():
            print(f"check_doc_counts: input path does not exist: {a}", file=sys.stderr)
            return 2
        paths.extend(sorted(p.glob("*.out")) if p.is_dir() else [p])
    if not paths:
        print("check_doc_counts: no .out files in the given input", file=sys.stderr)
        return 2

    expected = read_expected(problems)
    actual = read_actual(paths, problems)
    listed = {
        ln.strip() for ln in SUITE_LIST.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    }
    seen = {p.stem for p in paths}

    for name in sorted(listed - seen):
        problems.append(f"{name}: listed in asserting_suites.txt but has no output file")
    for name in sorted(seen - listed):
        if name.startswith("test_"):
            problems.append(f"{name}: has an output file but is not in asserting_suites.txt")
    for name in sorted(listed - set(expected)):
        problems.append(f"{name}: runs in CI but has no row in CLAUDE.md's table")
    for name in sorted(set(expected) - listed):
        problems.append(f"{name}: has a row in CLAUDE.md's table but is not run")

    for name in sorted(set(expected) & set(actual)):
        want, want_broken = expected[name]
        got, got_broken = actual[name]
        if (want, want_broken) != (got, got_broken):
            problems.append(
                f"{name}: CLAUDE.md says {want}"
                f"{f' + {want_broken} broken' if want_broken else ''}, "
                f"the suite printed {got}"
                f"{f' + {got_broken} broken' if got_broken else ''}"
            )

    if problems:
        print("CLAUDE.md's assertion table does not match what the suites printed:\n")
        for pr in problems:
            print(f"  - {pr}")
        print(
            "\nIf a count differs, update the table in CLAUDE.md to the printed "
            "count. Do not adjust the counts to match the table."
        )
        return 1

    print(f"CLAUDE.md assertion table matches all {len(actual)} suites.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
