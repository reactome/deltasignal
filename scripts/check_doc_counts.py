#!/usr/bin/env python3
"""Fail if CLAUDE.md's assertion table disagrees with what the suites printed.

Why this exists. Three counts in that table were wrong from the day they were
written and survived months of edits: `test_and_curves.jl` was listed at 17
when it had 71, `test_cli_observations.jl` at 20 when it had 39, and
`test_cycle_handling.jl` at 39+2 when it had 34+2. The table even claimed the
figures were CI-verified. Nothing checked, so nothing caught it.

The counts are already printed by every run. This turns them into an assertion.

Reads the concatenated suite output (Julia's `Test Summary:` blocks) and the
table in CLAUDE.md, and exits non-zero on any disagreement, any file present in
one and not the other, or any suite whose summary reports failures or errors —
because a `Pass` count understates a file that did not fully pass.

Usage:
  julia ... test/foo.jl | tee out.txt          # for each suite
  python scripts/check_doc_counts.py out.txt   # or several files, or a dir
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CLAUDE_MD = REPO / "CLAUDE.md"
SUITE_LIST = REPO / "test" / "asserting_suites.txt"

# "configuration guard rails |  192    192  17.7s"  ->  name, pass, total
# "loop elasticity           |  150   1  151  4.2s" ->  name, pass, broken, total
SUMMARY_ROW = re.compile(
    r"^\s*(?P<name>\S.*?)\s*\|\s*(?P<pass>\d+)"
    r"(?:\s+(?P<broken>\d+))?\s+(?P<total>\d+)\s+[\d.]+s\s*$"
)
BAD_COLUMN = re.compile(r"\b(Fail|Error)\b")

# | `test/test_and_curves.jl` | 71 | |
# | `test/test_cycle_handling.jl` | 34 | + 2 `@test_broken` |
TABLE_ROW = re.compile(
    r"^\|\s*`test/(?P<file>test_[a-z_]+)\.jl`\s*\|\s*(?P<count>\d+)\s*\|"
    r"(?P<note>[^|]*)\|\s*$"
)
BROKEN_NOTE = re.compile(r"\+\s*(\d+)\s*`@test_broken`")


def read_expected() -> dict[str, tuple[int, int]]:
    """file stem -> (assertions, broken) as CLAUDE.md claims."""
    out: dict[str, tuple[int, int]] = {}
    for line in CLAUDE_MD.read_text().splitlines():
        m = TABLE_ROW.match(line)
        if not m:
            continue
        broken = BROKEN_NOTE.search(m.group("note"))
        out[m.group("file")] = (int(m.group("count")), int(broken.group(1)) if broken else 0)
    return out


def read_actual(paths: list[Path]) -> tuple[dict[str, tuple[int, int]], list[str]]:
    """file stem -> (pass, broken) from the suite output, plus any bad summary."""
    out: dict[str, tuple[int, int]] = {}
    problems: list[str] = []
    for p in paths:
        text = p.read_text(errors="replace")
        # The stem names the suite when the caller tees per-suite (test_foo.out).
        stem = p.stem if p.stem.startswith("test_") else None
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if not line.lstrip().startswith("Test Summary:"):
                continue
            if BAD_COLUMN.search(line):
                problems.append(f"{p.name}: summary reports failures or errors: {line.strip()}")
            for row in lines[i + 1:i + 3]:
                m = SUMMARY_ROW.match(row)
                if not m:
                    continue
                name = stem or _stem_from_context(lines, i)
                if name:
                    out[name] = (int(m.group("pass")), int(m.group("broken") or 0))
                break
    return out, problems


def _stem_from_context(lines: list[str], idx: int) -> str | None:
    """Find `test_foo` in the surrounding lines when the file is a concatenation."""
    for row in lines[max(0, idx - 40):idx + 6]:
        m = re.search(r"\b(test_[a-z_]+)\.jl\b", row) or re.search(r"\b(test_[a-z_]+) PASSED", row)
        if m:
            return m.group(1)
    return None


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    paths: list[Path] = []
    for a in sys.argv[1:]:
        p = Path(a)
        paths.extend(sorted(p.glob("*.out")) if p.is_dir() else [p])
    paths = [p for p in paths if p.exists()]
    if not paths:
        print("check_doc_counts: no readable input files", file=sys.stderr)
        return 2

    expected = read_expected()
    actual, problems = read_actual(paths)

    listed = {
        ln.strip() for ln in SUITE_LIST.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    }

    if not actual:
        problems.append("no Test Summary blocks found in the given output")

    for name in sorted(listed - set(actual)):
        problems.append(f"{name}: listed in asserting_suites.txt but produced no summary")
    for name in sorted(set(actual) - listed):
        problems.append(f"{name}: produced a summary but is not in asserting_suites.txt")
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
            "\nUpdate the table in CLAUDE.md to the printed counts. Do not adjust "
            "the counts to match the table."
        )
        return 1

    print(f"CLAUDE.md assertion table matches all {len(actual)} suites.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
