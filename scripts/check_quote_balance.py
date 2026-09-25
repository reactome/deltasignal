#!/usr/bin/env python3
"""Fail if any markdown paragraph has unbalanced double quotes.

Written because a bulk paraphrasing pass left five dangling closing quotes: it
rewrote an attribution and the FIRST line of a multi-line quote, leaving the
rest of the quote body and its closing punctuation stranded mid-paragraph. All
five read as normal prose to a reviewer and none was caught by eye.

What is checked is prose only. Before paragraphs are split, whole lines are
removed when they sit inside a fenced block (``` or ~~~, at any indentation),
inside an indented code block, or are table rows. Splitting first and skipping
"paragraphs that start with a fence" was wrong: a fenced block containing a
blank line then had its second half checked as prose.

Within prose, these legitimate lone quotes are ignored: an inches or
arc-seconds mark after a digit (12" rack, 39' 12"), and an escaped \\".
Code spans are stripped only when a paragraph's backticks are balanced; with an
odd count, stripping would swallow text, so quotes are counted unstripped.

Runs over every tracked .md from the repository root, whatever the current
directory. Exits 2 if it cannot list files or finds none, and reports any file
it cannot read rather than skipping it.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

FENCE = re.compile(r"^\s*(```|~~~)")
TABLE = re.compile(r"^\s*\|")
INDENTED = re.compile(r"^( {4,}|\t)")
LIST_ITEM = re.compile(r"^\s*([-*+]|\d+[.)])\s")
ESCAPED_QUOTE = re.compile(r'\\"')
# A quote right after a digit MIGHT be an inches or arc-seconds mark (12" rack),
# or might be an ordinary closing quote ("+12", width="100"). So it is only
# ever used to RESCUE a paragraph that is otherwise odd, never to break one that
# is already even. An earlier version stripped these unconditionally and turned
# every quoted number into a false positive -- 15 in this repo alone.
INCH_MARK = re.compile(r'(?<=\d)"')
CODE_SPAN = re.compile(r"`[^`]*`")


def repo_files() -> tuple[Path, list[Path]]:
    try:
        root = Path(subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True).stdout.strip())
        raw = subprocess.run(
            ["git", "ls-files", "-z", "*.md"],
            capture_output=True, check=True, cwd=root).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"check_quote_balance: cannot list tracked files: {exc}", file=sys.stderr)
        sys.exit(2)
    names = [n for n in raw.decode("utf-8", "surrogateescape").split("\0") if n]
    return root, [root / n for n in names]


def prose_lines(lines: list[str]) -> list[str | None]:
    """Same length as `lines`; None where a line is code, a fence or a table row."""
    out: list[str | None] = []
    in_fence = False
    in_indented = False
    prev_blank = True
    prev_listish = False
    for line in lines:
        if FENCE.match(line):
            in_fence = not in_fence
            out.append(None)
            prev_blank = False
            continue
        if in_fence:
            out.append(None)
            continue
        blank = not line.strip()
        if INDENTED.match(line) and (in_indented or (prev_blank and not prev_listish)):
            in_indented = True
            out.append(None)
            prev_blank = False
            continue
        if not blank:
            in_indented = False
        if TABLE.match(line):
            out.append(None)
        else:
            out.append(line)
        if not blank:
            prev_listish = bool(LIST_ITEM.match(line)) or (prev_listish and line.startswith(" "))
        prev_blank = blank
    return out


def unbalanced(para: str) -> bool:
    body = para
    if body.count("`") % 2 == 0:
        body = CODE_SPAN.sub("", body)
    body = ESCAPED_QUOTE.sub("", body.replace('"""', ""))
    curly = body.count("“") != body.count("”")
    odd = body.count('"') % 2 == 1
    odd_without_inches = INCH_MARK.sub("", body).count('"') % 2 == 1
    return curly or (odd and odd_without_inches)


def main() -> int:
    root, files = repo_files()
    if not files:
        print("check_quote_balance: no tracked .md files found", file=sys.stderr)
        return 2
    hits: list[str] = []
    unreadable: list[str] = []
    for path in files:
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            unreadable.append(f"{path.relative_to(root)}: {exc}")
            continue
        para: list[str] = []
        start = 0
        for n, line in enumerate(prose_lines(lines) + [""], 1):
            if line is None or not line.strip():
                if para and unbalanced("\n".join(para)):
                    hits.append(f"{path.relative_to(root)}:{start}: {para[0].strip()[:88]}")
                para = []
                continue
            if not para:
                start = n
            para.append(line)
    for h in hits:
        print(h)
    for u in unreadable:
        print(f"UNREADABLE {u}")
    print(f"\n{len(hits)} unbalanced paragraph(s) across {len(files)} files"
          + (f", {len(unreadable)} unreadable" if unreadable else ""))
    return 1 if hits or unreadable else 0


if __name__ == "__main__":
    sys.exit(main())
