#!/usr/bin/env python3
"""Fail if any markdown paragraph has unbalanced double quotes.

Written because a bulk paraphrasing pass left five dangling closing quotes: it
rewrote an attribution and the FIRST line of a multi-line quote, leaving the
rest of the quote body and its closing punctuation stranded mid-paragraph. All
five read as normal prose to a reviewer and none was caught by eye.

Checks every tracked .md. Code spans and triple-quote fences are ignored, and
table rows and fenced blocks are skipped.
"""
import re, subprocess, sys

paths = subprocess.run(["git","ls-files","*.md"], capture_output=True, text=True).stdout.split()
bad = 0
for p in paths:
    try:
        text = open(p, encoding="utf-8").read()
    except Exception:
        continue
    for para in re.split(r"\n\s*\n", text):
        if para.lstrip().startswith("```") or para.lstrip().startswith("|"):
            continue
        body = re.sub(r"`[^`]*`", "", para)           # ignore code spans
        body = body.replace('"""', "")                # ignore triple-quote fences
        if body.count('"') % 2 or body.count('“') != body.count('”'):
            line = text[:text.index(para)].count("\n") + 1
            print(f"{p}:{line}: unbalanced quotes -> {para.strip().splitlines()[0][:88]}")
            bad += 1
print(f"\n{bad} unbalanced paragraph(s)")
sys.exit(1 if bad else 0)
