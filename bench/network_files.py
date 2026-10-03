"""Read a pathway's network the way the solver does (specs/044).

The generator writes the curated network to logic_network.csv and the edges it
derives at root inputs and terminal outputs (assembly, dissociation) to
boundary_edges.csv beside it, same columns. Every script that reads a network
must read both, or it silently sees a different network from the one solved.
`open_network(d)` returns one text stream: logic_network.csv, then the rows of
boundary_edges.csv without its header. It drops in wherever
`open(d / "logic_network.csv")` was used (csv.reader, csv.DictReader, or line
iteration). A bundle from before the split has no boundary file.
"""
import csv
import io
from pathlib import Path


def _header(text: str):
    return next(csv.reader(io.StringIO(text.lstrip("\ufeff").partition("\n")[0])), [])


def open_network(pathway_dir) -> io.StringIO:
    d = Path(pathway_dir)
    with open(d / "logic_network.csv", newline="") as f:
        text = f.read()
    boundary = d / "boundary_edges.csv"
    if boundary.exists():
        with open(boundary, newline="") as f:
            header, _, rows = f.read().partition("\n")
        # Compared as parsed columns, as the solver does (quoting, a BOM or
        # trailing whitespace are not a mismatch).
        own = [c.strip() for c in _header(text)]
        if [c.strip() for c in _header(header)] != own:
            raise ValueError(f"{boundary}: header {header!r} differs from logic_network.csv {own!r}")
        if rows:
            if not text.endswith("\n"):
                text += "\n"
            text += rows
    return io.StringIO(text, newline="")
