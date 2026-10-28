#!/usr/bin/env python3
"""specs/121: read CPPA delta-p measurements, for one portion of the registered split only.

`load(unpacked, split_tsv, portion)` returns single-compound samples of that portion: "train", or one of the test
portions "A" (unseen drug), "B" (unseen cell line), "C" (both). A test portion is read only when the caller passes
`final=True`: the registered single scoring of every arm at the end (specs/121). Development code passes "train".
Each sample: UID, set, drug, cell line, time, dose, and {antibody: delta-p}. The values never leave this machine
(the CPPA terms); nothing here writes them to the repository.
"""
import csv
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import split_121  # noqa: E402

PORTIONS = ("train", "A", "B", "C")


def portion_of(drug_role, cell_role):
    if drug_role in ("train", "seen", "stimulus"):
        return "train" if cell_role == "train" else "B"
    return "A" if cell_role == "train" else "C"


def roles(split_tsv):
    out = {"drug": {}, "cell_line": {}}
    for r in csv.DictReader(open(split_tsv), delimiter="\t"):
        out[r["unit"]][r["name"]] = r["role"]
    return out


def design(unpacked):
    """{(set, UID): (drug, cell line, time, dose)} for single-compound samples, from the treatment summaries."""
    out = {}
    for f in sorted(glob.glob(f"{unpacked}/delta-p_TreatmentSummary/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t"):
            comps = [((r.get(f"compound_name_{i}") or "").strip(), (r.get(f"dosage_{i}") or "").strip())
                     for i in range(1, 5)]
            comps = [(c, d) for c, d in comps if c and c.upper() != "NA"]
            if len(comps) == 1:
                out[(r["set"].strip(), r["UID"].strip())] = (split_121.norm(comps[0][0]),
                                                            split_121.norm(r["cell_line_name"]),
                                                            r["time"].strip(), comps[0][1])
    return out


def load(unpacked, split_tsv, portion, final=False):
    if portion not in PORTIONS:
        raise ValueError(f"portion must be one of {PORTIONS}")
    if portion != "train" and not final:
        raise SystemExit(f"test portion {portion} is read only at the registered final scoring (specs/121)")
    rl, des = roles(split_tsv), design(unpacked)
    keep = {k: v for k, v in des.items()
            if portion_of(rl["drug"][v[0]], rl["cell_line"][v[1]]) == portion}
    out = []
    for f in sorted(glob.glob(f"{unpacked}/delta-p_RPPA/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        set_id = "set" + Path(f).stem.rsplit("_set", 1)[1]          # the summaries write "set101"
        with open(f, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                cells = line.rstrip("\n").split("\t")
                d = keep.get((set_id, cells[0].strip()))
                if d is None:               # another portion's sample: its values are skipped, not parsed
                    continue
                vals = {}
                for a, v in zip(header[1:], cells[1:]):
                    try:
                        vals[a] = float(v)
                    except ValueError:
                        pass
                out.append({"set": set_id, "uid": cells[0].strip(), "drug": d[0], "cell_line": d[1],
                            "time": d[2], "dose": d[3], "values": vals})
    return out
