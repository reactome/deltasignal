#!/usr/bin/env python3
"""specs/160 FR-001: derive the per-gene effective half-life table.

Source: Mathieson et al. 2018, Nat Commun 9:689 (CC BY 4.0), Supplementary
Data 2 (41467_2018_3106_MOESM5_ESM.xlsx), sheet "protein half lives high
qual". Per gene: the median of the finite half-lives across the eight HUMAN
replicate columns (B cells, NK cells, hepatocytes, monocytes; mouse neurons
excluded). t_eff = t_gene x (46 h / median over all genes of t_gene) — the
dataset supplies RELATIVE turnover, anchored at the dividing-cell median
(Schwanhäusser 2011), as registered in spec.md. Nothing fitted.

The xlsx is read with the standard library (zipped XML), so no third-party
code ever touches the downloaded file.

Usage: derive_halflives.py MOESM5_XLSX OUT_TSV
"""
import math
import re
import statistics
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
SHEET = "protein half lives high qual"
HUMAN = ("Bcells", "NK cells", "Hepatocytes", "Monocytes")
ANCHOR_H = 46.0


def col(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_sheet(path, wanted):
    z = zipfile.ZipFile(path)
    shared = [
        "".join(t.text or "" for t in si.iter("{%s}t" % NS["m"]))
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target")
            for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    for s in wb.find("m:sheets", NS):
        if s.get("name") != wanted:
            continue
        t = rels[s.get(RID)].lstrip("/")
        t = t if t.startswith("xl/") else "xl/" + t
        rows = []
        for row in ET.fromstring(z.read(t)).iter("{%s}row" % NS["m"]):
            vals = {}
            for c in row.findall("m:c", NS):
                v = c.find("m:v", NS)
                if v is None:
                    continue
                vals[col(c.get("r"))] = shared[int(v.text)] if c.get("t") == "s" else v.text
            if vals:
                rows.append([vals.get(i, "") for i in range(max(vals) + 1)])
        return rows
    raise SystemExit(f"sheet {wanted!r} not found")


def main():
    rows = read_sheet(sys.argv[1], SHEET)
    header = rows[0]
    cols = [i for i, h in enumerate(header)
            if h.endswith("half_life") and h.startswith(HUMAN)]
    if len(cols) != 8:
        raise SystemExit(f"expected 8 human half-life columns, found {len(cols)}")
    per_gene = {}
    for r in rows[1:]:
        gene = r[0].strip() if r else ""
        if not gene:
            continue
        vals = []
        for i in cols:
            try:
                v = float(r[i]) if i < len(r) else float("nan")
            except ValueError:
                continue
            if math.isfinite(v) and v > 0:
                vals.append(v)
        if vals:
            per_gene[gene] = (statistics.median(vals), len(vals))
    med = statistics.median(t for t, _ in per_gene.values())
    scale = ANCHOR_H / med
    with open(sys.argv[2], "w") as f:
        f.write("# Mathieson et al. 2018, Nat Commun 9:689 (CC BY 4.0), Supplementary Data 2,\n")
        f.write(f"# sheet '{SHEET}'. t_human_median: median of finite human replicate half-lives (h).\n")
        f.write(f"# t_eff = t_human_median x {ANCHOR_H} / {med:.4f} (the all-gene median), specs/160.\n")
        f.write("gene\tt_human_median_h\tn_replicates\tt_eff_h\n")
        for g in sorted(per_gene):
            t, n = per_gene[g]
            f.write(f"{g}\t{t:.4f}\t{n}\t{t * scale:.4f}\n")
    print(f"{len(per_gene)} genes; all-gene human median {med:.2f} h -> scale {scale:.4f}")


if __name__ == "__main__":
    main()
