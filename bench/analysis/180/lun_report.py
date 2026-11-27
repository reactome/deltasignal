#!/usr/bin/env python3
"""specs/180: the Lun 2019 report of record (Adam, 2026-11-27).

Every arm side by side (model + tyrosine-kinase rule, model alone, reachability, rule alone, never), per half and
pooled, and beside them the same arms split by first-hop engagement (specs/179): proteins that raise a direct
substrate, proteins that do not, and proteins with no measured direct pair. Direct pairs are left out of the split
rows, as korkut's own markers are (specs/138). The raw rows are the headline; the split describes and never gates.
Both halves are FIT: nothing here is a held-out test.

Usage: lun_report.py PREDICTIONS A_FILE B_FILE [CLASS_DIR] > report.json
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
for d in ("173", "176", "177", "179"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import engage_179 as E  # noqa: E402
import lun_a_173 as A  # noqa: E402
import score_177 as S  # noqa: E402

SPLITS = ("engaged", "not engaged", "no measured direct pair")


def split_of(gene, eng):
    return SPLITS[2] if gene not in eng else SPLITS[0] if eng[gene] else SPLITS[1]


def rows(cases):
    return {"pairs": len(cases), "proteins": len({c["gene"] for c in cases}), **S.arm_scores(cases)} if cases else None


def report(m, pred, direct, tk, tyr):
    cases = S.build_cases(m, pred, tk, tyr)
    eng = E.engagement(m, direct)
    out = {"all (raw)": rows(cases)}
    rest = [c for c in cases if (c["gene"], c["readout"]) not in direct]
    for s in SPLITS:
        out[s] = rows([c for c in rest if split_of(c["gene"], eng) == s])
    return out


def main():
    import dataset_roles
    pred_file, a_file, b_file = sys.argv[1:4]
    tk, tyr = S.read_classes(sys.argv[4]) if len(sys.argv) > 4 else S.read_classes()
    genes, direct = A.marker_genes(), E.direct_table()
    res, m_all, p_all = {}, {}, {}
    for half, path in (("A", a_file), ("B", b_file)):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        m, pred = S.measured(path, genes, half), S.predictions(pred_file, half)
        res[f"half_{half}"] = report(m, pred, direct, tk, tyr)
        m_all.update(m)
        p_all.update(pred)
    res["pooled"] = report(m_all, p_all, direct, tk, tyr)
    json.dump(res, sys.stdout, indent=1, default=str)
    print()


if __name__ == "__main__":
    main()
