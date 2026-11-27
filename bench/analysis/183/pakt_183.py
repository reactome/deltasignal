#!/usr/bin/env python3
"""specs/183: who lowers p-AKT in Lun 2019, and what moves with it. Both halves are FIT.

Usage: pakt_183.py PREDICTIONS OUT_JSON A_FILE B_FILE [CLASS_DIR]
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
for d in ("173", "176", "177"):
    sys.path.insert(0, str(HERE.parents[1] / d))
import budget_176 as B  # noqa: E402
import lun_a_173 as A  # noqa: E402
import score_177 as S  # noqa: E402

TARGET = "p-AKT"
BANDS = ((0, 0), (1, 3), (4, 9), (10, 99))


def band(n):
    return next(f"{lo}-{hi}" if hi < 99 else f"{lo}+" for lo, hi in BANDS if lo <= n <= hi)


def by_protein(calls):
    """{gene: {readout: call}} from {(gene, readout): call}."""
    out = defaultdict(dict)
    for (g, r), c in calls.items():
        out[g][r] = c
    return out


def co_movement(per, target=TARGET):
    """For proteins that lower `target` against the rest: UP and DOWN rates at every other readout."""
    low = [g for g, d in per.items() if d.get(target) == "down"]
    rest = [g for g, d in per.items() if target in d and d[target] != "down"]
    out = {}
    for r in sorted({r for d in per.values() for r in d} - {target}):
        row = {}
        for lab, gs in (("lowers", low), ("rest", rest)):
            cs = [per[g][r] for g in gs if r in per[g]]
            row[lab] = {"n": len(cs), "up": sum(c == "up" for c in cs) / len(cs) if cs else None,
                        "down": sum(c == "down" for c in cs) / len(cs) if cs else None}
        out[r] = row
    return {"proteins_lowering": len(low), "proteins_rest": len(rest), "readouts": out}


def activity(per, target=TARGET):
    """Share of proteins lowering `target`, by how many OTHER readouts each changes."""
    tab = defaultdict(Counter)
    for g, d in per.items():
        if target in d:
            tab[band(sum(B.changed(c) for r, c in d.items() if r != target))][d[target]] += 1
    return {b: {"proteins": sum(c.values()), "lower": c["down"] / sum(c.values()), "raise": c["up"] / sum(c.values())}
            for b, c in sorted(tab.items())}


def main():
    import dataset_roles
    pred_file, out, a_file, b_file = sys.argv[1:5]
    tk, _ = S.read_classes(sys.argv[5]) if len(sys.argv) > 5 else S.read_classes()
    genes = A.marker_genes()
    lib = {r["gene"]: r["library"] for r in A.read_tsv(A.MAPS / "pois.tsv")}
    m, m10, pred = {}, {}, {}
    for half, path in (("A", a_file), ("B", b_file)):
        dataset_roles.require_role(f"lun2019_{half.lower()}", "fit")
        m.update(S.measured(path, genes, half))
        m10.update(S.measured(path, genes, half, timepoint="10"))
        pred.update(S.predictions(pred_file, half))
    calls = {k: B.call_of(e, d) for k, (e, d) in m.items() if k[0] != "AKT1"}
    per = by_protein(calls)
    low = sorted(g for g, d in per.items() if d.get(TARGET) == "down")
    falls = [(g, m[(g, TARGET)]) for g in low]
    rep = [B.call_of(*m10[(g, TARGET)]) == "down" for g in low if (g, TARGET) in m10]
    res = {"proteins_measured_at_target": sum(TARGET in d for d in per.values()), "lowering": len(low),
           "raising": sum(d.get(TARGET) == "up" for d in per.values()),
           "who": {"library": dict(Counter(lib.get(g, "?") for g in low)),
                   "library_all": dict(Counter(lib.get(g, "?") for g, d in per.items() if TARGET in d)),
                   "tyrosine_kinases_lowering": sum(g in tk for g in low),
                   "tyrosine_kinases_measured": sum(g in tk for g, d in per.items() if TARGET in d)},
           "by_activity": activity(per), "co_movement": co_movement(per),
           "size_in_delta": dict(Counter("1-1.5" if -e / d < 1.5 else "1.5-3" if -e / d < 3 else "3+" for _, (e, d) in falls)),
           "repeat_in_egf_well": {"n": len(rep), "share": sum(rep) / len(rep) if rep else None},
           "model_at_lowering": dict(Counter((pred[(g, TARGET)][0] if (g, TARGET) in pred else "not scored") + "/reach:" +
                                             (pred[(g, TARGET)][1] if (g, TARGET) in pred else "-") for g in low)),
           "lowering_proteins": low}
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
