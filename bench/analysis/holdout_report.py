#!/usr/bin/env python3
"""Report accuracy split into the tuning pathways and the held-out rest.

The MP-BioPath paper tuned on ten pathways and reported on the others. This
project drifted off that protocol: the ten-pathway set kept reversing decisions
that held at scale, so decisions moved to the wide curator set -- which is the
set we then report accuracy on. That fixed the misleading and silently turned
the test set into the training set.

This tool restores the split at REPORTING time, over an existing case dump, so
no benchmark has to be re-run. Use it for any number that leaves the building.

Two things it will show you:

  * the headline all-pathways figure is dragged down by the tuning ten, which
    are hard for everyone (TP53, WNT, PIP3, cell cycle);
  * whether a config decision made on the wide set also wins on pathways it
    never saw -- pass --compare with a second dump to check that directly.
    A decision that only wins on the tuning half is overfitting.
"""
from __future__ import annotations

import argparse
from collections import Counter
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import mcnemar_exact  # noqa: E402

# The paper's tuning set. Names as they appear in the curator files.
TUNING_PATHWAYS = {
    "Cell_Cycle_Checkpoints",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA_",
    "Mitotic_G2-G2_M_phases",
    "Mitotic_Prophase",
    "Mitotic_G1-G1_S_phases",
    "PIP3_activates_AKT_signaling",
    "RAF_MAP_kinase_cascade",
    "Signaling_by_ERBB2",
    "Signaling_by_WNT",
    "S_Phase",
    "Transcriptional_Regulation_by_TP53",
}

LABELS = ("0", "1", "2")


def load(path: Path) -> dict:
    with path.open(newline="") as fh:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(fh, delimiter="\t")
                if r["expected"] in LABELS}


def macro_f1(rows) -> float:
    f1s = []
    for lab in LABELS:
        tp = sum(1 for r in rows if r["predicted"] == lab and r["expected"] == lab)
        fp = sum(1 for r in rows if r["predicted"] == lab and r["expected"] != lab)
        fn = sum(1 for r in rows if r["predicted"] != lab and r["expected"] == lab)
        p = tp / (tp + fp) if tp + fp else 0.0
        rc = tp / (tp + fn) if tp + fn else 0.0
        f1s.append(2 * p * rc / (p + rc) if p + rc else 0.0)
    return sum(f1s) / len(f1s)


def report(name: str, cases: dict) -> None:
    groups = {"TUNING (the paper's ten)": [], "HELD-OUT (report this)": []}
    for k, r in cases.items():
        key = "TUNING (the paper's ten)" if k[0] in TUNING_PATHWAYS else "HELD-OUT (report this)"
        groups[key].append(r)
    print(f"\n=== {name} ===")
    print(f"{'split':<26}{'pathways':>9}{'cases':>8}{'accuracy':>10}{'macro-F1':>10}")
    for g, rows in groups.items():
        if not rows:
            continue
        pw = len({r["pathway"] for r in rows})
        ok = sum(1 for r in rows if r["predicted"] == r["expected"])
        print(f"{g:<26}{pw:>9}{len(rows):>8}{ok/len(rows):>10.4f}{macro_f1(rows):>10.4f}")
    allrows = [r for rs in groups.values() for r in rs]
    ok = sum(1 for r in allrows if r["predicted"] == r["expected"])
    print(f"{'(all pathways)':<26}{len({r['pathway'] for r in allrows}):>9}"
          f"{len(allrows):>8}{ok/len(allrows):>10.4f}{macro_f1(allrows):>10.4f}")
    # specs/025: cases whose entity is not in this Reactome release (or whose
    # gene name matches nothing) cannot measure the generator or the solver.
    # They stay in the rows above; these rows leave them out.
    if any(r.get("exclusion") for r in allrows):
        for g, rows in groups.items():
            kept = [r for r in rows if not r.get("exclusion")]
            if kept:
                ok = sum(1 for r in kept if r["predicted"] == r["expected"])
                label = g.split(" (")[0] + ", in release"
                print(f"{label:<26}{len({r['pathway'] for r in kept}):>9}{len(kept):>8}"
                      f"{ok/len(kept):>10.4f}{macro_f1(kept):>10.4f}")
        kept = [r for r in allrows if not r.get("exclusion")]
        ok = sum(1 for r in kept if r["predicted"] == r["expected"])
        print(f"{'(all), in release':<26}{len({r['pathway'] for r in kept}):>9}{len(kept):>8}"
              f"{ok/len(kept):>10.4f}{macro_f1(kept):>10.4f}")
        print("excluded:", dict(Counter(r["exclusion"] for r in allrows if r.get("exclusion"))))


def main(argv=None) -> int:
    # `argv` is for tests; the CLI passes nothing and argparse reads sys.argv.
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--compare", type=Path,
                    help="A second dump. Reports whether the change wins on the "
                         "HELD-OUT half as well as the tuning half -- a change "
                         "that only wins on tuning is overfitting.")
    a = ap.parse_args(argv)

    base = load(a.cases)
    report(a.cases.name, base)

    if a.compare:
        arm = load(a.compare)
        report(a.compare.name, arm)
        # Condition on the perturbation set, exactly as compare_wide_cases.py
        # does. A case key can be present in both arms while the two arms
        # resolved a different number of gene or knockout uuids for it -- then
        # the arms answered DIFFERENT questions and the pairing is invalid.
        # Unconditioned, a boundary-removal arm read +265; conditioned on the
        # same 7,168 cases it read +64. n_gene_uuids/n_ko_uuids are the
        # experiment's fingerprint.
        paired = set(base) & set(arm)
        shared = sorted(
            k for k in paired
            if base[k]["n_gene_uuids"] == arm[k]["n_gene_uuids"]
            and base[k]["n_ko_uuids"] == arm[k]["n_ko_uuids"]
        )
        dropped = len(paired) - len(shared)
        print(f"\n=== does it generalise?  {a.compare.name} vs {a.cases.name} ===")
        if dropped:
            print(f"  {dropped} of {len(paired)} shared cases dropped: the arms "
                  f"resolved a different perturbation set for them.")
        print(f"{'split':<26}{'cases':>8}{'baseline':>10}{'arm':>10}{'net':>8}"
              f"{'fixed':>7}{'broke':>7}{'p':>9}{'pathways':>10}{'readouts':>10}{'genes':>7}")
        concentrated = []
        for label, keys in (
            ("TUNING (the paper's ten)", [k for k in shared if k[0] in TUNING_PATHWAYS]),
            ("HELD-OUT (report this)", [k for k in shared if k[0] not in TUNING_PATHWAYS]),
        ):
            if not keys:
                continue
            b = sum(1 for k in keys if base[k]["predicted"] == base[k]["expected"])
            m = sum(1 for k in keys if arm[k]["predicted"] == arm[k]["expected"])
            # Discordant pairs and their exact McNemar p-value. A net figure
            # alone does not say whether the sign means anything: "-2 of 849"
            # was 0 fixed / 2 broke, p = 0.50, indistinguishable from a coin
            # flip, while "-16 of 18,808" was 11 fixed / 30 broke, p = 0.0043.
            # As percentages those look comparable. They are not.
            fixed = sum(1 for k in keys
                        if arm[k]["predicted"] == arm[k]["expected"]
                        and base[k]["predicted"] != base[k]["expected"])
            broke = sum(1 for k in keys
                        if base[k]["predicted"] == base[k]["expected"]
                        and arm[k]["predicted"] != arm[k]["expected"])
            pval = mcnemar_exact(fixed, broke)
            # HOW CONCENTRATED is the movement? McNemar assumes the discordant
            # pairs are independent. Cases sharing a pathway, a readout and a
            # mechanism are not: a coverage fix that makes ONE entity
            # addressable produces a cluster of "fixed" cases that reads as a
            # distributed gain and is not one.
            #
            # LNG #89 is the case that forced these columns. It reported
            # "+14 held-out, 16 fixed / 2 broke, p = 0.0013" -- which was all
            # 18 cases on a SINGLE readout (HSP90B1) in a single pathway, 9
            # genes x 2 directions. Effective n is 1, not 18, and the p-value
            # was meaningless.
            disc = [k for k in keys
                    if (arm[k]["predicted"] == arm[k]["expected"])
                    != (base[k]["predicted"] == base[k]["expected"])]
            n_pw = len({k[0] for k in disc})
            n_ro = len({(k[0], k[3]) for k in disc})
            # Distinct PERTURBATIONS, not just readouts. A loop fix once showed
            # 19 held-out discordant cases over 18 distinct readouts -- and 13
            # of them were one gene (DOK1) read at 13 places. Readouts alone
            # said "distributed"; genes said "two perturbations".
            gene_counts = collections.Counter((k[0], k[1]) for k in disc)
            n_ge = len(gene_counts)
            # Dominance, not just count: 19 discordant cases over 5 genes still
            # had 13 from ONE gene. The count said "five"; the share says "68%".
            top_gene, top_n = (gene_counts.most_common(1)[0] if gene_counts else (None, 0))
            top_share = top_n / len(disc) if disc else 0.0
            scored_pw = len({k[0] for k in keys})
            concentrated.append((label, n_pw, n_ro, n_ge, scored_pw, top_gene, top_share))
            print(f"{label:<26}{len(keys):>8}{b/len(keys):>10.4f}{m/len(keys):>10.4f}"
                  f"{m-b:>+8d}{fixed:>7}{broke:>7}{pval:>9.4f}"
                  f"{f'{n_pw}/{scored_pw}':>10}{n_ro:>10}{n_ge:>7}")
        print("\nA decision that wins on TUNING but not HELD-OUT is overfitting.")
        print("p is two-sided exact McNemar on the discordant pairs. p >= 0.05 means"
              "\nthe net figure's SIGN is not established, whatever its magnitude.")
        print("pathways = how many moved / how many scored. readouts = distinct"
              "\n(pathway, readout) pairs among the discordant cases. genes = distinct"
              "\n(pathway, gene) perturbations among them.")
        for label, n_pw, n_ro, n_ge, scored_pw, top_gene, top_share in concentrated:
            if n_ro == 0:
                print(f"\n  !! {label}: NO comparable case moved. Every discordant "
                      f"case was dropped by the\n     perturbation-set conditioning "
                      f"above, which is what happens when the change is a\n     "
                      f"COVERAGE fix -- it makes cases answerable that were not, so "
                      f"there is no\n     like-for-like pair. Score it as coverage "
                      f"(how many cases became answerable,\n     and how well they "
                      f"are answered), not as fixed/broke.")
                continue
            if top_share >= 0.5 and n_ro > 1:
                print(f"\n  !! {label}: one perturbation ({top_gene[1]} in {top_gene[0]}) "
                      f"accounts for {top_share:.0%} of the discordant cases\n     "
                      f"(across {n_ro} readouts). The readouts differ but the cause is "
                      f"shared, so the\n     independent evidence is closer to one case "
                      f"than to {n_ro}. Report the gene, not just the count.")
            if n_ro > 1 and n_pw > 1 and n_ge <= 2:
                print(f"\n  !! {label}: the movement spans {n_ro} readouts but only "
                      f"{n_ge} distinct (pathway, gene) perturbation(s).\n     Those cases "
                      f"share their cause, so McNemar's independence assumption fails "
                      f"even though\n     the readouts differ. Report it as a "
                      f"{n_ge}-perturbation result, and treat the p-value as overstated.")
            if n_ro <= 1 or n_pw <= 1:
                print(f"\n  !! {label}: the movement spans {n_pw} pathway(s) and "
                      f"{n_ro} readout(s).\n     These cases are NOT independent, so "
                      f"McNemar does not apply and the p-value above is\n     "
                      f"meaningless. Report this as a single-{'readout' if n_ro <= 1 else 'pathway'} "
                      f"fix, not a distributed gain.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
