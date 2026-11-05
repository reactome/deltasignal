#!/usr/bin/env python3
"""specs/138: did each korkut drug engage its target? Its direct marker on the array (named by mechanism below, before
the values of the drugs not yet looked at were read; AK, 901, P6, SR, ST had been) must move at |log2| >= log2 1.5 in the
expected direction at either dose. Drugs with no direct substrate on the array (HDAC, PKC) are `no marker`.

Usage: engagement.py > table.tsv
"""
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import korkut_bench as k  # noqa: E402

# drug code -> [(marker, expected sign)]: the target's own activation mark or its direct substrate
MARKERS = {
    "901": [("MAPKpT202", -1)],                                  # MEK -> ERK T202/Y204
    "AK": [("AKTpS473", -1), ("AKTpT308", -1)],                  # allosteric AKT inhibitor blocks both
    "HN": [],                                                    # HDAC: no acetylation readout
    "NT": [("p53", +1), ("p21", +1)],                            # MDM2 blocked -> p53 stabilised, p21 induced
    "P6": [("STAT3pY705", -1), ("STAT5pY694", -1)],              # JAK -> STAT Y705 / Y694
    "PLX": [("MEKpS217", -1), ("MAPKpT202", -1)],                # BRAF V600E -> MEK S217/221
    "RO": [],                                                    # PKC: no direct substrate readout
    "RY": [("RbpS807", -1)],                                     # CDK4 -> RB S807/811
    "SR": [("SRCpY416", -1)],                                    # SRC activation-loop autophosphorylation
    "ST": [("STAT3pY705", -1)],                                  # Stattic blocks STAT3 activation
    "Tm": [("p70S6KpT389", -1), ("S6pS235", -1), ("4EBP1pS65", -1)],  # mTORC1 substrates
    "ZS": [("AKTpS473", -1), ("AKTpT308", -1)],                  # PI3K -> PIP3 -> AKT
}
T = math.log2(1.5)


def table():
    names = [l.strip() for l in open(k.DATA / "node_Index.csv")]
    vals = {}
    for (i, ds), row in zip(k.conditions(), csv.reader(open(k.DATA / "expr.csv"))):
        if len(ds) == 1:
            vals.setdefault(ds[0], []).append(dict(zip(names, map(float, row))))
    out = []
    for d, marks in MARKERS.items():
        seen = [(m, s, [v[m] for v in vals[d]]) for m, s in marks]
        hit = any(s * x >= T for _, s, xs in seen for x in xs)
        out.append({"drug": d, "markers": ";".join(f"{m}:{'+' if s > 0 else '-'}:" + "/".join(f"{x:.2f}" for x in xs)
                                                    for m, s, xs in seen),
                    "engaged": "no marker" if not marks else "yes" if hit else "no"})
    return out


if __name__ == "__main__":
    w = csv.DictWriter(sys.stdout, ["drug", "engaged", "markers"], delimiter="\t")
    w.writeheader()
    w.writerows(table())
