# specs/082: how many cases read the same off-baseline value under knockout and
# overexpression (a switch, not a response), and how many of them are wrong.
# Usage: python3 switch_census.py <results dir>
import csv, sys
from collections import Counter, defaultdict
R = sys.argv[1]
for ax in ("curator", "experimental"):
    rows = [r for r in csv.DictReader(open(f"{R}/{ax}_cases.tsv"), delimiter="\t") if r["valid"] == "1"]
    by = defaultdict(dict)
    for r in rows: by[(r["pathway"], r["gene"], r["key_output"])][r["direction"]] = r
    sw = []; 
    for k, d in by.items():
        if "0" in d and "2" in d:
            a, b = float(d["0"]["pred_ui"]), float(d["2"]["pred_ui"])
            if abs(a - 1) > 0.15 and abs(b - 1) > 0.15 and (a > 1) == (b > 1):
                sw.append((k, d))
    cases = [x for _, d in sw for x in d.values()]
    wrong = [r for r in cases if r["predicted"] != r["expected"]]
    allwrong = sum(r["predicted"] != r["expected"] for r in rows)
    print(f"{ax}: {len(sw)} (pathway, gene, readout) pairs read the same direction under KO and OE; "
          f"{len(cases)} cases, {len(wrong)} wrong (of {allwrong} wrong in all)")
    print("   by pathway (wrong):", Counter(r["pathway"] for r in wrong).most_common(8))
    print("   by gene (wrong):", Counter((r["pathway"][:20], r["gene"]) for r in wrong).most_common(8))
    # a rail reads the SAME value both ways (PTEN 83.5 / 83.5); a U-shape reads two values
    rail = [(k, d) for k, d in sw if abs(float(d["0"]["pred_ui"]) - float(d["2"]["pred_ui"])) < 1e-3]
    rw = [x for _, d in rail for x in d.values() if x["predicted"] != x["expected"]]
    print(f"   of these, identical value both ways (a rail): {len(rail)} triples, {len(rw)} wrong;",
          Counter((r["pathway"][:20], r["gene"]) for r in rw).most_common(6))
