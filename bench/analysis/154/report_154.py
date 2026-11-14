#!/usr/bin/env python3
"""specs/154 FR-006: the coverage report — does the detector explain the
recorded rails, and do the two traced cases pass? Serving bar: both named
checks PASS; the rails share is reported whatever it is.

Usage: report_154.py SWITCH_CASES_TSV
"""
import csv, sys
from collections import defaultdict

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
sb = lambda v: v == "True"

scored = [r for r in rows if r["fold_served"] not in ("", "None")]
mismatch = sum(1 for r in scored if r["fold_recorded"] not in ("", "None")
               and abs(float(r["fold_served"]) - float(r["fold_recorded"])) > 5e-4 * max(1.0, float(r["fold_recorded"])))
print(f"cases {len(rows)}, served-valued {len(scored)}, served-vs-recorded fold mismatches {mismatch}")

railed = [r for r in scored if sb(r["railed_recorded"])]
flagged = [r for r in railed if sb(r["switch"])]
print(f"\nrailed recorded folds covered: {len(railed)}; flagged as switch: {len(flagged)} "
      f"({100*len(flagged)/max(1,len(railed)):.1f}%)")
by_pw = defaultdict(lambda: [0, 0])
for r in railed:
    by_pw[r["pathway"]][0] += 1
    sb(r["switch"]) and (by_pw[r["pathway"]].__setitem__(1, by_pw[r["pathway"]][1] + 1))
print("per pathway (railed, flagged), worst coverage first:")
for pw, (n, f) in sorted(by_pw.items(), key=lambda x: x[1][1]/max(1,x[1][0]))[:12]:
    print(f"  {pw[:48]:48} {f}/{n}")

# --- the serving bar ---
RAS_RAF = {"HRAS", "KRAS", "NRAS", "ARAF", "BRAF", "RAF1"}
raf = [r for r in railed if r["pathway"] == "RAF_MAP_kinase_cascade" and r["gene"] in RAS_RAF]
raf_ok = raf and all(sb(r["switch"]) for r in raf)
states_ok = True
for r in raf:
    if not sb(r["switch"]):
        continue
    hi = float(r["fold_served"]); lo = float(r["other_state_fold"])
    if r["state"] == "low":
        hi, lo = lo, hi
    if not (hi >= 20.0 and lo <= 1.5):
        states_ok = False
        print(f"  RAF state off-record: {r['gene']} dir{r['direction']} ko {r['key_output']}: "
              f"served {r['fold_served']} / other {r['other_state_fold']} ({r['state']})")
print(f"\nRAF check (specs/054): {len(raf)} railed RAS/RAF MAPK folds, all flagged: {raf_ok}, "
      f"states match the record (high >=20x, low <=1.5x): {states_ok} -> "
      f"{'PASS' if raf_ok and states_ok else 'FAIL'}")

pten = [r for r in scored if r["pathway"] == "PIP3_activates_AKT_signaling" and r["gene"] == "PTEN"
        and sb(r["railed_recorded"])]
pten_ok = bool(pten) and any(sb(r["switch"]) for r in pten)
print(f"PTEN/MKRN1 check (specs/082): {len(pten)} railed PTEN folds in PIP3, "
      f"flagged: {sum(1 for r in pten if sb(r['switch']))} -> {'PASS' if pten_ok else 'FAIL'}")

print(f"\nSERVING BAR: {'PASS' if raf_ok and states_ok and pten_ok else 'FAIL'}")
