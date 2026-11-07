#!/usr/bin/env python3
"""specs/146 FR-004: any lookup difference fails (exit 1). An 'untouched'
(baseline-bin) case must come back null from the engine (amendment 1) and is
checked as exactly that; every other case must match confidence and stratum."""
import json
import sys
from pathlib import Path

here = Path(__file__).parent
fitter = json.loads((here / "fitter_conf.json").read_text())
engine = json.loads((here / "engine_conf.json").read_text())
assert len(fitter) == len(engine), (len(fitter), len(engine))
diffs = []
agree = untouched = 0
for f, e in zip(fitter, engine):
    if f["mag"] == "baseline":
        if e["confidence"] is None and e["stratum"] == "untouched":
            untouched += 1
        else:
            diffs.append(f"baseline case served {e}")
        continue
    if e["confidence"] != f["confidence"] or e["stratum"] != f["stratum"]:
        diffs.append(f"{f['call']}/{f['railed']}/{f['basis']}/{f['mag']}/{f['conv']}: "
                     f"fitter {f['confidence']}@{f['stratum']} engine {e['confidence']}@{e['stratum']}")
    else:
        agree += 1
print(f"{len(fitter)} lookups compared; agreeing {agree}, untouched-as-null {untouched}, "
      f"differences {len(diffs)}")
for d in diffs[:10]:
    print("  DIFF:", d)
if diffs:
    print("THE GATE FAILS")
    sys.exit(1)
print("no differences — the gate passes")
