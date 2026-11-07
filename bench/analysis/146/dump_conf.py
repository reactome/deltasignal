#!/usr/bin/env python3
"""specs/146 FR-004, fitter side: every validation case's features and the
fitter's own lookup answer, for conf_146.jl to replay through the engine."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parent))
from fit_confidence import cases, load_basis, fit, lookup, TUNING_PATHWAYS  # noqa: E402

rd = Path(sys.argv[1])
basis = load_basis(HERE.parents[1] / "141" / "bench_readouts.json")
curator = cases(rd / "curator_cases.tsv", basis)
experimental = cases(rd / "experimental_cases.tsv", basis)
counts = fit([c for c in curator if c["pathway"] in TUNING_PATHWAYS])

out = []
for split, cs in (("held_out", [c for c in curator if c["pathway"] not in TUNING_PATHWAYS]),
                  ("experimental", experimental)):
    for c in cs:
        conf, lv = lookup(counts, c)
        # the engine's vocabulary for the normal class, so the replay
        # crosses the same seam the handler does
        engine_call = "normal" if c["call"] == "not" else c["call"]
        out.append({"split": split, "call": engine_call, "railed": c["railed"],
                    "basis": c["basis"], "mag": c["mag"], "conv": c["conv"],
                    "confidence": round(conf, 4), "stratum": lv})  # the table serves 4dp
Path(HERE.parent / "fitter_conf.json").write_text(json.dumps(out))
print(f"{len(out)} validation lookups dumped")
