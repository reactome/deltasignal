#!/usr/bin/env python3
"""specs/131: run a bench with the surrogate in place of the solver.

Usage: run_surrogate_131.py korkut OUT_DIR
       run_surrogate_131.py mpbiopath [benchmark_vs_mpbiopath.py arguments]
The server named by DELTASIGNAL_BASE still parses the networks; every solve is the surrogate's.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import surrogate_131  # noqa: E402

surrogate_131.install(bm)
which, *rest = sys.argv[1:]
if which == "korkut":
    import korkut_bench
    sys.argv = ["korkut_bench.py", *rest]
    korkut_bench.main()
elif which == "mpbiopath":
    sys.argv = ["benchmark_vs_mpbiopath.py", *rest]
    bm.main()
else:
    raise SystemExit("first argument: korkut or mpbiopath")
