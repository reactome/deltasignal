# specs/093 research: arms mvd (control) and mvk

Variant build `20261023-1115_55b5464_mv2` (canonical content plus the `moves` column), solver
`be1d73e`. `mvd` (defaults) reproduces canonical's defaults (`8162f5f/dflt091`) on every case
and value, both axes. `mvk` adds `--server DS_CARRIER_LOOP_MOVES=keep`. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **−10 (10 / 20)** | 0.099 |
| curator tuning | 0 | 1 |
| experimental (RAF set aside) | 0 | 1 |

Only Signaling by MET moves (L1CAM interactions and the GTSE1 pathway move no class).
**The gate fails** (held-out −10 < 0; MET −10).

## Predictions

1. Only MET, L1CAM and R-HSA-453274 move: **holds** (only MET moves a class).
2. At least 6 of the 10 CBL cases fixed: **holds**, all 10 (CBL knockout 1.0× → 33–100× UP,
   overexpression 1.0× → 0.13–0.78× DOWN, as the curators say).
3. No prediction for L1CAM or GTSE1: neither moves a class.

## The losses: HGF and MET knockout read UP

All 20 broken cases are HGF or MET knockout, at 10 readouts each: 0× (DOWN, right) → 16–100×
(UP). Both solves converge (residual 5.6e-7), so this is a stable state, not the iteration cap.
With the recycling release kept, the HGF:p-MET dimer at the membrane is fed by its own return
from the recycling endosome: dimer → adaptor complexes → endocytosis → RAB4 recycling → dimer.
That loop has gain at least 1. When HGF or MET is knocked out, its outside supply is 0 but the
loop holds itself up, and the readouts rail high. This is the loop behaviour specs/077 dropped
open releases to prevent: the same release that carries CBL's effect also closes a
self-sustaining receptor cycle.

## Reading

Not recommended. Keeping the recycling release restores CBL's 10 cases but breaks the
receptor's own 20 knockouts. A rule that fixes CBL without the rail would need the recycled
receptor to be bounded by its own supply (a receptor cannot recycle more than was made), as a
pool or a closed-twin bound does. That is a different rule, not registered here.

## Decision

Adam, 2026-10-23: "yes to option one", i.e. merge as a not-adopted record.
`DS_CARRIER_LOOP_MOVES` stays `drop`; the generator's `moves` column (branch
`carrier-moves-093`) is not merged, so canonical tables have no such column and `keep` would
act on nothing until it is.
