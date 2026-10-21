# specs/085: silo bridges re-measured on today's networks (replication)

**Status:** pre-registered 2026-10-21, before the arm. No code change.

Adam, 2026-10-21: "continue", on replicating before closing the two-copy direction.

## Why

specs/080 and its traces showed that the "two-copy" breaks (e.g. phospho-TP53 tetramer made
by ATM in one copy, read by the CASP1-gene binding from a separate starting copy) are the
positional silo that `src/core/silo_bridges.jl` bridges. Silo bridges were measured four
times and lost: at cap 50, −73 and −77 on two catalogs (September), with 126 of 146
losses being NORMAL turned into a change call. Since then the networks and solver have
changed a lot (variant nodes, pools, carrier rules, recycle pools, component twins). This
is one replication at the same cap, not a sweep.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `silo50` | canonical `20261020-2038_23ef3e1` | main `b28d5d5`, `--server DS_SILO_BRIDGE_MAX_REACH=50` | `results/1b52a97/cct2` (= main's defaults) |

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`, on both axes.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10.

## Predictions

1. **The gate fails:** held-out net negative.
2. Most losses (at least 60%) are NORMAL turned into a change call, as in September.
3. If it fails, the two-copy direction is closed as an accuracy lever on today's networks,
   and the next work is the switch loops (specs/082).
