# specs/065 research: measured 2026-10-16

Variant build `20261016-0924_8efc1b1_sigcache` (LNG `8efc1b1` = main
`da08370` + the fix, LNG PR #115) against canonical `20261016-0006_da08370`.

## Networks

- Same totals: 211,767 nodes, 500,762 edges; every CSV in every pathway has
  the same row count.
- Files differ byte for byte only because every build mints new uuids. With
  each node labelled by its `nodes.csv` content (kind, diagram entity,
  compartment, member leaves, source sets, chosen members), the edge
  multiset of **all 92 pathways is identical**; no edge endpoint lacks a node
  row.

## Predictions

Solver main `41f3a5a`, `results/41f3a5a/sigcache` against `results/0406a4b`:

| axis | rows | predicted class changed | value changed |
|---|---|---|---|
| curator | 24,100 | 0 | 15 |
| experimental | 849 | 0 | 0 |

The 15 value changes are all in solves that do not converge (DSB Repair:
ATM, RAD50, MDC1, KAT5, NBN, RNF8, MRE11, KPNA2, RNF168 overexpression;
EGFR: SH3GL1 ×5, PLCG1 knockout), at most 0.44% relative. With identical
content they come from the re-minted uuids: tied positions in the sweep
keep label order (CLAUDE.md, `DS_SCC_ORDER=flow`), which only shows where
the iteration stops at the cap. Headline 20,894 / 24,100 in both.

## Decision (by the registered rule)

Nothing moves: the stale-signature path did not occur in this build order,
so the fix changes no network today and guards the next pathway list or
build order. Merge LNG #115; no new canonical build.
