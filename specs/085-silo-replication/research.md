# specs/085 research: arm silo50

Solver `ae1a74b` (main `b28d5d5` + this spec), `--server DS_SILO_BRIDGE_MAX_REACH=50`,
canonical `20261020-2038_23ef3e1`, against `results/1b52a97/cct2`, results
`results/ae1a74b/silo50`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **−90** (44 / 134) | < 1e-4 |
| curator tuning | +48 (67 / 19) | < 1e-4 |
| experimental (RAF set aside) | +1 (3 / 2) | 1 |

By pathway (curator, 27 moved):
- **Losses:** AP-2 −22 (CREBBP, EP300), intrinsic apoptosis −17 (34 broken / 17 fixed), RUNX2
  −14 (PPM1D, ESR1), ROBO −14, Base Excision Repair −10, IGF2BP −8, S phase −8, IQGAPs −6.
- **Gains:** TP53 +32, Mitotic G1 +25, GPVI +10 (SYK, PTPN11, the case silo bridges were
  built for), Insulin receptor +6.

**The gate fails** (held-out −90).

## Predictions

1. The gate fails: **holds**.
2. At least 60% of the losses are NORMAL turned into a change call: **holds**, 112 of 153
   (73%). The gains are mostly the same transition (100 of 111).
3. The two-copy direction is closed: **yes**.

## Reading

This is the fifth silo measurement and the first on today's networks. The result is the same
in sign and in shape as September's −73 and −77. Bridging copies is a machine for change
calls, and the curators call most of these cases NORMAL. The tuning pathways gain (TP53,
Mitotic G1) while held-out loses: the bridges fit the pathways traced most and generalise
worst. **The positional silo is not an accuracy lever on these networks**, and the "two-copy"
breaks of specs/080 are closed with it. Individual copies can still be wrong in specific,
traced cases (as specs/081's p130 and specs/084's ATM were); those are fixed one mechanism
at a time, never by bridging wholesale.
