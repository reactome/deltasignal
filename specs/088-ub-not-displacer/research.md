# specs/088 research: arm fsu

Solver `5783a3e`, `--server DS_CARRIER_LOOP_OPEN=fed_single_ub`, canonical
`20261021-2218_f99c4ab`, against `results/e4e7192/dflt086b`, results `results/5783a3e/fsu`.
Logs: 12,922 curator and 5,434 experimental release edge-solves dropped.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (6 / 6), all Signaling by EGFR | 1 |
| curator tuning | −4 (2 / 6) | 0.29 |
| experimental (RAF set aside) | 0 (2 / 2) | 1 |

By pathway:
- **Cell Cycle Checkpoints:** +2 curator, +2 experimental. MDM2 knockout at PolyUb-TP53
  (100× → 0.000001) and PolyUb-MDM2:MDM4 (100× → 0), both DOWN as expected.
- **EGFR:** 0 (6 / 6). SH3GL1 knockout and SRC overexpression gain; SRC knockout and EPS15
  overexpression lose.
- **Mitotic G1:** −6 curator, −2 experimental.

**The gate passes** (held-out 0, experimental 0, no pathway below −10). The net value is 0.

## Predictions

1. At least 1 of the 2 MDM2 experimental cases fixed: **holds**, both.
2. Only EGFR, Cell Cycle Checkpoints and Mitotic G1 move: **holds**.
3. No prediction for EGFR or held-out: EGFR 0, held-out 0.

## The Mitotic G1 loss (localised, not traced end to end)

All 8 changed G1 cases are at CCNA1 and CCNA2, the readouts that run through the
SCF(Skp2):Cks1 ligase:
- CDK4 and CCND1 knockout go from 0.84× (DOWN, right) to 14× and 2.7× (UP, wrong);
- CDKN2A and PPP2R1A knockout go from 18× and 14× (UP, right) to about 0.95× (NORMAL).

The rule drops SCF(Skp2)'s release after "Ubiquitination of p130 (RBL2) by SCF (Skp2)"
(R-HSA-1363331), so the ligase reads its supply. That supply is drained by the ubiquitin
pass's edges on free CKS1B, the defect specs/077 traced. specs/078's generator fix for it
lost the CBL-class cases. So the G1 loss is most likely a second, known defect exposed, not this rule
being wrong. specs/077 saw the same interaction.

## Reading

Correct where it acts (MDM2), net zero overall, because the SCF(Skp2) release it removes was
hiding the CKS1B drain. **Not recommended for adoption now.** It is worth revisiting together with a
fix for the CKS1B drain. That drain needs the ligase's own subunits excluded from the
ubiquitin pass without losing the CBL class (specs/078), and is not designed here.

## Decision

Adam, 2026-10-22: "go with your recommendation. we could revisit later as it has positive
attributes". **Not adopted.** `fed_single` stays the default; `fed_single_ub` stays reachable as
an opt-in, to re-measure once the CKS1B drain (ubiquitin pass on a ligase's own subunits) is fixed.

## Amendment 1: arm `fsu2` (2026-10-23)

Canonical `20261022-2011_233665c`, solver `b5beb78` (main) with
`--server DS_CARRIER_LOOP_OPEN=fed_single_ub`, against `results/b61b365/dflt095`. Results
`results/b5beb78/fsu2`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (6 / 6), Signaling by EGFR as in `fsu` | 1 |
| curator tuning | +2 (7 / 5) | 0.77 |
| experimental (RAF set aside) | **+1 (3 / 2)** | 1 |

- Fixed: both MDM2 knockout experimental cases in Cell Cycle Checkpoints (100× → 0×), ATM
  overexpression (a 50× rail → 0×, curator and experimental), 4 TP53 curator cases.
- Broken: CDK4 and CCND1 knockout in Mitotic G1 (0.82× → 0.88×, 0.80× → 0.86×, across the 0.85
  DOWN threshold). specs/089 removed the CKS1B drain, so the CDKN2A and PPP2R1A losses `fsu` showed
  are gone.

Predictions: (1) MDM2 knockout fixed: **holds**. (2) Mitotic G1 net ≥ −1 on each axis: **fails**
on experimental (−2), holds on curator (0). (3) Held-out: 0.

**The gate passes.** Adam, 2026-10-23: "ok", adopt. `fed_single_ub` is the default.
