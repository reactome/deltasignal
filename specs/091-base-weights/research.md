# specs/091 research: arm pbw

Canonical `20261022-2011_233665c`, solver `318c3cf` with `--server DS_POOL_BASE_WEIGHT=split`,
against `results/8a94f2e/dflt` (main defaults). Results `results/318c3cf/pbw`; 82 multi-base
pool-solves split. Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade
--pair-on-valid`; cases with `specs/089-ub-ligase-unit/diff_cases.py`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (0 / 0) | 1 |
| curator tuning | **+35 (40 / 5)** | < 1e-4 |
| experimental (RAF set aside) | +1 (2 / 1) | 1 |

By pathway: **Mitotic G1 +37** (curator), +2 (experimental); **WNT −2** (curator), −1
(experimental). Nothing else moves; L1CAM interactions has no scored change. Converged solves
unchanged (1,724 / 1,725, 244 / 244).

**The gate passes** (held-out 0, experimental +1, worst pathway −2).

## Predictions

1. Only Mitotic G1, WNT and L1CAM move: **holds** (L1CAM moves no class).
2. APC/AMER1 knockout no longer rail at β-catenin targets, class unchanged: **holds at the
   targets** (100× → 10×, all still UP), but the cytosolic β-catenin readout itself flips (below).
3. The p107/p130 kinase and phosphatase perturbations move further: **holds**, and in the right
   direction: CCND1 knockout 5, CDK4 knockout 6, CDKN2A knockout 8 and overexpression 2, PPP2R1A
   knockout 9 and overexpression 9 fixed (e.g. CDK4 knockout at CCNA2 0.90× → 0.73×, DOWN as the
   curators say). PPP2R1A at 8848441 breaks twice (threshold, 1.07× → 1.21×, 0.94× → 0.84×).
4. No net predicted: tuning +35.

## The losses

- **APC and AMER1 knockout at cytosolic β-catenin (R-HSA-448839)**, curator and experimental
  (APC): 9.09× → 0.91×. The readout is the mean of 44 copies of the entity. Under the
  knockouts 40 of them read 0 and the 4 pool copies rise (`trace_apc.jl`): to 100× under `full`
  (the 181× ceiling of 20 bases, clamped), to 10× under `split`. 4 × 100 / 44 = 9.09 (UP),
  4 × 10 / 44 = 0.91 (NORMAL). The case was right only because the 4 pool copies railed at the
  clamp. Why 40 copies of free β-catenin read 0 under an APC knockout is not traced here.
- **TCF7L2 overexpression at 2130284**: 0.95× → 0.50× (curators NORMAL), while the same
  perturbation at 448839 is fixed (1.80× → 0.95×, NORMAL).
- **PPP2R1A at 8848441**: two threshold shifts.

## Reading

Recommended for adoption. The rule removes a dependence of pool magnitudes on variant
multiplicity, gains 37 Mitotic G1 cases, all CDK4/6 and PP2A perturbations of the p107/p130
pools in the curated direction, and loses 3 WNT cases, 2 of which were right only through the
20-copy ceiling.

## Decision

Adam, 2026-10-23: "lets go with option 1", i.e. adopt. `DS_POOL_BASE_WEIGHT=split` is the
default; `full` restores specs/081's weights. The specs/081 test fixture keeps its full-weight
numbers by setting `full` explicitly.
