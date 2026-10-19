# specs/075 research: results of arms rvc and rvcd

Solver `3177881`, run 2026-10-18 with `scripts/run_arm.sh`.

| arm | build | against | results |
|---|---|---|---|
| rvc | canonical `20261017-1246_5e663d9` | `results/501e480/hsc` | `results/3177881/rvc` |
| rvcd | faithful `20261018-0815_ddbd828_dsub` (+ specs/071, 072 rules) | `results/ae9c56f/pljd` | `results/3177881/rvcd` |

## rvc: canonical

Identical to `hsc` on every case, both axes, down to `pred_ui`; the log
reports 0 pool-carrier pairs. **The registered gate passes, trivially.**

## rvcd: faithful build, against pljd

Only Transcriptional Regulation by TP53 moves (941 curator and 151
experimental readouts). The log reports 538 pool-carrier pairs over the
solves.

| | pljd | rvcd |
|---|---|---|
| TP53 curator correct (of 1,836) | 1,214 | **1,466** (canonical: 1,607) |
| curator, net | | **+252** (294 fixed, 42 broken) |
| experimental, net | | **+15** (20 fixed, 5 broken) |
| TP53 non-converged curator cases | 102 | 306 |
| converged solves, curator / experimental | 1,705 / 241 | 1,701 / 237 |

TP53 curator cases by perturbation (fixed, broken, newly non-converged):

| perturbation | fixed | broken | lost convergence |
|---|---|---|---|
| AKT1, AKT2 knockout and overexpression | 186 | 0 | 0 |
| MDM2 overexpression | 49 | 1 | 0 |
| MDM4 knockout | 48 | 0 | 51 |
| ATM, BRCA1, DAXX, MDM4 overexpression | 11 | 3 | 102 |
| CHEK2 knockout | 0 | 38 | 51 |

Four perturbations lose convergence (CHEK2 knockout, DAXX overexpression,
MDM4 knockout and overexpression; residuals 0.006 to 0.034). Their 48 fixes
and 38 losses are where the iteration stopped, not answers.

## Predictions

1. `rvc` identical to canonical: **holds**.
2. `rvcd`: only TP53 moves, and at least 30 of the 46 MDM2 overexpression
   cases read DOWN or NORMAL: **holds**. All 46 read DOWN, all converged.
3. `rvcd`: TP53 non-converged within 50 of pljd's 102: **fails**, 306.

## Reading

- The traced defect was the cause: with USP7 read by the carrier rule in the
  drive, MDM2 overexpression reads TP53's targets DOWN on the faithful build,
  as on canonical and as the curators do. AKT1 and AKT2 (186 cases) go the
  same way.
- Excluding the four perturbations that lose convergence, 241 cases are
  fixed and 3 broken.
- The faithful build's TP53 is now 1,466 / 1,836 against canonical's 1,607.
- Four perturbations newly fail to converge, and the joint LM (specs/073)
  does not rescue them. That is the next thing to trace in TP53.
