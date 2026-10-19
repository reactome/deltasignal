# specs/076 research: results of arms rst, rstd and rstd2k

Solver `bf4b688`, run 2026-10-19 with `scripts/run_arm.sh`, `DS_NONCONV_RESTARTS=4`.

| arm | build | against | results |
|---|---|---|---|
| rst | canonical `20261017-1246_5e663d9` | `results/3177881/rvc` | `results/bf4b688/rst` |
| rstd | faithful `20261018-0815_ddbd828_dsub` (+ specs/071, 072 rules) | `results/3177881/rvcd` | `results/bf4b688/rstd` |
| rstd2k | as rstd, `DS_MAX_ITERS=2000` | `rstd` | `results/bf4b688/rstd2k` |

## rst: canonical

- **No class changes on either axis.** Curator 23,511 paired valid cases, experimental 845.
  `holdout_report.py` (RAF set aside or not) finds no discordant pair. **The
  registered gate passes: held-out 0, experimental 0, no pathway moves.**
- 10 curator `pred_ui` values move by more than 1e-3, all in DSB Repair (9) and DSB Response (1). None crosses a class boundary.
- Converged curator cases rise from 23,117 to 23,284 (+167, all newly converged, none lost). Experimental stays at 824.
- Restarts: curator 30 used and 10 rescued; experimental 0 used.

## rstd: faithful build, against rvcd

Only Transcriptional Regulation by TP53 changes class (277 curator, 25
experimental). Every TP53 case now converges: non-converged goes from 306 to **0**. Curator
converged goes from 22,973 to 23,428, and experimental from 811 to 839. Restarts: curator 22 used and 17 rescued; experimental 9 used and 8 rescued.

TP53 curator cases correct: rvcd 1,466, rstd **1,384** (net −82: 56 fixed,
138 broken). Experimental: 166 to 162 (net −4: 5 fixed, 9 broken).

All 302 changes are in the six perturbations that rvcd left non-converged.
The two AKT knockouts were already non-converged in `pljd`; the other four
are the ones specs/075 lost.

| perturbation | curators | rvcd (cap answer) | rstd (fixed point) | canonical | fixed / broken |
|---|---|---|---|---|---|
| AKT1 knockout | UP 48 | mostly UP | NORMAL 49 | NORMAL 51 | 0 / 47 |
| AKT2 knockout | UP 48 | mostly UP | NORMAL 49 | NORMAL 51 | 0 / 47 |
| MDM4 knockout | UP 48 | mostly UP | NORMAL 44 | NORMAL 46 | 0 / 41 |
| MDM4 overexpression | DOWN 48 | mostly UP | DOWN 46 | DOWN 44 | 47 / 0 |
| CHEK2 knockout | DOWN 41 | mostly UP | NORMAL 38, DOWN 7 | DOWN 47 | 9 / 0 |
| DAXX overexpression | NORMAL 51 | mostly UP | DOWN 46 | NORMAL 51 | 0 / 3 |

(Readouts are counted by class. "Mostly UP" means UP was the largest class in
rvcd's non-converged answer.)

## rstd vs rstd2k

**No class changes on either axis.** At cap 2,000, 188 curator and 18
experimental cases end non-converged, and 5 curator cases newly converge.
Every TP53 case converged in both arms reads the same class. 32 curator
`pred_ui` values move by more than 1e-3.

## Run time

These are the wall times for both axes. rst and rstd ran in parallel, as rvc and rvcd did.

| arm | time | reference | time |
|---|---|---|---|
| rst | 24 min | rvc | 26 min |
| rstd | 29 min | rvcd | 30 min |
| rstd2k | 26 min | | |

Restarts add no measurable time.

## Predictions

1. In `rst`, class changes occur only in solves that were non-converged, and converged solves do not decrease: **holds**. There are no class changes, and 167 more cases converge.
2. In `rstd`, at least 3 of the 4 perturbations converge, and TP53 has at most 102 non-converged cases: **holds**. All four converge, and so do both AKT knockouts, leaving 0 non-converged.
3. Between `rstd` and `rstd2k`, no TP53 case that converged in both changes class: **holds**.
4. Run time: reported above; unchanged.

## Reading

- On canonical, the rule changes no class, and it converges 167 more DSB Repair cases.
- On the faithful build, TP53 now reaches its fixed point under every
  perturbation. The fixed point is worse against the curators than the
  cap's answer (−82). The cap had read AKT1, AKT2 and MDM4 knockout UP, and the curators agree with UP.
  - At the fixed point, those three read NORMAL.
  - Canonical, which converges, also reads them NORMAL.
  - In `pljd`, before specs/075, the AKT knockouts were non-converged and read
    NORMAL. In `rvcd` they were still non-converged, now reading UP. 98 of the 186 AKT
    cases that specs/075 counted as fixed were therefore the cap's answer, not
    the model's.
- Where the faithful build now differs from canonical at a fixed point:
  - CHEK2 knockout reads NORMAL where canonical and the curators read DOWN.
  - DAXX overexpression reads DOWN where both read NORMAL.

  In `pljd` (before specs/075), CHEK2 knockout converged and read DOWN (38
  correct), so specs/075's carrier rule is what moves it. These are the two
  TP53 differences left to trace.
- Faithful TP53: 1,384 of 1,836 correct (canonical 1,607).

## Adoption note (2026-10-19)

When the default became 4, `test/test_loop_elasticity.jl` failed. Its "knife-edge
at eps = 1" test pins the specs/014 defect, where the iteration budget picks the
answer, as a climb with the budget (`a50 < a500 < a5000`). The fixture is a
positive AND loop driven by U = 1.01. Probed values of A (UI):

| cap | 50 | 500 | 750 | 1,000 | 5,000 |
|---|---|---|---|---|---|
| restarts 0 | 0.058, no | 0.141, no | 0.015, yes | 0.025, yes | 100, yes |
| restarts 4 | 0.017, yes | 0.015, yes | 0.015, yes | 0.025, yes | 100, yes |

The specs/066 fallback already drew this loop to near zero from the cap-750
iterate, where the absolute tolerance (1e-6 on the 0.01 baseline) accepts a
state that is not a fixed point (the loop's only fixed points are 0 and the
cap). Restarts bring the same answer to caps 50 and 500. The defect stands,
with a new shape: the budget now chooses between collapse and rail. The test
pins both shapes, the restart-free climb under `DS_NONCONV_RESTARTS=0`
(175 assertions).
