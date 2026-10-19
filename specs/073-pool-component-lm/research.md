# specs/073 research: results of arms plj, pljd, pljd2k and control crsd74

Solver `ae9c56f` (includes the specs/074 fix), run 2026-10-18 with
`scripts/run_arm.sh`.

| arm | build | overrides | results |
|---|---|---|---|
| plj | canonical `20261017-1246_5e663d9` | `DS_NONCONV_POOLS=joint` | `results/ae9c56f/plj` |
| crsd74 | faithful `20261018-0815_ddbd828_dsub` | `DS_POOL_SUPPLY_REACH=material`, `DS_CYCLE_CARRIER_RELEASE=reaction` | `results/ae9c56f/crsd74` |
| pljd | faithful | as crsd74, + `DS_NONCONV_POOLS=joint` | `results/ae9c56f/pljd` |
| pljd2k | faithful | as pljd, + `DS_MAX_ITERS=2000` | `results/ae9c56f/pljd2k` |

Note on reading `holdout_report.py`: its "baseline" column is the `--cases`
file and "arm" is `--compare`. Every figure below is signed arm minus control
and was checked case by case (scratch `cmp074.py`).

## plj: canonical

No pool component is non-converged on canonical, so nothing is re-solved
(log: 0 re-solved). Every case is identical to `hsc` (`results/501e480/hsc`)
down to `pred_ui`. **The registered gate passes, trivially.**

## pljd: faithful build, against crsd74

| | crsd74 | pljd |
|---|---|---|
| pool components skipped / re-solved, curator | 11 / 0 | 0 / 11 (11 improved, 11 rescued) |
| pool components skipped / re-solved, experimental | 11 / 0 | 0 / 11 (11 improved, 11 rescued) |
| converged solves, curator | 1,694 / 1,725 | 1,705 / 1,725 |
| converged solves, experimental | 230 / 244 | 241 / 244 |
| TP53 non-converged curator cases | 663 | 102 |

Classes (RAF set aside): held-out +0; tuning **−11** (1 fixed, 12 broken,
p 0.0034); experimental +0 (3 class changes, none between right and wrong).
Only TP53 moves.

All 36 curator class changes were non-converged in crsd74 and are converged
in pljd:

| perturbation | changes | fixed | broken | crsd74 → pljd (expected) |
|---|---|---|---|---|
| MDM2 overexpression | 10 | 0 | 6 | DOWN → UP ×6, NORMAL → UP ×4 (DOWN) |
| MDM4 overexpression | 7 | 1 | 6 | DOWN → UP ×6 (DOWN) |
| AKT1, AKT2 overexpression | 12 | 0 | 0 | NORMAL → UP (DOWN) |
| DAXX overexpression | 7 | 0 | 0 | DOWN → UP ×6 (NORMAL) |

The 12 broken cases were right only where the cap stopped the oscillation.
At the fixed point the faithful network reads them UP.

Against canonical (`hsc`), pljd is held-out −171 and tuning −306, as before:
the gap to canonical is network construction (specs/070), not convergence.

## Predictions

1. `plj` identical to canonical, 0 re-solved: **holds**.
2. Only TP53 moves: **holds**, both axes.
3. At least three quarters of the skipped components improved; TP53
   non-converged at most half of crsd74's: **holds** (11/11 on each axis;
   663 → 102).
4. No TP53 case converged in both pljd and pljd2k changes class: **holds**.
   6 cases change class (MDM2 overexpression); all are converged in pljd and
   not in pljd2k.
5. At most 15 of the 46 MDM2 overexpression cases read DOWN: **holds**,
   0. All 46 converge and read UP (crsd74: 36 UP, 6 DOWN, 4 NORMAL).

## pljd2k: the cap still matters, through the start point

At cap 2,000, fewer solves converge than at cap 500: curator 1,696 against
1,705, and 257 cases converged in pljd are not converged in pljd2k (28 the
other way). Of 6 pool components re-solved, 4 are rescued. LM starts from
the last iterate, and from the cap-2,000 iterate it stalls in a local minimum
more often (the specs/073 probe showed the same for one case). Where both
converge, the classes agree.

## Reading

- On canonical the rule changes nothing.
- On the faithful build it does what it was registered for: TP53 now reaches
  a fixed point in 561 more curator cases. The −11 is the cap's lucky
  answers removed, not a new error.
- At that fixed point, MDM2 and MDM4 overexpression read TP53's targets UP
  where curators read DOWN. The specs/072 trace places the cause in the PLK3
  loop (p-S15,S20-TP53 induces PLK3, which phosphorylates TP53). That is now
  a network question about the faithful build, and it can be traced on
  converged solves.
- The start-point dependence (pljd2k) is a remaining limitation of a local
  root finder.
