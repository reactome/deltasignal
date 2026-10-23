# specs/090: a recycle pool reads its protein's supply, not its partners' mean

**Status:** pre-registered 2026-10-22, before the implementation and the arm. Arm `psa`: gate passes (held-out 0, experimental 0), but tuning −6 (0 / 6, p 0.03, all Mitotic G1) and prediction 3 refuted (RBL2 overexpression reads UP through DREAM repressing RBL1). Not recommended for adoption.

Adam, 2026-10-22: "ok" to candidate I1 of the 2026-10-22 review
(docs/reviews/2026-10-22-review.md, finding 2).

## The defect (review finding 2, traced on canonical)

`update_pool!` (src/core/cycles.jl) sets a pool's supply `s` to the plain mean of the fold of
every producer of every state, other than the pool's own steps. That is right for the states
of ordinary pools, whose several producers are OR alternatives (copies of one reaction). But a
specs/081 recycle pool's base states are complexes AND-assembled from the protein and its
partners: RBL1:E2F4:TFDP1 is fed by RBL1, E2F4 and TFDP1. The mean averages the protein with
its partners, so in Mitotic G1:
- an RBL1 knockout leaves the p107 pool at (0+1+1+0+1+1)/6 = 0.667× (it read 0 before specs/081,
  when the complex was an AND assembly);
- RBL1 overexpression at 100× gives s = 34, and p-p107 reads 17.9× instead of 51.6×.

## Census (`census.py` → `census.txt`, canonical `20261022-2011_233665c`)

States with more than one producer outside the pool's steps:
- ordinary pools: 17 states in 16 pools, every producer an OR `output` edge (62 edges): unchanged
  by the rule;
- recycle pools: 6 of 8 states in both pools (p107 `recycle1`, p130 `recycle2`, Mitotic G1),
  every producer an AND `assembly` edge (18 edges);
- release pools: none.

## The rule (`DS_POOL_SUPPLY_AND`, default `off`; `on`)

In the pool supply, a state's AND-joined producers count as **one** term: the product of
their folds, capped at 100× as the `hill_sat` AND is. Each OR producer stays its own term, as
now. `s` is the mean of the terms. A pool whose states have only OR producers (every ordinary
and release pool on canonical) is unchanged bit for bit. The solve reports `pool_supply_and`
and the number of AND groups it formed.

Biologically: the pool's resting complexes need all their subunits, so a knockout of the
protein (or of a partner) empties them, as the AND assembly did before specs/081.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `psa` | canonical `20261022-2011_233665c` | `--server DS_POOL_SUPPLY_AND=on` | `results/8a94f2e/dflt` (main, defaults) |

`dflt` is main's defaults on canonical: the review fixes changed no class, but they moved values,
so the baseline is re-run on main. Read with `holdout_report.py --set-aside
RAF_MAP_kinase_cascade`, case by case.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.
Mitotic G1 is a tuning pathway with no recycle pool elsewhere, so held-out cannot move, and the
decision rests on faithfulness and on the tuning/experimental cases, reported both ways.

## Predictions

1. Only Mitotic G1 moves (the only pathway with recycle pools).
2. RBL1 and RBL2 knockout: the p107/p130 pool forms read 0, and the E2F-target readouts move
   further UP. The curator cases (48 of 50 right) stay right. The 20 experimental cases stay
   wrong: the experiments read NORMAL (p107 and p130 back each other up), the curators UP. This
   rule does not address that disagreement.
3. RBL1 and RBL2 overexpression read lower at the E2F targets (stronger repression). No
   prediction on their net (45 of 50 right).
4. Net on all axes within ±5. **The rule is proposed on faithfulness, not on a gain.**
