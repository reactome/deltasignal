# specs/073: re-solve a non-converged pool component jointly

**Status:** drafted 2026-10-18. Blocked on specs/074 (a numerical defect in
the default AND, found by this spec's probe); the arms wait for it. The
predictions below are revised from one probed case and say so.

Adam, 2026-10-18: "ok" to #2 of specs/071 research ("LM fallback for pool
components"), after the specs/072 report.

## The defect

The specs/066 fallback (`DS_NONCONV_FALLBACK=lm`, default) re-solves a cyclic
component that is still above tolerance at the iteration cap by
Levenberg-Marquardt on |F(x) − x|². It skips any component that holds
cycle-balance pool nodes (`nonconv_skipped_pools`). Those nodes are written by
their pool's balance (π = πP), not by F, so the specs/066 residual does not
describe them.

- Canonical `20261017-1246_5e663d9`, solver `6fb4329`: 0 component-solves
  skipped on either axis. The skip changes nothing there.
- Faithful build `20261018-0815_ddbd828_dsub` under arm `crsd`
  (`results/d1d8bb9/crsd`): 12 component-solves skipped on each axis. TP53
  has 663 non-converged curator cases.

**Traced, before this spec** (`~/deltasignal-catalogs/trace_073/`, scratch
`capsweep.jl`): MDM2 overexpression in Transcriptional Regulation by TP53
under `crsd`, solved at caps 500 to 4,000. The component never converges
(residual 0.03 to 0.2). The readout 5628829 (output `be474bf4`) reads 0.33×,
0.34×, 0.36×, 0.39×, 0.47×, 1.00×, 0.33×, 0.47×, 7.47× and 0.38× at caps
500, 501, 502, 503, 505, 510, 550, 1,000, 2,000 and 4,000. PLK3 swings from
0.70× to 9.6× and p-S15,S20-TP53 from 0.51× to 8.6×. The class is the cap's
choice, as in specs/066 (an oscillation, not slow convergence).

**Probed, before the arms** (same case, `DS_NONCONV_POOLS=joint`):

- On the solver as it stands, LM improves the residual (0.03 → 0.0017) but
  stops on a vanishing step. The Jacobian is wrong along row 9, an AND of
  node 3451 (1.21×) and pinned MDM2 (80×): the exact gradient is 70.8, the
  finite difference 0. The cause is the `hill_sat` cap, whose output between
  87.5× and 99.8× is cancellation noise (specs/074). LM iterations 400, 2,000
  and 10,000 end at the same point.
- With the specs/074 fix in the working tree:
  - `skip` still oscillates (residual 0.03 to 0.2; readout 0.33×, 0.47×,
    7.47× at caps 500, 1,000, 2,000). The AND defect is not the oscillation.
  - `joint` converges (residual 5e-7) from the cap-500 and cap-1,000
    iterates, to the same state. The readout reads **1.76× (UP)**, PLK3
    2.84×, p-S15,S20-TP53 2.84×. From the cap-2,000 iterate LM stalls
    (residual 0.011, readout 4.46×).

## The rule

`DS_NONCONV_POOLS` = `skip` (default, unchanged) | `joint`. It applies only
under `DS_NONCONV_FALLBACK=lm`; otherwise it is reported as `skip`.

Under `joint`, a non-converged component that holds pool nodes is re-solved
by Levenberg-Marquardt over all its free nodes at once:

- A node a pool writes (state, intermediate, step copy, carrier) contributes
  the row G(x) − x, where G is the pool update (`update_pools!`). Any other
  node contributes F(x) − x, as in specs/066. No prior.
- That is the residual the final convergence check already reads, so a
  solution with residual below tolerance is converged in the usual sense.
- F's rows keep their exact per-reaction gradients. G's columns are forward
  differences over the inputs the pool update reads (`pool_update_inputs`:
  the pools' nodes and carriers, the inputs of their step copies, re-evaluated
  drives and own updates, and the supply producers).
- Warm start from the last iterate; `DS_NONCONV_LM_ITERS` (default 400). The
  result is kept only if its residual is lower; otherwise the iterate is
  restored.
- Converged components are never touched.
- The solve reports `nonconv_pools`, `nonconv_pool_components`,
  `nonconv_pool_improved` and `nonconv_pool_rescued`.
- Tests: `test/test_cycles.jl` "a non-converged pool component is re-solved
  jointly". A pool with a loop from its readout back into one of its steps,
  stopped after 3 sweeps, recovers the 500-sweep solution to 1e-7 under
  `joint` (residual ≤ 1.4e-11), on six fixtures; `skip` leaves it
  non-converged.

## Arms (solver only; no rebuild)

- `plj`: canonical, `--server DS_NONCONV_POOLS=joint`, against the specs/074
  arm on canonical.
- `pljd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction --server DS_NONCONV_POOLS=joint`, against
  the same overrides without `joint` on the specs/074 solver (`crsd74`, run
  with these arms), against `crsd` and against canonical.
- `pljd2k`: as `pljd` with `--server DS_MAX_ITERS=2000`, to test whether the
  cap still decides any class (specs/066 amendment 2).

All three arms run on a solver with the specs/074 fix.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`, on paired
cases and on all cases valid in both runs.

**Gate (registered).** `joint` is adopted as the default if `plj` against
canonical gives curator held-out net ≥ 0, experimental net ≥ −5, and no
pathway below −10 (RAF set aside). `pljd` and `pljd2k` inform the
faithful-network work and are not gated.

**Predictions:**
1. `plj`: no prediction changes; the log reports 0 pool components
   re-solved.
2. `pljd`: only Transcriptional Regulation by TP53 moves, on both axes.
3. `pljd`: of the component-solves `crsd` skipped (12 curator, 12
   experimental), at least 9 on each axis are improved, and TP53's
   non-converged curator cases fall from 663 to at most 330.
4. `pljd` vs `pljd2k`: no TP53 case that is converged in both changes class.
5. `pljd`: of the 46 MDM2 overexpression cases that `psrd` read UP against
   curator DOWN, at most 15 read DOWN (`crsd`: 6). Revised from the probe
   above, which found one of them UP (1.76×) at the fixed point: on the
   faithful network the PLK3 loop holds p-S15,S20-TP53 up under MDM2
   overexpression.
6. Run time is reported. No prediction.
