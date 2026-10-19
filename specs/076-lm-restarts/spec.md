# specs/076: retry the non-converged fallback from later iterates

**Status:** pre-registered 2026-10-19, before the arms; results in `research.md`. Canonical: no class changes (gate passes, trivially), converged curator cases +167. Faithful build: every TP53 case converges (306 -> 0) and TP53 curator net −82, the fixed point reading AKT1/AKT2/MDM4 knockout NORMAL as canonical does. Predictions 1-3 hold. Adoption is Adam's call.

Adam, 2026-10-18: "ok" to this, after the specs/075 report.

## The defect (probed)

Faithful build `20261018-0815_ddbd828_dsub`, Transcriptional Regulation by
TP53, solver `c271d1d` (specs/073 and 075 defaults) with
`DS_POOL_SUPPLY_REACH=material` and `DS_CYCLE_CARRIER_RELEASE=reaction`.
Under specs/075 four perturbations stopped converging: CHEK2 knockout, DAXX
overexpression, MDM4 knockout and overexpression (204 curator cases). Probed
with scratch `capsweep.jl` at several iteration caps (readout 5628829):

| perturbation | cap 500 | 750 | 1,000 | 1,500 | 2,000 |
|---|---|---|---|---|---|
| MDM4 overexpression | 4.28×, no | 0.61×, yes | 0.61×, yes | 0.61×, yes | 0.61×, yes |
| DAXX overexpression | 3.47×, no | 0.54×, yes | 0.54×, yes | 0.54×, yes | 0.54×, yes |
| CHEK2 knockout | 3.05×, no | | 1.40×, yes | | 0.11×, no |
| MDM4 knockout | 15.8×, no | | 1.08×, yes | | 15.2×, no |

The specs/073 joint LM starts from the last iterate. From the cap-500 iterate
it stalls; from later iterates it reaches a fixed point, the same one at
every cap where it converges. specs/073's `pljd2k` arm showed the same start
dependence (fewer converged at cap 2,000 than at 500).

## The rule

`DS_NONCONV_RESTARTS` = N (default 0, unchanged), with
`DS_NONCONV_RESTART_SWEEPS` (default 250).

- When the specs/066 fallback (plain or specs/073 joint) leaves a component
  above tolerance, the component is restored to its last iterate (not the LM
  state), iterated `DS_NONCONV_RESTART_SWEEPS` more sweeps, and the fallback
  is run again. This repeats up to N times, stopping at the first attempt
  below tolerance. The state with the lowest residual over all attempts is
  kept.
- Components the fallback skips, bounded ones (`DS_SCC_NEG_MODE`) and
  converged ones are untouched. N = 0 is bit-identical to the current solver.
- The solve reports `nonconv_restarts_max`, `nonconv_restarts` (used) and
  `nonconv_restart_rescued`. The specs/066 and 073 counters now describe the
  first attempt only.
- Tests: `test/test_limit_cycles.jl` "restarts". The specs/066 Goodwin ring
  with an 8-iteration LM budget: from the cap-500 iterate LM stalls; one
  restart (the cap-750 start) converges to the analytic fixed point
  D = 80^(1/10). Restarts stop at the first rescue. When none rescues, all
  are used and the best state is kept.

## Arms

N = 4 (sweeps 250, so starts at 500, 750, 1,000, 1,250 and 1,500), chosen from
the probe above before any arm, and not tuned afterwards.

- `rst`: canonical, `--server DS_NONCONV_RESTARTS=4`, against
  `results/3177881/rvc` (identical to canonical defaults).
- `rstd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction --server DS_NONCONV_RESTARTS=4`, against
  `rvcd` (`results/3177881/rvcd`).
- `rstd2k`: as `rstd` with `--server DS_MAX_ITERS=2000`.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade` (its
"baseline" is the `--cases` file), paired and on all cases valid in both, and
checked case by case.

**Gate (registered).** `DS_NONCONV_RESTARTS=4` is adopted as the default if
`rst` against canonical gives curator held-out net ≥ 0, experimental net
≥ −5, and no pathway below −10 (RAF set aside). `rstd` and `rstd2k` inform
the faithful-network work and are not gated.

**Predictions:**
1. `rst`: every class change is in a solve that `rvc` reports non-converged,
   and converged curator solves do not decrease.
2. `rstd`: at least 3 of the 4 perturbations above converge, and TP53
   non-converged curator cases are at most 102 (`pljd`; `rvcd` 306).
3. `rstd` vs `rstd2k`: no TP53 case converged in both changes class.
4. Run time per axis is reported; no prediction.
