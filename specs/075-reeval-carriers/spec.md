# specs/075: a carrier re-evaluated in another pool's drives reads by the carrier rule

**Status:** pre-registered 2026-10-18; results in `research.md`. Canonical identical (gate passes, trivially). Faithful build: TP53 curator +252, experimental +15, MDM2 overexpression 46/46 DOWN; predictions 1-2 hold, 3 fails (TP53 non-converged 102 -> 306). Adoption is Adam's call.

Adam, 2026-10-18: "yes" to tracing why the faithful build reads MDM2
overexpression UP after specs/073, and "ok" to this fix.

## The defect (traced)

Faithful build `20261018-0815_ddbd828_dsub`, Transcriptional Regulation by
TP53, MDM2 overexpression (UI 80), solver `6c2269b` with
`DS_POOL_SUPPLY_REACH=material` and `DS_CYCLE_CARRIER_RELEASE=reaction`.
The solve converges (specs/073). Dumps and probes are in
`~/deltasignal-catalogs/trace_073/`.

1. The readouts (5628829, 69685, 6799246) read TP53 tetramer
   (R-HSA-3209194, `52dc6033`) at 1.76×. Canonical reads them at 0.60×.
2. That node is a state of pool2: TP53 tetramer → MDM2:TP53 → PolyUb-TP53
   (base) → USP7 deubiquitination (R-HSA-3215310) → TP53 tetramer. At the
   solution: supply s = 0.352, π0 = (0.9, 0.1), π = (0.50, 0.50), drives
   USP7 step 50.5 and MDM2 binding 4.96 and 6.28. TP53 tetramer =
   0.352 × 0.499 / 0.1 = 1.76.
3. The USP7 step's inputs are water, PolyUb-TP53 (held at baseline in the
   drive) and USP7. Inside the drive USP7 reads **50.5×**, while its own value
   is 1×.
4. USP7 (R-HSA-3215278) is the carrier of pool4 (the MDM2 cycle). The carrier
   rule writes it from its producers other than the release steps: none left,
   so 1×. But USP7 is also a pool-fed input of pool2's steps, so pool2's drives
   re-evaluate it from its whole reaction (specs/039 amendment 5). That
   reaction is the OR of its two release copies (R-HSA-3215295). One reads a
   DAXX complex as input and catalyst, and saturates at 100× under MDM2
   overexpression. The mean is 50.5.

So the release flux that specs/049 amendment 7 and specs/072 keep out of a
carrier comes back through the drive re-evaluation. Canonical solves no pool
in TP53, so the case does not arise there.

## The rule

`DS_CYCLE_REEVAL_CARRIERS` = `full` (default, unchanged) | `carrier`.

- Under `carrier`, when a pool's drives re-evaluate a pool-fed input that is
  the carrier of a pool in the same component, it reads by the carrier rule:
  the mean fold of its producers other than the release steps of every pool
  listing it, times its modifier (`carrier_value`, shared with the carrier
  update). Everything else is re-evaluated as before.
- It applies only while the carrier rule is on. The solve reports
  `cycle_reeval_carrier_rule` and `cycle_reeval_carriers` (pool-carrier
  pairs, counted under both modes).
- Tests: `test/test_cycles.jl` "a carrier re-evaluated in another pool's
  drives reads by the carrier rule". Two pools in one component, the second's
  return step catalysed by the first's carrier K. Under `full` its state B
  moves with the first pool's supply (B 0.025× at G 80×) and drags the first
  pool down (Ss 2.2× instead of 80×). Under `carrier` A and B stay at 1× and
  Ss follows G exactly.

## Structural count, before any arm

`count_reeval_carriers.jl` (this directory): solve every pathway at rest and sum
`cycle_reeval_carriers`:

- canonical `20261017-1246_5e663d9`: **0** pairs in every pathway.
- faithful `20261018-0815_ddbd828_dsub` (with the specs/071 and 072 rules):
  **15** pairs, all in R-HSA-3700989 (Transcriptional Regulation by TP53, 5
  pools).

## Arms (solver only; no rebuild)

- `rvc`: canonical, `--server DS_CYCLE_REEVAL_CARRIERS=carrier`, against
  `results/501e480/hsc` (same solver behaviour as `6c2269b` on canonical:
  specs/073 moves nothing there).
- `rvcd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction --server DS_CYCLE_REEVAL_CARRIERS=carrier`,
  against `pljd` (`results/ae9c56f/pljd`, the same overrides under the
  specs/073 default).

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade` (its
"baseline" column is the `--cases` file), paired and on all cases valid in
both, and checked case by case.

**Gate (registered).** A defect fix. `carrier` is adopted as the default if
`rvc` against canonical gives curator held-out net ≥ 0, experimental net
≥ −5, and no pathway below −10 (RAF set aside). `rvcd` informs the
faithful-network work and is not gated.

**Predictions:**
1. `rvc`: identical to canonical on every case (the count is 0).
2. `rvcd`: only TP53 moves. In TP53, the 46 MDM2 overexpression cases that `pljd` reads UP
   against curator DOWN: at least 30 read DOWN or NORMAL. USP7 at 1× in
   the drive removes the 50× deubiquitination, so pool2's split moves back
   toward PolyUb-TP53. The PLK3 loop may still hold some readouts up.
3. `rvcd`: TP53 non-converged curator cases do not rise above `pljd`'s 102
   by more than 50.
