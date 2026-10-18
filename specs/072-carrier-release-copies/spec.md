# specs/072: every copy of a carrier's release reaction is a release step

**Status:** pre-registered 2026-10-18; arms not yet run.

Adam, 2026-10-18: "ok" to fixing first the defect that specs/071 found.

## The defect

The carrier rule (`DS_CYCLE_CARRIERS=on`, specs/039 and specs/049 amendment 7)
makes an enzyme's free form read its producers OTHER than the release steps
that `pool_carriers.csv` lists. Those steps close the enzyme's own gain-1
loop. The table names release reaction NODES. When the generator emits a
second variant copy of a release reaction, and that copy reads an
intermediate outside the pool, the copy is not listed. It then feeds the
carrier the flux through that intermediate.

The case comes from specs/071 research. On the faithful build in Transcriptional
Regulation by TP53, USP7 (R-HSA-3215278) is a pool4 carrier released by "USP7
deubiquitinates MDM2" (R-HSA-3215295):

- That reaction has two nodes. `8e405807` is listed. `33469e3e` reads a DAXX
  complex variant (`bb11841f`) that is not a pool4 member.
- Under MDM2 overexpression that variant reads 100×, so free USP7 reads 100×.
  USP7 then deubiquitinates TP53, and the readouts go UP.
- 46 MDM2 overexpression cases read UP in arm `psrd` where canonical and the
  curators read DOWN.

## The rule

`DS_CYCLE_CARRIER_RELEASE` = `listed` (default, unchanged) | `reaction`.

- Under `reaction`, a producer of a carrier's free form whose stable id is the
  stable id of one of its listed release nodes is a release step too. It is a
  copy of the same Reactome reaction.
- It applies only while the carrier rule is on. The solve reports
  `cycle_carrier_release` and `cycle_carrier_release_added`.
- Tests: `test/test_cycles.jl` "every copy of a listed release reaction is a
  release step".

## Structural count, before any arm

`count_release_copies.py` (this directory) lists every pool carrier and its
other positive producers that share a stable id with a listed release:

- canonical `20261017-1246_5e663d9`: 71 carriers, **0** copies to add.
- faithful `20261018-0815_ddbd828_dsub`: 82 carriers, **2** copies in 2
  carriers, both in R-HSA-3700989 (Transcriptional Regulation by TP53).

## Arms (solver only; no rebuild)

- `crs`: canonical build, `--server DS_CYCLE_CARRIER_RELEASE=reaction`,
  against `results/6fb4329`.
- `crsd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction`, against `psrd` (`results/dc7f56a/psrd`)
  and against canonical.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`, on paired
cases and on all cases valid in both runs.

**Gate (registered).** This is a defect fix. `reaction` is adopted as the
default if `crs` against canonical gives curator held-out net ≥ 0,
experimental net ≥ −5, and no pathway below −10 (RAF set aside). Given the
count above, `crs` is expected to be identical to canonical. `crsd` informs the
faithful-network work and is not gated.

**Predictions:**
1. `crs`: no prediction changes, and the log reports 0 copies added.
2. `crsd`: only TP53 moves (curator and experimental). The log reports copies
   added in TP53 solves only.
3. `crsd`: at least 40 of the 46 MDM2 overexpression cases that `psrd` reads UP
   against curator DOWN read DOWN.
4. `crsd`: at least 45 of the 49 MDM2 knockout cases `psrd` got right stay
   right.
5. Non-convergence is reported, not predicted. The CDKN2A knockout classes are
   decided by the iteration cap in `psrd` (specs/071), and this rule does not
   address that.
