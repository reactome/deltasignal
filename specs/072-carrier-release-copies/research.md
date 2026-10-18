# specs/072 research: results of arms crs and crsd

Solver `d1d8bb9`, run 2026-10-18 with `scripts/run_arm.sh`.

| arm | build | overrides | results |
|---|---|---|---|
| crs | `20261017-1246_5e663d9` (canonical) | `DS_CYCLE_CARRIER_RELEASE=reaction` | `results/d1d8bb9/crs` |
| crsd | `20261018-0815_ddbd828_dsub` | `+ DS_POOL_SUPPLY_REACH=material` | `results/d1d8bb9/crsd` |

## The counter over-reported (amendment 1, diagnostics only)

The logs report 2,016 copies added on canonical and 2,088 on the faithful
build. The structural count was 0 and 2.

- The difference is a carrier shared by several pools. In RAF, one GAP is a
  carrier of four pools, and each pool lists a different node of the same
  release reaction (R-HSA-5658231).
- `index_pools` runs per pool, so it "added" the other pools' listed nodes.
  The solver already leaves those out (specs/049 amendment 7).
- So no behaviour changed on canonical, but the count did.

Amendment 1 counts only copies that no pool lists for the carrier.
Re-checked on the bundles: RAF (canonical) now reports 0, TP53 (faithful) 2.
The new count matches the structural count. The test is the shared-carrier
testset in `test_cycles.jl`: under `reaction` it adds 0 and every activity is
unchanged. Predictions are untouched, so the arms were not re-run.

## crs: canonical

Every case is identical to `results/6fb4329`, curator and experimental, down
to `pred_ui`.

- Prediction 1: **holds** for predictions. The log figure was the counter
  artefact above.
- The gate passes, but only because nothing moves. The decision is Adam's.

## crsd: faithful build

`holdout_report.py --set-aside RAF_MAP_kinase_cascade`:

| comparison | held-out | tuning | experimental |
|---|---|---|---|
| crsd vs psrd | +0 | **−127** (84 / 211, TP53 only) | −3 (12 / 15, p 0.70, TP53 only) |
| crsd vs canonical | −171 | −302 | −20 (p 0.007) |

- Prediction 2 (only TP53 moves): **holds**.
- Prediction 3 (≥ 40 of the 46 MDM2 overexpression cases read DOWN):
  **fails**. 6 read DOWN, 4 read NORMAL, 36 still read UP.
- Prediction 4 (≥ 45 of the 49 MDM2 knockouts stay right): **holds**, 49.

TP53 cases, crsd vs psrd:

| | fixed | broken |
|---|---|---|
| BRCA1 knockout | 39 | 0 |
| ATM knockout | 36 | 0 |
| MDM2 overexpression | 7 | 0 |
| MDM4 knockout / overexpression | 0 | 50 / 39 |
| AKT1, AKT2 knockout | 0 | 100 |
| AKT1, AKT2 overexpression | 0 | 12 |
| CHEK2 knockout | 1 | 10 |
| CDKN2A knockout | 1 | 0 |

### Trace: MDM2 overexpression (cap 500, `trace_071/xbase`, `xmoe`)

The rule does what it should at USP7:

- Free USP7 reads 1× (100× in psrd), and so does DAXX.
- Free TP53 tetramer falls to 0.33× and PolyUb-TP53 to 0.09×: MDM2 now lowers
  TP53.

The readouts still read UP because of a positive loop that never converges.
p-S15,S20-TP53 tetramer (R-HSA-3222171) has three producers:

- CHEK1 (R-HSA-6799246) and CHEK2 (R-HSA-69685) phosphorylation read 0.395×.
- PLK3 phosphorylation (R-HSA-6805285) reads 2.32×. Its input, p-S15-TP53,
  reads 0.35×, but its catalyst PLK3 reads 4.25×.
- PLK3 is induced by p-S15,S20-TP53 itself (R-HSA-6802165 → R-HSA-6802170).

The OR mean is 1.29×, and 26 of the 36 UP readouts read exactly 1.29×. The
component does not converge (`nonconv_skipped_pools=1`).

### Most of the movement is non-converged

| | psrd non-converged | crsd non-converged | crsd converged |
|---|---|---|---|
| broken | 211 | 151 | 60 |
| fixed | 84 | 82 | 2 |

Every case that moved was non-converged in psrd, and 233 of 295 are
non-converged in crsd. TP53 non-convergence is 663 (it was 765 in psrd).
Exceptions: the 50 MDM4 knockouts converge in crsd and are wrong; not traced.

## Reading

- The rule fixes the defect it was registered for: no carrier is fed by an
  unlisted copy of its own release.
- On canonical it moves nothing. On the faithful build it touches two carriers
  in TP53, and that pathway's classes are decided in a non-converged component.
- TP53 tuning numbers on the faithful build therefore cannot judge a rule until
  the component converges. That makes the next step the LM fallback for
  components that hold pools (specs/066 skips them), before any further TP53
  rule.
