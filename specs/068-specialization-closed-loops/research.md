# specs/068 research: arm `sctwin`

Generator: LNG branch `specialization-closed` at `b251d86` (LNG PR, default
off). Build `20261017-0916_b251d86_sctwin` (`LNG_SPECIALIZATION_CLOSED=twin`),
92 of 92 pathways, 500,965 edges (canonical `20261017-0654_d3ebf38`: 500,965).
Arm: `run_arm.sh sctwin --catalog …_sctwin`, solver `ca5c3b4` (`src/` and
`bench/` identical to main `a177c6f`), against canonical `results/fa2c451`.
Results: `builds/20261017-0916_b251d86_sctwin/results/ca5c3b4/sctwin/`. Case
comparison: `~/deltasignal-catalogs/probe_spec/cmp.py <arm> <baseline>`.

## Content

- `logic_network.csv`: the same edge count in every pathway, 500,965 in all.
- `carrier_loops.csv`: against canonical, exactly the 10 PDGF free forms
  (copies of R-HSA-186811) differ, all `closed` True; against the specs/066
  build (`20261016-0006_da08370`, before the bridges) no row differs.

## Results

| axis | net | fixed | broken | McNemar p |
|---|---|---|---|---|
| curator held-out | **+18** | 18 | 0 | 7.6e-6 |
| curator tuning | 0 | 0 | 0 | |
| experimental | 0 | 0 | 0 | |

All 18 are PLCG1 in Signaling by PDGF, knockout and overexpression, 9
readouts each: exactly the 18 specs/067 broke. No other case moves, so RAF
has none and both readings are the same.

Curator held-out in release: 89.81% / 0.8695 → **89.91% / 0.8706**. All
pathways, every case: 86.95% → 87.02% (20,972 / 24,100), macro-F1 0.8403 →
0.8411; MP-BioPath 83.34%, so +3.68pp. Experimental 606/849 (71.38%),
unchanged.

## Predictions, scored

| prediction | outcome |
|---|---|
| `logic_network.csv` identical; only the 10 PDGF rows of `carrier_loops.csv` differ, all back to closed | **held** (edge counts identical; content compared by count, not byte) |
| class changes only in PDGF | **held** |
| PLCG1: at least 16 of 18 fixed | **held**: 18 |
| at least 75 of specs/067's 80 PDGF fixes kept | **held**: 80 (nothing broken) |
| held-out net +10 to +20; tuning 0; experimental 0 | **held**: +18, 0, 0 |

## Decision (by the registered rule)

Passes: held-out +18 at p 7.6e-6, experimental 0 (≥ −5), worst pathway 0
(≥ −10). **Adopted as the generator default 2026-10-17 (Adam).**
