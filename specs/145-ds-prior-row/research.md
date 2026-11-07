# specs/145 research: the DS response row for the unseen drug, measured

## Verdict

**H1 fails as registered; so do H1′ and H2.** But unlike specs/144 this one
leaves a real, positive trace that a larger test could confirm.

| arm | LODO mean per-drug corr | single-drug corr | balanced acc (log2 1.5) |
|---|---|---|---|
| lr (control, = specs/118 rebuild) | 0.7697 | 0.000 | 0.5219 |
| lr+ds_fold (the 118 feature arm) | 0.7686 | −0.020 | 0.5216 |
| **lr+dsrow** | **0.7794** | **0.167** | **0.5387** |
| lr+reachrow (H1′ control) | 0.7717 | 0.126 | 0.5281 |
| lr+dialrow (H2) | 0.7660 | 0.167 | 0.5234 |

- H1 (dsrow − lr): mean +0.0097, bootstrap **[−0.0134, +0.0369]** — includes 0
  at n = 12 drugs. Fails as registered.
- H1′ (dsrow − reachrow): **[−0.0015, +0.0190]** — misses by a hair; about
  two-thirds of the gain survives the reachability control (contrast
  specs/113, where the control ate essentially everything).
- H2 (dial on top, markers excluded): **[−0.0026, +0.0023]** — nothing, again
  consistent with specs/144: shrinkage does not help correlation.

## Traces

- **The negative-c scare resolved**: the fitted scale c ≈ −0.04 in every fold
  because `pert.csv` encodes drug presence negatively (the CellBox
  convention; every nonzero D entry is in [−3.48, −0.27]), flipping the
  fitted R rows; the injection multiplies back through D, so the signs
  cancel. DS's rows themselves correlate **positively** with the measured
  single-drug responses for 8 of 12 drugs (ZS +0.60, AK +0.57, PLX +0.43,
  901 +0.34, SR +0.30, HN/NT/P6 ≈ 0; Tm −0.28 the one real inversion),
  mean +0.17 — reproducing specs/118's +0.157 independently.
- **Where the gain lives**: exactly where the design says it must — the LODO
  model predicts literally zero for the held-out drug's direct effect
  (single-drug corr 0.000); the DS row supplies 0.167 of correlation there
  and +0.017 balanced accuracy overall. The specs/118-style *feature* arm
  (lr+ds_fold) adds nothing; injecting the row *through the design matrix*
  is the right construction.

## What this redirects

1. n = 12 drugs is the wall: a +0.01-mean effect cannot clear a drug-level
   bootstrap at this size. The construction transfers unchanged to
   **HPN-DREAM** (already-scored FIT data, 30 contexts × 48 readouts, a
   leaderboard to place against) — the natural next measurement, same
   registration pattern.
2. Tm (temsirolimus) is the one anti-correlated row — the known
   mTORC1→IRS1 feedback miss (specs/118 traced it). A single structural fix
   there is worth more to this construction than any coefficient.
3. The dial composes with nothing (144, 145-H2): engagement conditioning is
   a false-alarm tool for the product mode, not an accuracy lever. Stop
   pulling it.

## Provenance

Harness: `bench/analysis/lr_118.py` imported whole (lodo split, fit, boot;
control reproduces its 0.770). DS rows: frozen record
`20261026-1349_8b6527e / a490cdc / k134_abund_a3`, single-drug clipped log2
folds, the harness's own `ds_fold` matrix. α: the specs/144 pooled estimator
imported from `bench/analysis/144/engagement_dial.py`. Script:
`bench/analysis/145/hybrid_145.py`; korkut2015 is FIT
(`dataset_roles.require_role` enforced). Registered at 7ae5c75 before
measurement.
