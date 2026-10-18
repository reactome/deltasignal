# specs/074: the hill_sat AND is wrong just below its 100× cap

**Status:** pre-registered 2026-10-18, before the arms.

Found by specs/073's probe (Adam's "ok" to #2, 2026-10-18): Levenberg-Marquardt
stalled on the TP53 pool component because one AND's gradient was wrong. The
AND was right; its output was not.

## The defect

`DS_AND_MODE=hill_sat` (default) caps the fold product at 100× with a smooth
minimum. With d = raw − 1 on the internal scale and ε = `DS_HILL_SAT_EPS`
(default 1e-9 since `7b57ae3`, 2026-09-17), the code computes

    denom = d + sqrt(d² + ε²)
    hi = denom > 0 ? 1 − ε² / (2·denom) : (raw + 1 − sqrt(d² + ε²)) / 2

The first form was introduced to stop a cancellation for raw ≫ 1 (wide
reactions; the record is in `test_and_curves.jl`). But `denom > 0` holds
almost everywhere, so raw below the cap also takes it. There d < 0 and
`d + sqrt(d² + ε²)` cancels:

- For raw below about 0.875 (87.5×), d² + ε² rounds to d², denom is exactly
  0, and the second form runs. No error.
- For raw between 0.875 and 0.998, denom is a few ulps of noise, so `hi` is
  off by up to 0.089 (8.9 on the UI scale), quantized and not monotone.
  Measured on 200,001 points (scratch `hs.jl`): 12,197 points off by more
  than 1e-6, worst at raw 0.875.
- At ε = 1e-3 (the default before `7b57ae3`) the worst error was 2e-10.

Traced instance (specs/073): an AND of p-S15,S20-TP53-side node 3451 (1.21×)
and pinned MDM2 (80×) should read 96.88×. It reads 96.40×, 97.12× or 97.60×
depending on the last digits of its input. Its finite-difference gradient is
0 where the exact one is 70.8.

A value in that band is UP either way, so the error alone rarely changes a
class. What it breaks is any loop that carries such a node: the map has no
fixed point to tolerance there, so the iteration cannot converge and the
specs/066 fallback works on a wrong Jacobian.

## The fix

Branch on the sign of d. Above the cap keep `1 − ε² / (2(d + √(d² + ε²)))`.
Below it use the mirror form `raw − ε² / (2(√(d² + ε²) − d))`, which is
algebraically identical to `(raw + 1 − √(d² + ε²)) / 2` and has no
cancellation. Worst error after the fix: 1.1e-16 on the same 200,001 points.

No flag: this is a numerical defect in the default model, not a modelling
choice. `test/test_and_curves.jl` "exact and monotone just below the ceiling"
fails 16 of its assertions on the old code and passes on the fix.

Out of scope: the inhibition-side `hill_sat` (`DS_INHIBITION_MODE=hill_sat`,
not the default) uses the textbook form, which is exact below its cap and
loses precision only for h_raw beyond about 1e15.

## Exposure, before the arms

Canonical `20261017-1246_5e663d9`, `results/6fb4329`: 110 of 23,511 valid
curator readouts and 38 of 845 experimental readouts have a predicted value
in the band. Non-converged: 581 curator, 41 experimental. Interior nodes are
not counted.

## Arms

- `hsc`: canonical, solver with the fix, no overrides, against
  `results/6fb4329`.
- `hscd`: faithful build `20261018-0815_ddbd828_dsub`, no overrides, against
  `results/40adace/dsub`. Informs specs/070–073; not gated.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`, paired and
on all cases valid in both.

**Gate (registered).** A defect fix. Kept if `hsc` gives curator held-out net
≥ 0, experimental net ≥ −5, and no pathway below −10 (RAF set aside). If it
fails, the failing cases are traced before anything is decided, because a
loss here means some prediction relied on the noise.

**Predictions (`hsc` against `6fb4329`):**
1. At most 1% of valid curator cases (235) change class.
2. Converged curator solves do not decrease; non-converged valid cases do not
   increase (581).
3. Cases that change class are concentrated in solves that were
   non-converged in `6fb4329`: at least half of them.
4. Readouts in the band keep their class (UP), except where a case also
   changes convergence.
