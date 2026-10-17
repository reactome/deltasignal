# specs/066 research: arm `nclm`

Arm `nclm`: `run_arm.sh nclm --server DS_NONCONV_FALLBACK=lm`, solver
`c04371e`, canonical build `20261016-0006_da08370`, against
`results/0406a4b` (same `src/` and `bench/` apart from this rule, default
off). Results: `builds/20261016-0006_da08370/results/c04371e/nclm/`.

## Convergence

| axis | converged solves, `0406a4b` → `nclm` | components the fallback ran on | improved | rescued (below tolerance) | skipped (pools) |
|---|---|---|---|---|---|
| curator | 1,654 → 1,682 of 1,725 | 64 | 64 | 33 | 0 |
| experimental | 222 → 234 of 244 | 17 | 17 | 12 | 0 |

Per case: curator 680 cases now come from a converged solve (0 the other
way), leaving 597; experimental 108, leaving 41. TP53 falls from 612
non-converged curator cases to 51, and from 108 experimental to 6. What is
still above tolerance: DSB Repair 162 (the RNF8/RNF168 component the probe
could not rescue), WNT 111, Mitotic G1 100, TP53 51, TGF-beta 42.

## Predictions

| axis | split | net | fixed | broken | p |
|---|---|---|---|---|---|
| curator | held-out | +1 | 2 | 1 | 1.0 |
| curator | tuning | 0 | 0 | 0 | — |
| experimental | all (with or without RAF) | 0 | 0 | 0 | — |

Held-out accuracy 88.87% in both; macro-F1 0.8559 → 0.8561 (corrected
2026-10-16: first written reversed).

Canonical re-run on main (solver `0709a1e`, `results/0709a1e`) reproduces
`nclm` on every prediction and value, both axes. RAF has no
moved case, so the set-aside reading is the same.

The three class changes are exactly the three the spec predicted:
- FANCA and FANCC overexpression, readout 6785124: 0.0066 → 2.00 (DOWN →
  UP, curator UP). The Fanconi Goodwin ring is rescued; the value is the
  analytic fixed point and no longer depends on the cap.
- SH3GL1 knockout, readout 167679: 1.31 → 0.0048 (UP → DOWN, curator UP).
  The EGFR component is rescued; at its fixed point every SH3GL1 readout
  reads DOWN. The previous UP was a phase of the limit cycle (cap 500: UP;
  1000: UP; 2000 and 4000 differ, week review section 5), not a stable
  answer. The ring runs through CBL ubiquitination edges that deplete free
  EGF-family ligands (the specs/057 generator class), so the fixed point
  inherits that defect.

No other prediction changed, though 552 curator and 59 experimental values
moved (all in solves the fallback touched).

## The 23 cap-decided held-out cases

- **5 now converge** and are decided by the fixed point: FANCA and FANCC at
  6785124 (DOWN → UP, right); SH3GL1 at 167679 (UP → DOWN, wrong); SH3GL1
  and ARHGEF7 at 179882 (DOWN at cap 500 and at the fixed point, wrong both
  ways).
- **18 still do not converge**: RNF8 ×15 and RNF168 ×3 (DSB Repair and DSB
  Response). The 454-node component improves but stays above tolerance, so
  these are still the cap's choice; their class is unchanged (16 right,
  2 wrong at the default cap) and they are credited to neither arm.

## Decision (by the registered rule)

The gate passes: held-out +1 (≥ −15), experimental 0 (≥ −15), both read
with RAF set aside, and converged solves rise on both axes. By the rule,
`DS_NONCONV_FALLBACK=lm` is adopted as a faithfulness fix. Making it the
default was Adam's call: **adopted as the default 2026-10-16**.
