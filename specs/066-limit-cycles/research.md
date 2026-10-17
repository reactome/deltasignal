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
  Response). The 454-node component improves but stays above tolerance;
  their class is unchanged (16 right, 2 wrong at the default cap) and they
  are credited to neither arm. *Corrected by amendment 2: under `lm` their
  class no longer depends on the cap.*

## Decision (by the registered rule)

The gate passes: held-out +1 (≥ −15), experimental 0 (≥ −15), both read
with RAF set aside, and converged solves rise on both axes. By the rule,
`DS_NONCONV_FALLBACK=lm` is adopted as a faithfulness fix. Making it the
default was Adam's call: **adopted as the default 2026-10-16**.

## Amendment 2 (2026-10-16, after adoption): is anything still cap-decided?

Above tolerance is not the same as decided by the cap. Under `off`, the 23
cases were cap-decided because their class changed between caps (week
review section 5). Re-tested under the adopted default, on main
(`0709a1e`, build `20261016-0006_da08370`): every case whose solve did not
converge in `results/0709a1e` (597 curator, 41 experimental) was re-solved at
`DS_MAX_ITERS` 500, 1,000, 2,000 and 4,000 (probe
`~/deltasignal-catalogs/probe_osc/nccap.jl`, output `nccap.tsv`).

- Cap 500 reproduces `results/0709a1e` exactly: 0 class and 0 value
  mismatches.
- **0 of 638 cases change class across the four caps.**
- 600 values move by less than 0.1%. Of the 38 that move more (Mitotic G1 14,
  DSB Repair 9, Caspase via death receptors 6, VEGF 5, DSB Response 3,
  Fanconi 1), the one nearest a class boundary is FANCD2 at 6785124
  (0.38×, 0.12×, 0.023×, 0.38×; DOWN at every cap). Next nearest is RNF168
  at 5693527 in DSB Response (56×, 5.6×, 45×, 45×; UP at every cap).
- An LM budget of 2,000 instead of 400 changes no RNF8 class either
  (`dsbcap.tsv`).

So the 18 RNF8/RNF168 held-out cases, and every other non-converged case,
now have a class the cap does not choose. Their values are not fixed
points (the residual stays above tolerance), so magnitudes in those solves
are approximate. The open item is narrower than recorded at adoption:
nothing is cap-decided, and 43 curator and 10 experimental solves are still
inexact in magnitude.
