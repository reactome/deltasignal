# specs/165 research: the measurement (2026-11-16)

Scorer `bench/analysis/165/score_165.py` (fixture-checked, committed before
the read). Change calls on reachable cases; labels truth-blind on every
scored axis.

| axis | F0 all: precision (n) | F3 robust + confident: precision, coverage | Δ [95% boot] | gate |
|---|---|---|---|---|
| curator held-out | 0.888 (6,021) | **0.920**, 77% | +0.032 [+0.017, +0.049] | PASS |
| experimental (MP-BioPath) | 0.876 (687) | **0.911**, 79% | +0.035 [+0.009, +0.070] | PASS |
| korkut (drug RPPA) | 0.109 (1,097) | **0.152**, 58% | +0.043 [+0.026, +0.062] | PASS |

Each label helps on its own (F1 robust, F2 confident), and the two combine:
the combination is the best on every axis.

## The reading

1. **The served labels work as a filter on every axis.** Keeping only robust,
   confident calls raises precision with a modest loss of coverage. On korkut
   that is a 39% relative gain (0.109 → 0.152) while keeping 58% of the calls.
2. **The low precision is specific to drug readouts.** On gene perturbations
   (the curator expectations and MP-BioPath's measured knockouts and
   over-expressions) about 88% of the model's change calls are right, rising
   to about 91–92% under the filter. The drug RPPA axis is where precision is
   low (0.11 → 0.15). This agrees with specs/164: drug responses behave like
   weaker, saturating propagation than the gene-level axes the model was
   built on.

## Serving

Per FR-003, a pass documents the labels as a usable filter, with these
measured figures (docs/API.md). No engine change.

## Note, 2026-11-20 (specs/167)

The gates and figures above stand. Reading 2 does not: specs/167 shows that on every axis precision is the share
of reached readouts that changed (88–94% on the gene axes, 14% on korkut) plus at most 0.024, and that the recorded
axes confound the perturbation type with how readouts were chosen. See `specs/167-reached-base-rate/research.md`.
