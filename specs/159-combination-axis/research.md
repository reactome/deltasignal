# specs/159 research: the measurement and its verdict

Run 2026-11-16, entirely offline on the recorded 98 korkut TSV
(`results/cff961d/r98fix_korkut`), scorer `bench/analysis/159/combo_159.py`
(fixture-checked on computed known answers before the real read; role guard
korkut2015 FIT; registered nulls, seed 159).

## The verdict — GATE FAILS, and the anatomy is the finding

1,442 joinable pair cases over 50 pairs (8,017 pair rows had no scored single
counterpart for the same (pathway, readout) and are counted, not used):

| predictor | accuracy | balanced |
|---|---|---|
| recorded pair solve | 0.7809 | 0.4892 |
| product null (fold_a × fold_b) | 0.7836 | 0.4903 |
| extreme null (larger \|log2\|) | 0.7836 | 0.4903 |

Gate: pair − product = **−0.0028 [−0.0087, +0.0000]** over 50 pairs — FAIL.

**The anatomy**: the pair solve agrees with the product null on **99.7%** of
cases — the engine's simultaneous-pin combination is, in practice,
multiplicative composition of the single effects. On the 5 cases where the
solve DOES model an interaction beyond the product, it is right **0** times
and the null 4. On the 41 combination-specific measured cases (the pair's
measured call differs from BOTH singles' — the class only a genuine
interaction model can win), all three predictors score an identical 11/41.

## What this means, honestly

- Simultaneous pins are a CONVENIENCE, not a modelled interaction: a user
  gets the same answer (to 99.7%) by multiplying single-drug predictions.
  That is coherent with the solver's design (specs/010: AND multiplies
  fold-changes; the propagation is built from products), and it is now
  measured rather than assumed.
- Per FR-004 the serving decision is: **nothing is served** — no scope line,
  no combination-validity claim. The product keeps accepting simultaneous
  pins exactly as before, with no new accuracy implied.
- The 41-case combination-specific class is the honest target any future
  interaction modelling must beat (currently 11/41 for everything); it is on
  the record for that purpose.

This closes improvement #256 and with it the whole approved backlog
(items 1–6): four served features (switches, intervals, soft observations,
fragility), one served endpoint (inverse queries), and two measured negatives
recorded at full strength (influence pre-ranking, combination interactions).
