# specs/163 research: the measurement (2026-11-16)

Scorer `bench/analysis/163/score_163.py` (fixture-checked, committed before
the read). Recorded predictions only; own markers excluded.

## korkut (98 record; 78 drug sets)

| | precision@10 | precision@20 |
|---|---|---|
| model ranking (|predicted log2 fold|) | **0.105** | 0.060 |
| reachability ranking | 0.086 | 0.072 |
| random | 0.027 | 0.027 |

Gate (precision@10, model vs reachability, 95% bootstrap over drug sets):
**+0.019 [+0.007, +0.032], PASS.**

## CPPA (specs/121 final frozen run; numbers private under CPPA terms)

**Gate: FAIL.** The model's ranking does not beat reachability at
precision@10 (the interval includes zero), and it trails reachability at
precision@20.

## The reading

1. **The usable ranking signal is mostly reachability.** On both datasets,
   ranking by "does a route reach this readout" puts real, correctly
   directed changes several-fold above random. That is the curated
   structure working.
2. **Magnitude ordering on top of reachability barely helps, and does not
   replicate.** On korkut the model's top 10 adds a little over
   reachability, but by the top 20 it falls below; on CPPA it adds nothing.
3. **A correction to the motivating claim.** The AUROC "consistent ranking
   signal" (korkut 0.65–0.67; the CPPA figure, removed 2026-11-20 under FR-003, is in the private notes) that motivated this spec scores unreachable readouts at
   the bottom (predicted fold 1), so it largely measures reachability too.
   It was not evidence that the magnitudes rank well.

## Serving

The registration allows a ranked view where the gate passes (korkut only).
Given the CPPA failure and the reversal at precision@20, **no magnitude-ranked
view is recommended**. The product already says what the record supports: a
reachable readout is several-fold more likely to change than an unreachable
one ("not reached here", specs/146).
