# specs/179 research (2026-11-27)

Registration `7d59d85c`. Code `bench/analysis/179/` (`direct_179.py`, `engage_179.py`, tests). Predictions: the
specs/175 file (Release 97 build, abundance rule); the Release 98 predictions were still running. Output outside the
repo (`~/deltasignal-catalogs/analysis/179/engage_r97.json`). Both halves of Lun 2019 are FIT.

## Who is engaged

85 proteins have a measured direct pair. **36 raise a direct substrate (engaged), 49 do not.** Test pairs (direct and
self pairs left out): 1,852, from 85 proteins.

## E1. The model is far better for engaged proteins

| model's change calls | engaged | not engaged | difference |
|---|---|---|---|
| pooled: right direction | **36%** of 196 | **6%** of 201 | +0.30 [+0.21, +0.40] |
| half A | 33% of 144 | 4% of 57 | +0.29 [+0.16, +0.41] |
| half B | 46% of 52 | 7% of 144 | +0.39 [+0.24, +0.53] |
| pooled, with the specs/177 rule | 47% of 253 | 6% of 201 | +0.41 [+0.30, +0.51] |
| pooled, tyrosine-specific kinases left out | 54% of 37 | 6% of 190 | +0.48 [+0.24, +0.72] |

Engaged proteins also change far more of what they have a route to: 44% of pairs with a route against 15% for
not-engaged proteins (32% and 13% without a route).

## E2. It is partly general activity, and partly the first hop

The placebo (the protein raises a readout it has no route to) separates precision too: +0.25 [+0.15, +0.35] pooled.
But inside each placebo group the direct substrate still separates:

| | engaged | not engaged |
|---|---|---|
| placebo yes (active elsewhere) | 39% of 155 calls | 15% of 72 |
| placebo no | 24% of 41 | 1% of 129 |

A protein that neither moves its direct substrate nor anything off-route gets 1 of 129 change calls right. The two
measurements add to each other; the first hop is not just a stand-in for general activity.

## E3. The gate fails, on the recall guard

Withholding the model's change calls for not-engaged proteins:

| | withheld | fixed | broken | recall | gate |
|---|---|---|---|---|---|
| half A | 57 | 53 | 2 | 0.221 → 0.212 | passes |
| half B | 144 | 122 | 10 | 0.176 → 0.124 | **fails** (recall falls 0.052) |
| pooled | 201 | 175 | 12 | 0.200 → 0.171 | fails |

The registration required both halves. Fixed outnumber broken about fifteen to one (p < 1e-30), but recall is low to
begin with, and on half B the ten broken calls cost more than the registered 0.01. Balanced accuracy: gated 0.396
against 0.378 for the model (+0.018 [−0.002, +0.035]).

## Reading

- **Per the registration: E3 fails, so no gated model is proposed.**
- **What E1 and E2 establish:** whether the protein moved its direct substrate is a strong, replicated indicator of
  whether the model's other calls for that protein can be trusted (36% against 6%), beyond general activity. That is a
  confidence statement, not a correction: for a not-engaged protein the model's change calls are almost always wrong,
  but a few real changes still happen.
- It fits specs/176: the model's errors on this dataset are mostly about whether the protein is active, and one
  measured first hop recovers much of that.
- **Limits.** 85 of 428 scored proteins have a curated direct substrate on this panel. It needs a measurement from the
  experiment. Both halves had been seen before this test; the definitions were fixed before any engagement call.

## For Adam

A possible use that needs no gate: report Lun (and serve over-expression results) by engagement beside the raw score,
as korkut is reported (specs/138), and state the measured precision for each group. Not done here.
