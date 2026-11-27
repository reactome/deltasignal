# specs/176: why Lun 2019 half A scores so poorly — a failure budget, and what is fixable

**Status:** registered 2026-11-26, before any of the tests below was computed.
**Requested by:** Adam, 2026-11-26 ("can you figure out why the dataset is failing so poorly. I feel like we should be
able to fix this"). Diagnosis on the development half; any fix it suggests is registered separately and judged once on
half B. **Half B is never read.**

## Starting point (specs/173–175)

Half A, 3,127 pairs: 744 measured changes, of which the model catches 143 in the right direction (recall 0.19); 405
of its 548 change calls are wrong. Balanced accuracy 0.384, reachability 0.380. specs/174 classified the cases; this
spec asks how much each cause costs and whether it can be fixed.

## Registered tests (half A only; FIT)

**T1. Is the measured call reliable?** Each construct has one basal well, and each plate pair has about two FLAG-GFP
control wells.
- a. The same construct's EGF-10-minute well, called the same way (controls at 10 minutes), as a second measurement:
  share of basal change calls that repeat (same direction), by readout.
- b. Effect size: pairs binned by |E| / δ (1–1.5, 1.5–2, 2–3, > 3). Route enrichment, model precision and repeat rate
  per bin, and the specs/173 scores recomputed with the call at 2δ and 3δ.
- c. Bin choice: the call from bins (1,2,3) vs (6,7,8) against the registered (0,1,2) vs (7,8,9).

**T2. How much of the truth is "which protein" and "which readout", not "which route"?** A two-factor score, fitted
leave-one-out: (the protein's change rate on its other readouts) × (the readout's change rate under other proteins).
AUROC for "changed" of that score, of route (reach), and of the model's call, on the same pairs.

**T3. Missing routes.** For the measured changes with no route (N1): is the protein in the readout's network at all
other than as a pinned root; does any form of the protein reach the readout in any network; and for receptor tyrosine
kinases, where does the canonical route (receptor → RAS/ERK, PI3K/AKT, STAT, PLCγ) break?

**T4. Why does p-AKT fall?** If chronic signalling shuts AKT S473 down through S6K feedback, the proteins that lower
p-AKT should raise p-S6 / p-p70S6K / p-mTOR more often than proteins that do not. A 2×2 table, Fisher's exact test.

**T5. What do the known false-call sources cost?** Model scores on half A with change calls withheld (a) when the
shortest route passes through PIP3, (b) when it is longer than 12 steps, (c) both. Reported as fixed / broken with
recall, for sizing only; no rule is adopted here.

**T6. Is the over-expressed protein active at all?** Per protein, the share of all its measured readouts (the 8,640
half-A pairs, not only those in P) that change. Model precision and recall for proteins that change nothing else
(leave the scored pair out) against proteins that change something.

## Reading

Descriptive. The output is a budget: of the 744 measured changes and 405 wrong calls, how many fall to measurement
noise, to protein-level activity, to missing routes, to the switch, to long routes, and what is left. Candidate fixes
are listed with their half-A size, for separate registration.
