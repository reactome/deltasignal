# Feature Specification: The Engagement Dial — Conditioning Predictions on the Dataset's Own Markers

**Feature Branch**: `specs/144-engagement-dial`

**Created**: 2026-11-07

**Status**: Draft (pre-registered before any measurement)

**Input**: Adam, 2026-11-07: "keep moving forward with trying to make it so that we can do well with analyzing datasets", "we need to figure out what the next thing is". The record's answer: the model fails on drug-response datasets for traced cell-context reasons (specs/062/138: five korkut drugs never engaged their targets in these cells; drug target sets alone give 804 of 1,408 wrong change calls), three attempts to *predict* that context failed (specs/121/130/139), and specs/117's round-2 list registered — but never executed — "the drug's inhibition strength, fitted as a solver pin". Every fitted method that beats us on korkut (Long et al. LR r 0.784, CellBox 0.780; specs/118) conditions on the dataset's measurements. Analyzing a dataset means conditioning on it.

## Hypothesis (H1)

For each drug set, one engagement parameter α — fitted only on the drug set's own direct markers (the specs/138 mechanism-named markers, fixed before any of this) — improves direction calls on all *other* readouts: conditioned fold = fold^α. A drug whose own markers did not move gets α ≈ 0 and honestly predicts quiet downstream; an engaged drug keeps its prediction. H1: conditioned balanced accuracy over non-marker scored readouts exceeds the unconditioned record's, with a 95% by-drug-set bootstrap interval on the difference excluding 0, at the registered primary measured threshold (|log2| ≥ log2 1.5).

## Protocol (fixed before measurement)

- **Control**: the frozen record's korkut run — build `20261026-1349_8b6527e`, solver `a490cdc`, arm `k134_abund_a3` (`korkut_cases.tsv`, 11,626 rows, raw accuracy 0.812 on 10,469 scored+measured cases — verified against the file before registration). No re-solve; this is an output-side measurement model, like specs/130 in mechanism but conditioned per drug set on its own markers rather than fitted across drugs.
- **Markers**: `bench/analysis/138/engagement.py MARKERS` — imported, never re-listed. A drug set's markers are the union over its member drugs. Marker *readout rows are excluded from scoring in both arms* (identical denominators, the standing rule).
- **Measured marker folds**: mean log2 over the drug set's own condition rows in `expr.csv` (both doses), per marker, via the bench's own `conditions()`/`node_Index` readers.
- **α estimator**: least squares in log space over the drug set's markers, α = Σ(log m̂ · log f̂) / Σ(log f̂)², clipped to [0, 1], where f̂ is the model's predicted marker fold and m̂ the measured one.
  - Primary granularity: per (drug set, pathway), fitted on the markers that pathway scores.
  - Fallbacks, in order: the drug set's pooled α over all pathways' scored markers; α = 1 (no conditioning). Every fallback use is counted and reported.
  - Σ(log f̂)² below 1e-8 (the model predicts no marker signal) is undefined → fallback.
- **Conditioned call**: fold^α classified at the bench's served cutoffs (0.85 / 1.15), exactly `classify`.
- **Metrics**: balanced accuracy (mean per-class recall over DOWN/NORMAL/UP) at the primary threshold is the gate; raw accuracy and the secondary threshold (log2 1.15) reported beside it; RAF set aside and headline both reported (specs/054 D1); always-no-change stated alongside.
- **Interval**: 10,000 bootstrap resamples over drug sets of the arm−control difference.
- **Traces before any claim** (the standing rule): the five non-engaged drugs' changed calls; the largest fixed and broken classes; whether gains concentrate in a pathway or drug.
- **Independent check**: the CPPA training drugs, same recipe (its own direct markers), verdict recorded in the repo, numbers kept private under the data agreement.

## What this is and is not

- korkut2015 is FIT (specs/051/062): this is development, measured honestly, not a generalization claim. A passing gate motivates stage 2 (the mechanistic version: a fitted partial-inhibition *pin*, re-solving, so engagement propagates structurally) and, later, a sealed test (Lun 2019 remains sealed; opening it is Adam's alone).
- The fit uses only the drug set's own markers, which are excluded from scoring — nothing about the scored readouts leaks into α.
- This is the "analyze a dataset" mode in miniature: anchor on what the dataset measured, predict the rest.

## Success Criteria

- **SC-001 (the gate)**: H1 as stated.
- **SC-002**: The non-engaged drugs' false change calls fall, traced case-level, without a compensating loss class elsewhere (McNemar-style fixed/broken counts reported).
- **SC-003**: Negative or null results are committed with the same completeness as a pass.

## Non-goals

- No solver change (stage 2, separately registered, if stage 1 passes).
- No new datasets, no catalog build, no change to the bench's scoring files.
- No CPPA numbers in the repo.
