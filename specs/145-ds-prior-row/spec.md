# Feature Specification: DeltaSignal as the Unseen Drug's Response Row

**Feature Branch**: `specs/145-ds-prior-row`

**Created**: 2026-11-07

**Status**: Draft (pre-registered before any measurement)

**Input**: The next step after specs/144 (Adam 2026-11-07: "do well with analyzing datasets"; proceeding overnight on the standing proceed-on-judgment rule). specs/118 found DS's bare direction call adds a tiny but all-drugs-positive increment to the best published korkut predictor, and recorded "CellBox with DS prior" as the next option. specs/144 showed conditioning can only delete false calls; the open problem is supplying *correct change calls for a drug the fitted model has never seen*. The specs/118 harness makes the test exact: in its leave-one-drug-out, the held-out drug's design column is zeroed, so the fitted model predicts **zero direct effect** for it — whatever drug-specific signal appears must come from outside the data. DeltaSignal predicts a full response row for any drug from structure alone. This spec injects that row.

## Hypotheses

- **H1 (the prior works)**: adding the held-out drug's DS-predicted response row — its frozen single-drug clipped log2 folds per readout, scaled by one coefficient fitted on the *training* drugs only (how well DS's rows match their fitted rows) — raises the LODO mean per-drug correlation over the plain LR, 95% drug-bootstrap on the difference excluding 0.
- **H1′ (the structure is load-bearing)**: the DS row beats the same construction built from reachability's calls (the specs/113 lesson: a null that eats the gain means the wiring, not the propagation, carried it).
- **H2 (engagement composes)**: conditioning the injected row by the specs/144 dial (α from the held-out drug's own markers; markers excluded from scoring in **all** arms of this comparison) improves on the unconditioned row.

## Protocol (fixed before measurement)

- **Harness**: `bench/analysis/lr_118.py` — `load`, `fit`, `corr`, `calls`, the LODO split and the drug bootstrap — imported, never re-implemented. `dataset_roles.require_role("korkut2015", "fit")` as every fitting script must.
- **DS rows**: from the frozen record (`a490cdc/k134_abund_a3/korkut_cases.tsv`), the single-drug condition's per-readout mean over pathways of clipped log2 `pred_ui` — exactly the harness's existing `ds_fold` feature matrix at the single-drug rows.
- **Arms**:
  1. `lr` (control) and `lr+ds_fold` (the 118 arms, reproduced for continuity).
  2. `lr+dsrow`: prediction for held-out drug j's conditions = kept-column LR prediction + D[:, j] · (c · dsrow_j), with c = argmin Σ_{k∈train} |R_k − c · dsrow_k|² over the training drugs' fitted rows. No quantity touching drug j's measurements enters the fit.
  3. `lr+reachrow`: identical construction from `reach_predicted` calls (±1/0 rows) — the H1′ control.
  4. `lr+dsrow×dial`: dsrow_j scaled further by α_j, the specs/144 pooled estimator on drug j's own markers (`fit_alpha`, shared from `bench/analysis/144/engagement_dial.py`, never re-implemented); for this comparison (arm 4 vs arm 2) marker readouts are excluded from scoring in both arms.
- **Metrics**: LODO mean per-drug Pearson correlation (the published convention) is primary; single-drug-condition correlation and balanced accuracy at log2 1.5 reported beside it; per-drug table always shown (concentration check).
- **Intervals**: the harness's drug bootstrap (seed as shipped), 2,000 resamples.
- **Traces before any claim**: which drugs gain/lose; whether c is stable across folds; whether the dial's α values match specs/144's.

## What this is and is not

- korkut2015 is FIT: development, honestly measured. A pass motivates the same construction on HPN-DREAM's harness (already-scored FIT data) and ultimately a sealed test (Lun stays sealed; Adam's alone).
- This is the dataset-analysis product shape in miniature: a user's dataset supplies the fitted part; DeltaSignal supplies the never-measured perturbation's row.

## Success Criteria

- **SC-001**: H1 and H1′ as stated (both must hold for a "the prior works" claim; H1 without H1′ is recorded as "the wiring carried it").
- **SC-002**: H2 as stated, or its honest negative.
- **SC-003**: Negative results committed with full traces.

## Non-goals

- No solver or bench-scoring changes; no new datasets; no CellBox retraining (that is the follow-up if the LR form shows the prior is material); no CPPA numbers in the repo.
