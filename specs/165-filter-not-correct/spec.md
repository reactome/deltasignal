# Feature Specification: Filter, don't correct

**Feature Branch**: `specs/165-filter-not-correct`

**Created**: 2026-11-16

**Status**: Draft, registered BEFORE any scoring

**Input**: Adam, 2026-11-16 ("ok"), after asking whether the low specificity and precision could be improved by combining fixes. The honest answer was that stacking corrections is not shown to be additive. What can be measured is FILTERING: the product already serves truth-blind labels on every call, so what precision does a user get by keeping only the calls the labels vouch for, and at what coverage?

## The labels (all truth-blind on the axes scored)

- **reachable**: the predicted fold is not 1 (a route carries signal).
- **robust**: not fragile under the specs/157 registered ensemble (`specs/157-structural-fragility/runs/labels.tsv`). Computed from rule-variant flips, never from truth.
- **confident**: the specs/146 served confidence ≥ **0.85**, looked up in `src/api/confidence.json` at the record's own features (call, railed, magnitude bin, converged; the record carries no basis, so the walk starts at `no_basis`; min cell from the table). The table was fitted on the curator TUNING pathways, so the curator axis is scored on held-out pathways only.

## The measurement

Units are change calls (predicted down or up) on reachable cases, per axis:

- **korkut** (98 record; a drug set's own markers excluded);
- **experimental** (98 record, valid cases);
- **curator held-out** (98 record, valid cases, tuning pathways excluded).

**Filters**: F0 all; F1 robust; F2 confident; F3 robust AND confident.

**Reported per axis and filter**: precision (measured change in the called direction, among retained calls), coverage (retained / F0), and the chance reference (F0's precision, which a random filter keeps in expectation).

**Gate, per axis**: precision(F3) − precision(F0) > 0 with a 95% bootstrap over perturbation units (korkut drug sets; otherwise (pathway, gene, direction)), 10k resamples, seed 165, excluding zero. A claim ("keeping only robust, confident calls raises precision to X at coverage Y") is made only for axes that pass.

## Requirements *(mandatory)*

- **FR-001**: the scorer is committed and fixture-checked before the record is read; korkut role guard (FIT).
- **FR-002**: labels, threshold, units, filters, and gate are fixed above; nothing is fitted.
- **FR-003**: on a pass, the product's served labels are documented as a usable filter with the measured precision and coverage (docs only); on a fail, the negative is recorded.

## Non-goals

- No new labels, thresholds fitted to outcomes, or engine change; CPPA (no fragility labels exist for it) is out of scope.
