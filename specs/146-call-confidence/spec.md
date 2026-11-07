# Feature Specification: Per-Call Confidence From the Record

**Feature Branch**: `specs/146-call-confidence`

**Created**: 2026-11-07

**Status**: Draft (pre-registered before any measurement)

**Input**: Adam, 2026-11-07: predictions in the UI must "look right and be useful", and dataset analysis must give "the most useful answer possible — for example how certain we are of each result". The engine already serves, with every readout, the signals that predict whether a call is right (fold, direction, railed, basis; the solve serves converged). The 25,000-case validation record says how often calls with each combination of those signals were actually right. This spec measures that, validates that the numbers are honest on cases that did not shape them, and serves them.

## What a confidence is here

For a served call, `confidence` = the empirical share of validation-record calls *in the same stratum* that matched the curators — a measured frequency with provenance, never a model's self-estimate. The claim it makes is exactly: "calls that looked like this were right N% of the time on the named record." Its honesty is testable: on cases the table was not fitted on, the claimed N% must match the observed rate (calibration).

## Strata (fixed before measurement)

Features, all available at serve time with no truth:

- `call`: down / not / up (the engine's served call at the served cutoffs)
- `railed`: the specs/142 flag
- `basis`: direct / proxy / set_members (specs/141)
- `magnitude`: |log2 fold| binned at [0, exactly-baseline), then (0, 1], (1, 3], (3, ∞) — with fold == 1.0 exactly its own bin (an untouched readout, the record's `no_path`-shaped class)
- `converged`: the solve's flag (the readouts request gains an optional `converged` field the client passes from its solve; absent → the stratum marginalizes over it)

Hierarchy for small cells, applied innermost-first until a cell holds ≥ 50 fitting cases: drop `basis`, then `magnitude`, then `converged`, then `railed`; the root is per-`call` accuracy. Each served confidence names the level it resolved at. Laplace smoothing (α = 1) toward nothing further — the hierarchy is the shrinkage.

## Fitting and validation splits

- **Fit**: the frozen record's curator TUNING cases only (the 11 `TUNING_PATHWAYS` of `bench/analysis/holdout_report.py`, imported) — the one curator split that is already "spent". Build `20261026-1349_8b6527e`, solver `05dc867` (`curator_cases.tsv`); `basis` joined per (pathway, key_output) from the specs/141 dump (`bench/analysis/141/bench_readouts.json`); `railed` by the specs/142 function imported from `bench/analysis/score_124.py`; `converged` from the case row.
- **Validate calibration on**: (a) the curator held-out split (with the standing honesty note: SELECT — it chose rules, but never touched this table), and (b) the 849 experimental cases (the independent axis).
- **Gates (H1)**: expected calibration error (ECE, confidence-weighted by bucket size over the resolved strata) ≤ 0.05 on curator held-out AND ≤ 0.10 on experimental; and the ranking must carry information: AUC of confidence for predicting call-correctness > 0.5 with a bootstrap interval (by pathway, 2,000 resamples) excluding 0.5, on each validation split.

## Functional Requirements

- **FR-001**: The fitting script (`bench/analysis/146/`) MUST import the record's own code for every reused definition (TUNING_PATHWAYS, `railed`, the 141 dump for basis) — never re-list or re-implement.
- **FR-002**: The fitted table ships as a versioned, provenance-carrying file (`src/api/confidence.json`: the strata, their counts, the smoothed rates, the fit split, build, solver, date, spec), enforced by a repo test like the scope statement's.
- **FR-003**: `/api/readouts` results valued with activities MUST additionally carry `confidence` (and `confidence_stratum`, the resolved level, so a user can see *why*), computed by table lookup on the served call/railed/basis/fold and the request's optional `converged`; `null` exactly where `call` is `null`; absent on resolution-only responses. Additive only. **Amendment 1 (2026-11-07, after the registered measurement; research.md has the trace): an untouched readout — fold exactly 1.0 — is served `confidence: null` with `confidence_stratum: "untouched"`.** The measured rates for that stratum diverge across the truth axes (curator 0.96, experimental 0.30: the specs/080 lost-route class), so serving a single confidence there would mislead by precision; the client presents "the perturbation does not reach this readout in this network" instead.
- **FR-004**: A parity check holds the engine's lookup equal to the fitting script's answer on every validation-split case (the specs/140–142 pattern).
- **FR-005**: A missing or malformed table refuses readout confidence honestly (field served as `null` with a named `missing_tables`-style note), never a fabricated number; the endpoint itself still works.
- **FR-006**: Unit tests with fixtures that compute: lookup through each hierarchy level, the exact-baseline bin, marginalization without `converged`, provenance completeness on the shipped table, and absence/null rules.

## Success Criteria

- **SC-001**: H1 as stated (both calibration gates and both AUC gates), or the honest negative with traces.
- **SC-002**: A served call carries a confidence whose provenance a user can follow to the record.
- **SC-003**: Existing consumers unaffected; the full API suite passes apart from new assertions.

## Non-goals

- No model self-estimates, no fitted probabilistic model beyond stratum frequencies (isotonic or logistic layers are a later spec if the flat table's calibration fails).
- No per-cell-line or per-dataset context (the record has none).
- No UI work in this repo (the WebsiteAngular consumption follows once served).
- korkut/CPPA confidence: out of scope for the table (their truth is drug-response, a different question than the curator/experimental axes measure); the scope statement's drug-response limitation keeps covering that.
