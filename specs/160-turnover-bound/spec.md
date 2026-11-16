# Feature Specification: Per-protein turnover bound on abundance readouts

**Feature Branch**: `specs/160-turnover-bound`

**Created**: 2026-11-16

**Status**: Draft — registered BEFORE any scoring

**Input**: Issue #262 tier 2, approved by Adam 2026-11-16 ("ok" to fetching the per-protein dataset and registering this). Tier 1 (one median constant) is specs/135's stock rule and is call-neutral at the bench cutoffs (floor 0.70 < 0.85), recorded on the issue.

## The problem

On reachable korkut cases, railed predictions on TOTAL-protein readouts are wrong 221 times and right 7 (2026-11-15 route split). A drug pins a target's activity to zero and the solve propagates that zero to downstream proteins' levels, but the antibody measures abundance. Over the 24 h assay (Korkut 2015 Methods, as specs/135 cites), a protein's level can only approach its new steady state as fast as the protein turns over. One median constant cannot change a call (0.70 is still DOWN). Per-protein turnover can: a long-lived protein cannot fall below 0.85 in 24 h.

## The registered rule

For a readout valued in ABUNDANCE mode (total-protein antibodies; phospho and forms-mode readouts untouched), with predicted steady-state fold f, the predicted MEASURED fold at the assay time is the first-order approach

    L = f + (1 − f) · 2^(−T / t_eff),   T = 24 h

which bounds drops (f = 0 → L = 2^(−T/t_eff)) and rises symmetrically. t_eff is the protein's effective half-life:

- **Source**: Mathieson et al. 2018, Nat Commun 9:689 (open access, CC BY 4.0), Supplementary Data 2, sheet "protein half lives high qual". Per gene: the median of the finite half-lives across the eight HUMAN replicate columns (B cells, NK cells, hepatocytes, monocytes, ×2). Mouse neurons are excluded.
- **Scale**: those are non-dividing primary cells (medians 105–375 h), while korkut's line divides, and growth dilution shortens the effective half-life. The dataset supplies RELATIVE turnover only: t_eff = t_gene × (46 h / median over all genes of t_gene), anchoring the median at the dividing-cell 46 h (Schwanhäusser 2011, the constant specs/135 uses). This is registered and parameter-free; nothing is fitted.
- **Readout → genes**: the readout's own gene(s) from bench/korkut maps. Several genes take the median of their t_eff. A gene absent from the dataset leaves the readout UNCHANGED (counted).
- **Calls**: at the bench's existing cutoffs (0.85 / 1.15). The cutoffs are untouched (specs/137's refuted move is not made).

## User Scenarios & Testing *(mandatory)*

### User Story 1 - The registered measurement (Priority: P1)

Offline on the recorded 98 korkut TSV (`results/cff961d/r98fix_korkut`): recompute each abundance readout's call from its recorded `pred_ui` under the rule; phospho rows keep their recorded call. Gate: balanced accuracy beats the recorded calls with a 95% bootstrap interval over the 78 drug sets (10k resamples, seed 160) excluding zero; raw accuracy, the RAF-set-aside figure, and the per-class anatomy (railed-down fixes vs correct-down calls broken) reported whatever they say. Built-in control, reported beside it: the SAME rule with every protein at the 46 h median (tier 1). It is expected call-neutral, and if the per-protein rule's gain does not exceed it, the gain is not the per-protein information.

**Independent Test**: the scorer is fixture-checked on computed known answers before reading the record.

### User Story 2 - Nothing else moves (Priority: P1)

The curator and experimental axes read forms mode (untouched by construction). No solver change, no cutoff change, no new solves.

### User Story 3 - Serving, only per the gate and Adam's adoption (Priority: P2)

A gate pass implements the rule at the readout layer behind `DS_TURNOVER_BOUND=off|mixed2018` (default OFF: it changes served abundance folds and calls). The half-life table ships with provenance, and adoption as a default is Adam's call, as every rule adoption has been. A gate fail records the negative and serves nothing.

## Requirements *(mandatory)*

- **FR-001**: the derived half-life table (gene, t_human_median, t_eff, n_replicates) is committed with source attribution before scoring; the raw download stays out of the repo.
- **FR-002**: role guard, korkut2015 FIT, nothing fitted; T, the 46 h anchor, and the aggregation are fixed here.
- **FR-003**: the gate, the control, and the anatomy as US1; numbers committed whatever they show.
- **FR-004**: serving only per US3.

## Success Criteria *(mandatory)*

- **SC-001**: the measurement is on the record with intervals, the tier-1 control beside it, and the coverage (readouts with a half-life) stated.
- **SC-002**: no recorded figure on any axis moves unless Adam adopts the flag.

## Non-goals

- No change to the call cutoffs, the solver, or phospho readouts.
- No fitted half-lives, no cell-line-specific turnover, and no time-course model (the long-term direction this would feed).
- No CPPA scoring in this spec (its terms keep numbers out of the repo; a verdict-only check would be its own step if the gate passes).
