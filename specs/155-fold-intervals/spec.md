# Feature Specification: Measured fold intervals beside the point fold

**Feature Branch**: `specs/155-fold-intervals`

**Created**: 2026-11-14

**Status**: Draft

**Input**: Improvement backlog #252 (the 2026-11-13 improvement list, item 2; Adam: "I would keep this list somewhere so that you can work your way thought it"). Item 1 (switch states, specs/154) shipped end to end 2026-11-14; this is the next item in the approved order.

## The problem

The product serves a point fold — "0.40×" — as if it were a measured magnitude, but no magnitude claim has ever been validated: specs/021 defines what a magnitude claim can and cannot be, specs/124 showed magnitudes move with the model's ceiling while calls barely do, and specs/128–131 found the magnitude signal too weak to fit response curves or edge weights. The served confidence (specs/146) rates the CALL, not the number. A user reading 0.40× has no way to know whether measured reality lands anywhere near it. The honest, record-backed answer is an interval: "when the model answered like this, the measured fold landed between X× and Y× half the time, over n recorded cases" — the same epistemic footing as the served confidence, applied to the magnitude.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A prediction carries what the record measured for answers like it (Priority: P1)

A user solves a perturbation and reads a readout's predicted fold. Beside it, the product serves a measured range: for each recorded axis with enough pairs in this answer's class, the band the measured fold actually landed in (at a stated nominal share, e.g. the middle 50%), the number of recorded pairs behind it, and which record it comes from. The number stops pretending to be a measurement; the record speaks instead.

**Why this priority**: it is the whole feature — the served point fold is the last headline number without a measured honesty statement attached.

**Independent Test**: request readouts for a solved benchmark perturbation; a readout with a call whose class has recorded pairs carries the range object(s) with the stated share, n and axis; the values match the shipped table exactly (no server-side computation beyond lookup).

**Acceptance Scenarios**:

1. **Given** a solved perturbation and a readout with a call whose class has at least the minimum recorded pairs on an axis, **When** readouts are requested, **Then** the entry additionally carries, per such axis: the interval's low and high fold, the nominal share it covers, the recorded pair count, and the axis name — and every existing field is unchanged byte for byte.
2. **Given** the same request with the feature disabled or the table absent, **Then** the response is exactly today's, and the protocol says why no ranges are served.

---

### User Story 2 - Honest absence everywhere the record is silent (Priority: P1)

A readout the perturbation does not reach (stratum "untouched"), a readout without a call, a class with fewer recorded pairs than the registered minimum, an axis with no recorded pairs for the class, or a shipped table whose provenance does not match the served record — each serves no range at all, never a guess, and the protocol states the rule.

**Why this priority**: the feature's value is honesty; an interval extrapolated past the record is worse than a bare point fold.

**Independent Test**: fixture requests covering each absence case return entries without the range fields; a table with mismatched provenance is refused with the refusal named in the protocol, exactly as the specs/146 confidence table is.

**Acceptance Scenarios**:

1. **Given** an untouched readout or one without a call, **When** readouts are requested, **Then** no range fields are present on that entry.
2. **Given** a class below the registered minimum pair count, **Then** no range is served for it on that axis.
3. **Given** a shipped table whose recorded provenance does not match the serving record, **Then** no ranges are served anywhere and the protocol names the mismatch.

---

### User Story 3 - The calibration is measured before it is served, with registered gates (Priority: P2)

Before anything ships, a registered measurement decides WHAT may be served. The pairs are the recorded (predicted fold, measured fold) cases of the public fitting datasets with measured magnitudes — korkut 2015 (RPPA fold changes, the specs/134 abundance rule) and decryptM (phospho fold changes) — per axis, never pooled across axes without the axis named. CPPA is excluded entirely (its terms keep data and derived numbers outside this repo). Two conditioning sets are registered up front:

- **(a)** the measured fold's quantiles conditional on the answer's class alone (the specs/146 confidence stratum, which already encodes direction, size band, resolution and convergence) — what the record says about answers like this, regardless of the predicted number;
- **(b)** additionally conditional on the predicted-fold band through a monotone binning — narrower intervals, served only if the predicted magnitude actually carries information.

The gates are measured out-of-sample by WHOLE drug (leave-one-drug-out; the specs/051 discipline applied within fitting data, since no magnitude-bearing SELECT or TEST set exists):

- **Gate 1 (serve at all)**: the nominal-share intervals must be empirically valid out-of-sample — the held-out coverage's 95% bootstrap interval (resampling drugs) must include the nominal share, per axis.
- **Gate 2 (serve (b) over (a))**: (b) must pass gate 1 AND be materially narrower than (a) out-of-sample — registered: median interval-width reduction of at least 10% on the held-out pairs.

If (b) fails gate 2, the honest outcome is serving (a). If (a) fails gate 1 on an axis, that axis serves nothing. If every axis fails, the feature ships nothing and the research record says why — that too is a finding (the point fold would be proven uninformative of magnitude even as a class statistic).

**Why this priority**: P2 only because it precedes serving rather than being the served behaviour; nothing in US1 ships without it.

**Independent Test**: the fitting script runs from the recorded benchmark results alone, proves it reads only FIT-role datasets, writes the per-axis gate verdicts and the table, and a report states coverage and widths per axis and arm.

**Acceptance Scenarios**:

1. **Given** the recorded korkut and decryptM scorings, **When** the fitting script runs, **Then** it refuses any dataset not registered as FIT, and writes gate verdicts measured by leave-one-drug-out with the registered statistics.
2. **Given** a gate failure, **Then** the corresponding arm (or axis, or the whole feature) is not served, and the research record states the measured numbers.

---

### Edge Cases

- A class with recorded pairs on one axis and not the other serves the covered axis only.
- A readout whose predicted fold falls outside every recorded fold band under conditioning (b) serves the class-level (a) interval for that axis instead of extrapolating — the band structure is a refinement, never an extrapolation licence.
- The table is versioned against the record it was fitted on; serving a different catalog's record refuses the table (specs/143/146 pattern), so a stale table can never quietly describe a new model.
- Intervals are shares of the record, not probabilities, and the served wording must say "of recorded cases", never "chance" or "probability".

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001 (additive only)**: with the feature on, every existing response field is byte-identical; the specs/141 readout dump and specs/142 rails dump replay unchanged through the feature branch.
- **FR-002 (registered gates)**: the gates of US3, exactly as stated (gate 1: out-of-sample coverage valid per axis, 95% drug-bootstrap including the nominal share; gate 2: (b) valid and at least 10% median width reduction over (a) out-of-sample; minimum class size: 30 recorded pairs), are measured BEFORE any serving decision and their numbers recorded in research.md whatever they are.
- **FR-003 (provenance, the specs/146 pattern)**: the shipped table records what it was fitted on (catalog build, solver commit, fitting spec, axes, pair counts); the server refuses a table that does not match the serving record, and the response protocol names the table's version or the refusal.
- **FR-004 (honest absence)**: no range for untouched readouts, call-less readouts, classes under the minimum pair count, axes without pairs, or a refused/absent table; absence is field absence, never a null-shaped guess.
- **FR-005 (tests with fixtures that compute)**: unit tests cover the lookup (present and absent classes), the minimum-n refusal, the provenance refusal, field absence per FR-004, and FR-001 on the existing fixtures; the suite is registered in the asserting list and the documentation table.
- **FR-006 (role guard)**: the fitting script calls the dataset-role guard so it can only read FIT-role datasets; CPPA appears nowhere in its inputs or outputs.
- **FR-007 (axis named)**: every served interval carries the axis it was measured on; nothing is pooled across axes.

### Key Entities

- **Calibration table**: per (axis, class[, predicted-fold band]): lo fold, hi fold, nominal share, n pairs; plus provenance (record identity, fitting spec, date) and the gate verdicts it shipped under.
- **Pair**: one recorded benchmark case with both a predicted fold and a measured fold on a magnitude-bearing axis, keyed by (axis, drug/perturbation, readout).
- **Range object**: the served per-readout, per-axis statement: {axis, lo, hi, share, n}.

## Success Criteria *(mandatory)*

- **SC-001**: on the served catalog, every benchmark readout entry with a call and a covered class carries at least one measured range, and spot-checks against the shipped table match exactly.
- **SC-002**: replaying the full recorded readout and rails dumps with the feature on produces zero differences in existing fields (24,356 folds).
- **SC-003**: the serving decision (nothing / (a) / (b), per axis) is exactly what the registered gates measured, with the numbers in the research record before the decision.
- **SC-004**: out-of-sample, the served intervals' empirical coverage includes the nominal share within the registered bootstrap bound, per served axis.
- **SC-005**: a user-facing statement can be generated from served fields alone ("on the drug-RPPA record, answers like this landed between 0.3× and 0.9× half the time, over 214 cases") with no client-side computation.

## Assumptions

- The recorded korkut and decryptM scorings (both FIT since their one scoring) contain per-case predicted and measured folds sufficient to form pairs without re-running benchmarks; if a re-run is needed it uses the frozen protocols on the served catalog and changes no recorded figure.
- The nominal share is 50% (the middle 50% of recorded measured folds); a second, wider share (90%) may be shipped alongside under the same gates.
- korkut log2 measurements convert to folds for serving; the axis label makes the record's protocol unambiguous.
- The WebsiteAngular consumption of the range fields is its own spec once served, like specs/154 → WA specs/026.

## Non-goals

- No change to any call, fold, railed flag, confidence or scoring; no benchmark figure moves.
- No CPPA data or derived numbers anywhere in this repo.
- No refitting of solver parameters (specs/128–131 closed that at current data).
- No cross-axis pooling that hides the axis; no "probability" language.
- No UI work in this repo; no new catalog build.
