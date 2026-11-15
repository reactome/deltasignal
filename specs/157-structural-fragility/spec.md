# Feature Specification: Structural-fragility flag

**Feature Branch**: `specs/157-structural-fragility`

**Created**: 2026-11-15

**Status**: Draft

**Input**: Improvement backlog #254 (the 2026-11-13 approved list, item 4; items 1–3 shipped end to end, item 3 on 2026-11-15).

## The problem

Roughly twenty solver rules were measured against near-equivalent alternatives on the record (the specs/033–109 arms: depletion modes, carrier-loop variants, pool weights, set-pool aggregation, free-ubiquitin handling, metabolite pools, catalyst-loop holds, self-interaction weights), and the adopted choice often won by a handful of cases — or was adopted on faithfulness with no case moving at all. A served call that FLIPS when one of those contested choices is toggled rests on a single modelling decision the record could barely distinguish; a call that survives the whole ensemble is structurally robust. Today the product serves both identically. The specs/146 confidence rates a call's CLASS against the validation record; this flag rates the CALL's dependence on contested structure — orthogonal, and cheap, because every variant is a recorded `DS_*` setting and the arms re-run in hours.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - The registered ensemble (Priority: P1)

The exact variant list is a versioned registry in the repo (variant name, `DS_*` setting, the spec that measured it, its recorded delta), resolved in research from the record and FROZEN before the measurement runs. The inclusion rule, registered up front: each variant is a recorded measured arm that is either (a) a not-adopted near-equivalent — its gate failed with a small aggregate delta (candidates: specs/060 `suppress_only_entity`, specs/090 pool supply AND, specs/093 carrier moves, specs/098 inert ubiquitin, specs/103 depletion share, specs/104 drive twins, specs/107 metabolite pools, specs/109 catalyst loops, specs/136 arm L, specs/033 geomean/mean set pools) — or (b) the immediately-preceding default an adopted rule replaced (candidates: specs/075 full re-evaluation, specs/077 `fed_single`, specs/083 off, specs/091 full base weight, specs/095 off). REFUTED arms with large recorded losses (specs/063 sequestration, specs/056 identity context, specs/040 rule A, set-pool max) are EXCLUDED: a flip against a known-wrong rule is not fragility.

**Why this priority**: the ensemble's membership IS the claim's meaning; an unregistered or shifting list makes "robust" unfalsifiable.

**Independent Test**: the registry file exists with every row citing its measuring spec and recorded delta; the research record states, for each candidate considered, in or out and by which clause.

---

### User Story 2 - The registered measurement: does fragility predict error? (Priority: P1)

Before any serving decision, every recorded benchmark case (curator and korkut, on the served catalog) is solved under every ensemble variant through the served stack. Each (case, readout) is labelled ROBUST (call identical across all variants, the canonical included) or FRAGILE (any flip). The gate: on at least one axis, the robust class's accuracy exceeds the fragile class's, with a 95% bootstrap interval (over perturbations on curator, over drug sets on korkut) excluding zero. The specs/154-style informativeness statement — how much of the record is fragile, where the flips concentrate, which variants drive them — is recorded whatever it says.

**Why this priority**: the serving decision rests on it; a flag that does not separate right from wrong answers serves only noise, and the honest outcome is recording that.

**Acceptance Scenarios**:

1. **Given** the ensemble run completes, **Then** per-axis robust/fragile accuracies with bootstrap intervals are recorded in research.md whatever they show, with per-variant flip counts.
2. **Given** a gate failure on every axis, **Then** nothing is served; whether an unvalidated flag is worth serving is Adam's decision, recorded as his call — never a silent default (the specs/156 US3 scenario-2 pattern).

---

### User Story 3 - Serving, only per the measurement (Priority: P2)

If the gate passes, the flag serves per readout. The serving shape is resolved in research from the code: either a per-readout-class fragility-rate table shipped with provenance exactly like the specs/146 confidence table (the stratum-walk, min-cell, provenance-refusal pattern the codebase already carries twice), or a live per-request ensemble — the latter REJECTED if it needs per-request env mutation in the served container (the solver resolves its config from process ENV; research must check whether a variant config can be passed per-solve without touching ENV). Honest absence everywhere the record is silent: untouched readouts, readouts without a call, classes below a registered minimum cell.

**Independent Test**: a readout in a measured class serves the flag with its provenance; a readout outside every class serves nothing; a table fitted for another catalog is refused, never substituted.

---

### Edge Cases

- A case whose canonical call is NOT an error but flips under one variant: fragile — fragility is about dependence on contested structure, not about being wrong.
- A variant that fails to converge on a case where the canonical converges: the honest label is resolved in research before the run (registered there, not decided after seeing results).
- A readout served through a composite or lead (specs/149): out of scope for the flag (the lead class already carries its own rate).
- The ensemble's variants interact with the serving flag `DS_SOFT_OBS` and other serving env: each variant server runs the canonical serving env PLUS exactly one variant setting, verified positively in-container.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001 (registry frozen first)**: the ensemble registry is a versioned table in the repo (variant, `DS_*` setting, measuring spec, recorded delta, inclusion clause), committed before the measurement runs; the measurement cites its commit.
- **FR-002 (the measurement protocol)**: the ensemble runs on the served catalog through replay servers with each variant's env verified POSITIVELY in-container (the specs/156 lesson: never inferred from the env alone), pinned to a pushed commit; per-case-per-variant calls committed under the spec whatever they say; role discipline — curator held-out is SELECT, korkut is FIT, nothing is fitted.
- **FR-003 (the gate)**: as US2 — robust accuracy > fragile accuracy on at least one axis, 95% bootstrap over the axis's independent units, excluding zero; both axes always reported; RAF set aside read beside the headline per the standing rule.
- **FR-004 (additive serving, if served)**: no existing field changes; the specs/141/142 dumps replay unchanged; provenance refusal on catalog mismatch; honest absence below the registered minimum cell and for untouched readouts.
- **FR-005 (tests with fixtures that compute)**: the lookup (class present/absent, minimum-cell refusal), the provenance refusal, additive parity on the existing fixture; registered in the asserting list and the CLAUDE.md count table.
- **FR-006 (honest framing)**: the served flag's wording claims stability across recorded near-equivalent structure, never correctness; the robust/fragile accuracy figures travel with it.

## Success Criteria *(mandatory)*

- **SC-001**: the registry is on the record with every inclusion/exclusion justified by its clause before any measurement number exists.
- **SC-002**: the measurement is on the record with intervals and per-variant flip anatomy, whatever it shows.
- **SC-003**: the serving decision cites the gate; a fail serves nothing without Adam's recorded call.
- **SC-004**: if served, every existing figure and field is byte-identical with the feature on, and a user can tell a robust call from a fragile one with the measured rates attached.

## Assumptions

- The ~20 variants re-run in hours on the served catalog via the sharded replay-server pattern (the specs/156 arm took ~3 h for 31k solves across 4 servers; the ensemble is ~20 × the benchmark's solve count, so sharding and per-variant servers are expected to be needed, up to the 4-server limit).
- The recorded benchmark case lists (curator + korkut on the served 98 record) are the measurement population; no new cases are created.
- WebsiteAngular consumption is its own spec, only if the gate passes and the flag serves.

## Non-goals

- No new rule adoption or re-litigation of any recorded arm decision — the ensemble reads the record, it does not reopen it.
- No change to any call, fold, confidence or scoring; no CPPA data anywhere.
- No per-user ensemble configuration; no UI work in this repo.
- No claim that robustness means correctness — it is stability across recorded near-equivalent structure, stated as such.
