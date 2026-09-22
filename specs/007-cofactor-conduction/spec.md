# Feature Specification: Metabolic cofactors are participants, not conduits

**Feature Branch**: `007-measurement-validity`

**Created**: 2026-09-14

**Status**: Implemented — see Outcome

**Input**: small molecules such as ATP must not act as signal conduits, but
the logic networks must keep them.

## Context

Reactome curates the full chemistry of each reaction, so ATP, ADP, H2O, NAD+
and inorganic phosphate are real curated participants and appear as nodes in
the generated logic networks. They belong there.

As a single shared node, one of these molecules reaches a degree no signalling
protein does: ATP has degree 262 in `DNA_Double-Strand_Break_Repair`, against
16 across MP-BioPath's entire 85-network corpus. A perturbation can therefore
reach an unrelated readout by travelling through the cell's energy currency.
That path exists in the graph and every edge on it is a curated fact, but the
traversal is not biology.

The cost is measured. False change calls are the largest error bucket at
scale: 41.9% of 4,329 wrong predictions across 23,788 curator cases.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A perturbation does not travel through the energy currency (Priority: P1)

Someone solves a pathway with a gene knocked out. A readout on the far side of
the network, connected to the knockout only by a path that passes through ATP,
is reported unchanged rather than changed.

**Why this priority**: this is the entire feature. It addresses the largest
measured error bucket, and every other story here is a guard rail on it.

**Independent Test**: on the wide curator set, paired A/B against the current
behaviour on one shared catalog build, reporting macro-F1, per-pathway net and
the both-arms-converged count. The arm must not lose.

**Acceptance Scenarios**:

1. **Given** a network where the only route from a perturbed gene to a readout
   passes through a cofactor, **When** the pathway is solved, **Then** the
   readout is reported at baseline.
2. **Given** a cofactor participating in a reaction that also has a genuinely
   perturbed input, **When** the pathway is solved, **Then** the reaction still
   sees the cofactor as an input and its AND aggregation is unchanged — the
   cofactor contributes a fold of 1.0, not an absence.
3. **Given** the wide curator evaluation, **When** both arms are compared,
   **Then** the set of perturbed root inputs is identical in every shared case,
   because nothing has been removed from any network.

---

### User Story 2 - The networks stay faithful to what curators recorded (Priority: P1)

Someone else generates logic networks for their own purpose and finds every
participant Reactome records still present, including the cofactors.

**Why this priority**: equal-first with US1, and it is a constraint on *how*
US1 may be satisfied rather than a separate increment. Constitution principle I
makes this non-negotiable: the generator represents pathways as curators
intended, and whether a participant may carry a perturbation is a modelling
decision that belongs in the solver.

**Independent Test**: the generated networks are byte-identical before and
after; no generator code is touched.

**Acceptance Scenarios**:

1. **Given** a generated pathway network, **When** this feature lands,
   **Then** its node and edge files are unchanged.
2. **Given** a consumer other than this solver, **When** they read the
   networks, **Then** they see the cofactors and may model them however they
   choose.

---

### User Story 3 - A caller's own measurement is never silently overwritten (Priority: P2)

Someone running an ATP-depletion experiment supplies ATP as an observation.
Their value is used.

**Why this priority**: it is a correctness guard rather than an accuracy
lever, but the failure is silent, which is the class of defect this repo has
been burned by most (`DS_*` guard rails, DS #13).

**Independent Test**: a solve with an explicit observation on a cofactor
returns that observed value, not baseline.

**Acceptance Scenarios**:

1. **Given** no observation for a cofactor, **When** a pathway is solved,
   **Then** it is held at baseline.
2. **Given** an explicit observation for that same cofactor, **When** the
   pathway is solved, **Then** the observation stands.

---

### User Story 4 - The choice is visible, auditable and reversible (Priority: P3)

Someone reading the code can see exactly which molecules are treated this way,
why, and what it measured — and can turn the behaviour off in one step to
reproduce an earlier result.

**Why this priority**: the decision rests on a list that is a judgement call
about chemistry. An unreadable or unpinned list is how the wrong molecules get
excluded (see Edge Cases).

**Acceptance Scenarios**:

1. **Given** a misspelled mode selector, **When** a solve is attempted,
   **Then** it fails loudly rather than silently selecting a different model.
2. **Given** the list, **When** a reader inspects it, **Then** every entry is a
   release-pinned stable identifier grouped by molecule, and the exclusions are
   stated with their reasons.

### Edge Cases

- **A molecule that looks metabolic but is the signal.** Ca2+, PI(3,4,5)P3,
  PI(4,5)P2, cAMP, cGMP, DAG and I(1,4,5)P3 are second messengers: in a
  signalling pathway they *are* the message. An earlier attempt to exclude the
  whole `SimpleEntity` class cost 163 cases by deleting exactly these. A
  rule-based selection was then validated twice without ever checking
  PI(3,4,5)P3 — the molecule the worst-affected pathway is named after. Hence a
  hand-curated list with stated exclusions rather than a class or degree rule.
- **A molecule whose transfer is itself the regulatory event.** Ubiquitin and
  SUMO are excluded for the same reason.
- **A cofactor with no incoming edge.** It already sits at baseline and never
  moves, so the behaviour is a no-op there; the change only bites where a
  cofactor is produced by a reaction inside the pathway.
- **A cofactor that is the perturbed entity or the readout.** US3 covers the
  first. For the second: holding a cofactor at baseline would pin a scored
  readout flat, so this was left as a stated risk — then measured. Across all
  1,006 distinct curator readouts, **0 are themselves a listed cofactor** (all
  1,006 dbIds resolved, so this is not a lookup artifact). The risk does not
  arise on this ground truth. It would arise on a metabolic pathway, where the
  readout may well be ATP or NADH, and nothing currently detects that.
- **A network where cofactors carry much of the routing.** The effect is
  proportionally larger, in either direction — see Outcome.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST hold each listed cofactor at its baseline
  activity for the duration of a solve, so that it cannot carry a perturbation.
- **FR-002**: A held cofactor MUST remain an input to every reaction it
  participates in. Reaction input sets, and therefore AND completeness, MUST be
  unchanged — the cofactor contributes a neutral fold, not an absence.
- **FR-003**: No logic network, and no generator behaviour, may be altered by
  this feature. The networks MUST keep every curated participant.
- **FR-004**: An explicit observation on a cofactor MUST take precedence over
  the baseline hold.
- **FR-005**: The behaviour MUST be selectable, and an unrecognised selector
  MUST raise an error rather than fall back to a default.
- **FR-006**: There MUST NOT be a selectable mode that computes the cofactor
  set and then has no effect.
- **FR-007**: The cofactor set MUST be identified by Reactome stable
  identifier, pinned to a stated release, and MUST cover every compartment
  variant of each listed molecule.
- **FR-008**: The list MUST be derived from the connected release rather than
  typed, keyed by a chemical identifier (ChEBI) rather than by name or by
  stable id, so that a new compartment variant appears by itself and a retired
  stable id disappears by itself. Neither this list nor the generator's
  `_COFACTOR_STIDS` is required to contain the other: they answer different
  questions (may a perturbation travel through this node, versus may a diagram
  bridge be drawn across it), and the generator's copy is independently known
  to carry stale and mislabelled entries. An audit MUST be able to report the
  difference between them and explain each one.
- **FR-009**: The list MUST exclude second messengers (Ca2+, PI(3,4,5)P3,
  PI(4,5)P2, cAMP, cGMP, DAG, I(1,4,5)P3) and modifier tags whose transfer is
  the regulatory event (ubiquitin, SUMO), and MUST record that these exclusions
  are deliberate.
- **FR-010**: Any comparison between arms MUST be conditioned on the
  experiment being unchanged — same perturbed set, same readout set — so that
  two arms are never compared while answering different questions.
- **FR-011**: The ten-pathway experimental set MUST be reported as neutral and
  MUST NOT be cited as evidence for or against this feature.

### Key Entities

- **Cofactor set**: the collection of Reactome stable identifiers treated as
  non-conducting. Chemistry, not curation: energy and phosphate carriers,
  redox pairs, one-carbon donors, water, bulk ions and dissolved gases.
  Release-pinned; grouped by molecule so a reader can audit it.
- **Conduction mode**: which of the two behaviours a solve uses — hold at
  baseline, or the previous unrestricted behaviour.
- **Conditioned comparison**: a paired A/B restricted to cases where both arms
  perturbed the same entities and scored the same readouts.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: On the wide curator set, the arm does not lose: macro-F1 and
  correct-case count both at least hold, measured on one shared catalog build.
- **SC-002**: The perturbed set is identical in 100% of shared cases,
  demonstrating that no structure was removed.
- **SC-003**: Generated networks are unchanged — zero differing files.
- **SC-004**: An explicit cofactor observation is honoured in 100% of cases.
- **SC-005**: Every rejected mode selector raises an error; none is silently
  accepted.
- **SC-006**: Solver stability is not degraded — the non-converged count does
  not rise.
- **SC-007**: The list's size and shape are pinned by assertion, so widening it
  without updating its stated rationale fails.
- **SC-008**: The generator's list is fully covered by the solver's, checkable
  on demand whenever both checkouts are present.

## Assumptions

- "Baseline" is the same normal-activity level the solver already uses for any
  unobserved node; no new notion of a resting level is introduced.
- Cofactor readouts are not a scored case in either benchmark, so holding a
  cofactor cannot mask a readout that the evaluation would have scored.
- The wide curator set is the decision basis. Four hypotheses reversed sign
  between the ten-pathway set and the wide set during this investigation, so
  the small set is reported but never decisive.
- The list is a Release97 snapshot. A future release may rename or split
  entries; the list is expected to be refreshed per release, and is not
  expected to change much.

## Non-Goals

- Filtering or otherwise altering the logic networks, or any change in the
  generator.
- Extending the treatment to second messengers, or to an entity class such as
  `SimpleEntity`, or to a degree threshold.
- Tuning the list against the evaluation set.
- Patching the one regressing pathway with a pathway-specific rule.

## Outcome

Implemented. `DS_COFACTOR_MODE=inert` is the default; `propagate` restores the
previous behaviour.

**Wide curator set** — 71 pathways, 21,450 scored cases, Release97, one shared
catalog build, conditioned on the experiment being unchanged:

| | baseline (`propagate`) | arm (`inert`) |
|---|---|---|
| correct | 18,082 | **18,119** |
| accuracy | 0.8430 | **0.8447** |
| macro-F1 | 0.8095 | **0.8108** |

Net **+37** on 70 changed predictions: 53 spurious change calls removed against
16 true ones lost (27 UP→NO_CHANGE and 26 DOWN→NO_CHANGE correct; 7 and 7
wrong, plus 2 UP→DOWN). The experiment moved in **0 of 21,986** shared cases —
this is the first arm in the accuracy investigation immune to the root-input
confound, precisely because nothing leaves the network. Non-convergence is
unchanged at 200 in both arms.

Per-pathway, against cofactor share of edges:

| pathway | cofactor share | net |
|---|---|---|
| GPVI-mediated_activation_cascade | 52/269 = 19.3% | **−12** |
| DAP12_interactions | 70/469 = 14.9% | +12 |
| Signaling_by_SCF-KIT | 112/834 = 13.4% | +12 |
| Signaling_by_MET | 67/684 = 9.8% | +18 |
| Netrin-1_signaling | 33/366 = 9.0% | −2 |
| Intrinsic_Pathway_for_Apoptosis | 33/433 = 7.6% | +3 |
| Signaling_by_ROBO_receptors | 20/668 = 3.0% | +6 |

Cofactor density predicts the **magnitude** of the effect, not its sign.

### Why GPVI regresses: it is the UUID silo bug, not a cofactor problem

Root-caused rather than patched, per the non-goals. The readouts are
`VAV2_Rho/Rac_effectors:GTP` and `VAV3_...:GTP`. Under `inert` they sit at
**exactly 1.0000 for every perturbation** — PIK3CA, PTPN11 and SYK alike — so
the readout is not merely losing a route, it is losing every route.

The catalyst of the producing reaction, `R-HSA-442307`, is split into two
**disconnected** uuids:

| uuid | in-degree | out-degree | |
|---|---|---|---|
| `0f995db5…` | 0 | 5 | feeds reaction `442291`; nothing feeds it |
| `22fd0285…` | 1 | 0 | receives `SYK:p-VAV`; goes nowhere |

The signal arrives at one copy of the catalyst and the reaction reads the
other. So the perturbed genes have no working route to the readout through the
catalyst at all, and none of this is caused by cofactor handling — the split
predates it.

The only surviving route ran through the shared nucleotide: knockout → GDP
falls → `VAV2_Rho/Rac_effectors:GDP` falls (GDP is an `assembly` component of
that substrate complex) → the `:GTP` readout falls → DOWN, which matches the
truth. That is the **right answer for the wrong reason**: knocking out SYK does
not deplete cytosolic GDP. Pinning the nucleotides removes the accidental
compensation and exposes the pre-existing silo.

Confirmed by removing GTP and GDP from the list entirely: GPVI recovers
exactly (+12) and the overall net collapses from +37 to +3, because the same
pinning is worth +18 in MET, +12 in DAP12, +12 in SCF-KIT and +6 in ROBO.
The nucleotides stay in the list; the defect to fix is the silo.

**Ten-pathway experimental set: neutral and not evidence.** 13 predicted values
move; no case crosses a class boundary; macro-F1 identical to four decimals.
Those ten networks hold 147 cofactor nodes in total and cannot settle this
question either way. Reported here so the next reader does not mistake its
silence for a negative result.

### Negative result: deleting the nodes is worse (constitution III)

A generator-side filter removing these nodes from the networks was built,
measured and abandoned:

- **−84 cases** on 23,622 conditioned curator cases, macro-F1 0.7835 → 0.7766.
- Only 25% of the losses became unreachable. 137 kept a path and merely
  weakened — 74% of those moving toward baseline, median predicted 0.113 →
  1.028, mostly DOWN → NO_CHANGE.

Deletion removes real structure along with the spurious routing. The branch was
deleted; LNG `main` was never touched. (The −84 and the +37 come from different
pipelines and case sets and are **not** directly comparable to each other; each
is a paired A/B against its own baseline on one shared catalog build.)

### Negative result: excluding the `SimpleEntity` class

−163 cases. PI(3,4,5)P3, cAMP, cGMP, Ca2+, DAG and IP3 are all `SimpleEntity`,
and in a signalling pathway they are the signal. Do not re-attempt a
class-based or degree-based rule.

### Defects found and fixed during implementation

- A third mode, `drop`, computed the cofactor set and discarded it — a no-op
  wearing the name of the behaviour the evidence rejects. Removed; it is now an
  explicitly rejected selector so it cannot return silently (FR-006).
- The baseline hold took precedence over a caller's own observation, so
  declaring ATP inert would have silently discarded an ATP-depletion
  measurement. Fixed, with an assertion that fails against the previous code
  (FR-004).
- **The two repositories' lists had silently diverged, and BOTH were wrong.**
  The generator has carried `_COFACTOR_STIDS` since 2026-09-09, used to keep
  the diagram-bridge pass from asserting that a producer of ATP feeds a
  consumer of ATP. Audited against Release97, **6 of its 13 entries** are stale
  or mislabelled: `R-ALL-29438` commented "PPi" is GTP, `R-ALL-29390`
  commented "Pi variant" is PXLP (pyridoxal 5'-phosphate), `R-ALL-29360`
  commented "ADP variant" is NAD+, and three ids do not exist in the release at
  all. The solver's list, built independently from a Neo4j query by molecule
  name, was correct but incomplete: it silently missed compartment variants
  because the name match does not enumerate them.

  **Correction to an earlier claim in this record:** it was reported that the
  solver's list "omitted NAD+ as an entire molecule". That was wrong — the
  omission was an artefact of a `grep` for `# NAD` that did not match the
  `# NAD+` header. The solver's list always carried NAD+'s nine compartment
  variants. The real finding is the opposite of the one first reported: the
  generator's list contributed nothing valid the solver lacked, and the
  "superset of the generator's list" requirement briefly adopted here was
  itself wrong, because it would have imported four defective entries.

  Both lists are now **derived** rather than typed: every `SimpleEntity` whose
  `ReferenceMolecule` carries one of 28 curated ChEBI identifiers, resolved
  against the connected release. 253 species at Release97. ChEBI identity is
  the only stable handle — names miss compartment variants and match by
  substring, stable ids go stale between releases, and this list hit every one
  of those failure modes.

- The wide-set comparison tool computed macro-F1 over a hard-coded label set
  that matched the ten-pathway dump but not the wide one, silently dropping the
  UP class from the average and reporting 0.5522 where the correct figure is
  0.8095. Labels are now derived from the data. Accuracy and win/loss counts
  were never affected.
