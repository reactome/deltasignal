# Feature Specification: Long-range upregulation propagation

**Feature Branch**: `002-upregulation-propagation`

**Created**: 2026-09-10

**Status**: Draft

**Input**: Close the accuracy gap against MP-BioPath on root-input single-gene
perturbation cases. DeltaSignal must at minimum match MP-BioPath.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A knocked-in gene raises its downstream readout (Priority: P1)

A researcher overexpresses a gene at a pathway's root input and asks what
happens to a terminal output several reactions downstream. Today the readout
frequently reports no change at all, even when a directed path exists and the
curated biology says it should rise.

**Why this priority**: This single case class accounts for essentially the
entire accuracy gap. Downregulation already works — knockouts propagate
correctly — so nothing else needs to change first, and fixing it in isolation
is independently demonstrable.

**Acceptance Scenarios**

1. **Given** a gene at a root input and a terminal readout with a directed
   path between them, **When** the gene is overexpressed, **Then** the readout
   must move away from baseline in the direction the curators recorded.
2. **Given** the same network and a knockout instead, **When** the case is
   scored, **Then** downregulation accuracy must not fall below its current
   level. The fix must not trade one direction for the other.
3. **Given** a path that crosses one or more complexes, **When** one subunit
   is raised and the others stay at baseline, **Then** the complex must
   register some increase rather than reporting exactly baseline.

### User Story 2 - Results are trustworthy against a null model (Priority: P2)

Anyone reading a DeltaSignal result needs to know it beats guessing from
topology alone. Today a sign-of-shortest-signed-path traversal with no
propagation model scores higher than DeltaSignal on the same cases.

**Why this priority**: Matching MP-BioPath is the headline goal, but beating
the structural baseline is the floor beneath it. A propagation model that
loses to "follow the arrows and count signs" has not earned its complexity.

**Acceptance Scenarios**

1. **Given** the ten evaluation pathways with root-input perturbations,
   **When** DeltaSignal and the shortest-signed-path baseline are scored on
   identical cases, **Then** DeltaSignal must exceed the baseline.

### Edge Cases

- A complex whose subunit is knocked out must still collapse. Whatever
  replaces the current hard clamp must keep scarcity propagating, because
  that is the half that currently works.
- Cases with no directed path (~21% of scored cases, unchanged between arms)
  are out of scope here; they are a network-coverage problem, not a
  propagation one, and must not be counted as either a win or a loss.
- A readout reached by several paths of different lengths must not be
  dominated by the longest, most attenuated one.
- Perturbations that legitimately produce no change must still produce no
  change; loosening the clamp must not manufacture spurious movement.

## Requirements *(mandatory)*

### Functional Requirements

- **FR1** — An increase applied at a root input MUST remain detectable at a
  terminal readout across the path lengths present in the evaluation
  pathways, unless the curated biology says otherwise.
- **FR2** — A complex MUST NOT report exactly its baseline when one of its
  subunits is elevated and the rest are at baseline.
- **FR3** — A complex MUST still fall when one of its subunits is removed.
  Scarcity must continue to propagate.
- **FR4** — Every change MUST be measured against MP-BioPath and the
  structural baseline on the identical scored case set, per the project's
  constitution ("measure, then claim").
- **FR5** — Any change that alters predictions MUST report per-pathway
  attribution and a both-arms-converged count, because uuid relabelling makes
  non-converged solves differ between structurally identical networks and has
  twice produced spurious deltas in this project.
- **FR6** — Negative results MUST be recorded with their numbers rather than
  discarded.
- **FR7** — The evaluation MUST report a stable primary figure. Either TP53's
  non-convergence is resolved, or the primary comparison is restricted to
  converged solves and the non-converged remainder reported separately. A
  number that moves by 39 cases on rebuild cannot support any claim.

### Key Entities

- **Case** — one (pathway, gene, direction, readout) tuple from the
  MP-BioPath experimental set, with an expected class of DOWN, NO CHANGE or UP.
- **Root input** — a network node with no incoming edge; the only place a
  perturbation may be applied, matching how the cases were designed.
- **Complex aggregation rule** — how a bound species' level is derived from
  its subunits. Currently a hard minimum; this is the primary suspect.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC0** *(prerequisite)* — Rebuilding the catalog and re-running an
  unchanged arm reproduces the headline figure to within a handful of cases.
  Until this holds, SC1-SC6 cannot be evaluated.
- **SC1** — DeltaSignal matches or exceeds MP-BioPath on the identical scored
  case set under root-input perturbation. This is the goal the tool is
  judged on.
- **SC2** — DeltaSignal exceeds the sign-of-shortest-signed-path baseline on
  the identical scored case set (currently 393/564; DeltaSignal is 295 with
  the clamp on, 368 with it off).
- **SC3** — Correct upregulation calls rise from 76 of 246 toward the
  199-of-279 rate observed when the signal did not have to travel.
- **SC4** — Cases where a directed path exists but the readout sits at exactly
  baseline fall from 49 toward 0.
- **SC5** — Correct downregulation calls do not fall below 173 of 247.
- **SC6** — Macro-F1 improves, not just accuracy, so that "predict no change"
  is not rewarded.

## Assumptions

- The ten MP-BioPath experimental pathways remain the evaluation set, scored
  against experimental ground truth at the harness's default 0.85/1.15
  cutoffs.
- Perturbations are applied only at root inputs. Cases whose gene has no
  root-input occurrence are reported unscored rather than approximated.
- Reactome Release97 networks generated from current LNG `main`.
- The ~21% of cases with no directed path are a separate, structural problem
  and are excluded from the propagation judgement.
- MP-BioPath's per-case predictions are fixed published values and are not
  recomputed.

## Blocking Discovery: the instrument is not yet trustworthy *(2026-09-10)*

Re-establishing the baseline after a machine crash produced **334/564** where
the pre-crash run of the same arm produced **295/564** — a 39-case swing with
no generation-affecting change on either side (LNG differs only by credential
redaction; DeltaSignal only by wiring a flag that defaults to the same
behaviour).

The cause is one pathway:

| pathway | scored | non-converged |
|---|---|---|
| R-HSA-3700989 TP53 | 232 (41% of all cases) | **122 (52.6%)** |
| every other pathway | 332 | **0 (0.0%)** |

TP53 alone supplies 41% of the scored cases and over half of them do not
converge. Non-converged solves return different values for a structurally
identical network, because `uuid4` node ids are minted afresh on every
catalog build and Dict iteration order over those ids sets the Gauss-Seidel
sweep order inside an SCC. So a rebuild of the same networks moves tens of
cases.

**Consequence for this feature:** the headline accuracy carries roughly ±39
cases of build-to-build noise. That is larger than any effect worth measuring
except the assembly clamp itself, so no A/B here can be trusted until it is
addressed. Restricting to cases where the solve converged gives a stable
comparison:

| | all 564 | converged 442 |
|---|---|---|
| DeltaSignal | 334 (59.2%) | **237 (53.6%)** |
| MP-BioPath | 407 (72.2%) | **301 (68.1%)** |

The gap is real on the stable subset — **14.5 points** — so the feature's
premise stands. But the measurement must be fixed before the fix is measured.
This adds FR7 and SC0 below.

## Established Findings *(evidence carried into this spec)*

Measured on the ten evaluation pathways, root-input perturbation, identical
case sets:

| arm | correct | accuracy | UP correct | DOWN correct |
|---|---|---|---|---|
| DeltaSignal, clamp ON (default) | 295-334/564 | 52-59% | 76-119/246 | 169-173/247 |
| DeltaSignal, clamp OFF | 368/564 | 65.3% | 165/246 | 175/247 |
| shortest-signed-path baseline | 393/564 | 69.7% | — | — |

Failure decomposition for expected-UP cases, clamp ON:

| | count |
|---|---|
| no directed path | 55 (22.4%) |
| path exists, readout at exactly baseline | 49 (19.9%) |
| wrong direction | 42 (17.1%) |
| attenuated into the no-change band | 24 (9.8%) |
| correct | 76 (30.9%) |

Two explanations were tested and **eliminated**:

- **Not reachability.** `no_directed_path` is 21.1% with all occurrences
  pinned and 22.3% with root inputs only. Restricted to reachable cases,
  accuracy still falls from 85.1% to 61.0%.
- **Not readout-local structure.** Untouched readouts are 30.0%
  multi-producer versus 26.1% for correct ones, and 5.0% versus 4.0%
  assembly — indistinguishable. `max` masks decreases, not increases.

The confirmed mechanism is `DS_ASSEMBLY_LIMITING`. It aggregates a complex's
subunits with `min`, and `min(elevated, baseline)` is exactly baseline, so
every complex on the path hard-clamps any increase while transmitting scarcity
perfectly. Turning it off recovers 73 cases and collapses the
untouched-readout count from 49 to 2.

The rule was adopted as a default on the strength of a +0.9pp curator win
measured when perturbations were pinned across a median of 17 nodes, some
adjacent to the readout. Under those conditions an increase barely had to
cross a complex, so the clamp cost almost nothing. It is not a bug so much as
a rule whose cost was invisible under the old measurement.
