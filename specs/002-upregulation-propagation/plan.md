# Implementation Plan: Long-range upregulation propagation

**Branch**: `002-upregulation-propagation` | **Date**: 2026-09-10 | **Spec**: [spec.md](./spec.md)

## Summary

DeltaSignal scores 334/564 against MP-BioPath's 407/564 on root-input
perturbation cases, and loses to a model-free graph traversal (393/564). The
gap is almost entirely upregulation: an increase applied at a root input
frequently arrives at the readout as *exactly* baseline.

The plan has two stages, and the order is not negotiable. The headline figure
currently moves by 39 cases on a catalog rebuild, which is larger than most
effects worth measuring, so **the instrument is fixed before the model**.

## Technical Context

**Language**: Julia 1.10 (solver, API), Python 3.12 (benchmark harness)
**Primary surface**: `src/core/reaction_model.jl` — per-activator aggregation
**Evaluation**: `bench/benchmark_mpbiopath_cases.py`, ten MP-BioPath
experimental pathways, Reactome Release97, harness defaults (0.85/1.15)
**Comparators**: MP-BioPath and curator per-case published predictions;
`shortest_signed_path` structural baseline computed on identical cases

**Known cause (measured on one catalog build)**: `DS_ASSEMBLY_LIMITING`
aggregates a complex's subunits with `min`. `min(elevated, baseline)` is
exactly baseline, so every complex on a path hard-clamps an increase while
transmitting scarcity perfectly. Disabling it recovers **+21 cases (+42 on
upregulation, −3 on downregulation)**. An earlier +73 figure used a baseline
that has since been withdrawn.

**Scale check**: the clamp accounts for 29% of the 73-case gap to MP-BioPath.
It is the largest identified factor, not the majority of the problem. Stage 2
must not stop when it is fixed.

**Known measurement hazard**: `uuid4` node ids are reminted on every catalog
build; Dict order over them sets Gauss-Seidel sweep order inside an SCC; so
non-converged solves differ between structurally identical networks. TP53 is
41% of scored cases and 52.6% non-converged. Every other pathway is 0%.

**NEEDS CLARIFICATION**: none. The two open design questions (which softer
limiting rule; whether to fix TP53's SCC or report converged-only) are
resolved in Phase 0 by measurement rather than by asking.

## Constitution Check

| Principle | Status | Note |
|---|---|---|
| I. Processing is DeltaSignal's job | **PASS, and this is the case in point** | The clamp is a solver-side rule. The generator is representing complexes as curators define them; the propagator's handling of that representation is what fails. Fix belongs here, not upstream. |
| II. Measure, then claim | **PASS** (gate withdrawn — see research.md R1) | The ±39 alarm was wrong: two independent builds agree to ±1 case and the converged subset is bit-identical. The gate rested on a comparison against a catalog destroyed in a crash; that figure is withdrawn. A standing rule replaces it: an A/B is only valid across arms sharing one catalog build. |
| III. Negative results are results | PASS | Each candidate rule records its numbers in `research.md` whether it helps or not. |
| IV. New behaviour ships default-OFF | **JUSTIFIED DEPARTURE** | See Complexity Tracking. |
| V. Honest solver reporting | **PASS, and directly implicated** | 122 non-converged TP53 solves are currently reported and scored as if equivalent to converged ones. Stage 1 makes that visible. |
| VI. API boundary is a trust boundary | N/A | No API surface changes. |

**Gate verdict**: SC0 is satisfied (research.md R1). Stage 1 is complete.
Stage 2 may proceed, under the standing rule that comparisons share a catalog
build.

## Stage 1 — Make the measurement trustworthy (SC0, FR7)

Nothing else can be evaluated until a rebuilt catalog reproduces the headline.

1. **Quantify the instability directly.** Rebuild the same ten pathways twice
   with the seed pinned, run the same arm on both, and report the case-level
   delta and its per-pathway attribution. This turns "39 cases" from an
   inference about two runs into a measured property with a distribution.
2. **Establish why TP53 alone fails to converge.** Prior work characterises
   its giant SCC as Type II catalytic recycling — hundreds of virtual-reaction
   instances of one or two reaction stIds — rather than genuine feedback. If
   that holds, it is a network-structure problem, and the honest fix is
   upstream or in SCC handling, not in loosening the tolerance.
3. **Decide the primary figure.** Either TP53 converges, or the primary
   comparison is converged-only with the non-converged remainder reported
   separately. Converged-only is already stable: DeltaSignal 237/442 vs
   MP-BioPath 301/442.

**Exit criterion**: two independent catalog builds of an unchanged arm agree
to within a handful of cases, and the figure used for SC1–SC6 is named.

## Stage 2 — Restore long-range upregulation (SC1–SC6)

Only once Stage 1 exits.

1. **Characterise what the clamp should be.** `min` encodes something true —
   a complex cannot exceed its scarcest subunit in *absolute abundance*. It
   does not follow that a *fold-change* should be clamped to exactly zero.
   The candidate family is a soft minimum that preserves FR3 (scarcity still
   propagates; a knocked-out subunit still collapses the complex) while
   admitting FR2 (an elevated subunit moves the complex off baseline).
2. **A/B each candidate** on identical cases against MP-BioPath and the
   structural baseline, reporting per-pathway attribution and both-arms-
   converged counts per FR5.
3. **Check the original win survives.** The clamp was adopted for a +0.9pp
   curator gain (DSB +185). Any replacement must be measured against the
   curator ground truth too, not only the experimental set, or we trade a
   known win for an unknown one.
4. **Look past the clamp.** Clamp-off is 368/564 and MP-BioPath is 407 — so
   the clamp explains most of the gap but not all of it. The residual is
   Stage 2's second half: 42 wrong-direction and 24 attenuated-into-band
   cases remain unexplained.

## Project Structure

### Documentation (this feature)

```
specs/002-upregulation-propagation/
├── spec.md              # written
├── plan.md              # this file
├── research.md          # Phase 0 — candidate rules and their measurements
├── quickstart.md        # Phase 0 — how to reproduce every number here
└── checklists/
    └── requirements.md  # written, passing
```

No `data-model.md` or `contracts/`: this feature changes an internal
aggregation rule and adds no entities and no external interface.

### Source Code

```
src/core/reaction_model.jl          # assembly aggregation — the change site
bench/benchmark_mpbiopath_cases.py  # evaluation, unchanged by this feature
```

## Complexity Tracking

| Departure | Why needed | Simpler alternative rejected because |
|---|---|---|
| Principle IV says new behaviour ships default-OFF. A replacement limiting rule cannot ship OFF and also fix anything — the default path is what is broken. | The convention exists for modelling choices of *unknown direction*. Here the current default is measurably wrong for 73 cases, so leaving it on by default is itself the risk. | Shipping the new rule behind an opt-in flag would leave every default run with the clamp, i.e. would not close the gap. The mitigation is that the old behaviour stays reachable via the existing `DS_ASSEMBLY_LIMITING` switch, so the previous default is one env var away and the curator A/B can still be run both ways. |
