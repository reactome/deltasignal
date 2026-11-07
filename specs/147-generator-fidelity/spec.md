# Feature Specification: Generator Fidelity Audit — the Network vs Neo4j

**Feature Branch**: `specs/147-generator-fidelity`

**Created**: 2026-11-07

**Status**: Draft (analysis; registered before measurement)

**Input**: Adam, 2026-11-07: "do a very in-depth check with the logic network generator to see if you think it logically makes sense when comparing to what is in Neo4j. I feel like something has been going on with how the logic networks are being generated that is not correct. Looking for routes that should exist but don't would be a great start, and then making sure the pathways are being connected correctly."

The record's standing evidence agrees: half of all errors have no route in the scored network (specs/080); two-thirds of all validation readouts are untouched by the perturbation (specs/146); and the last specs/080 run found 392 curator `no_path` errors whose route exists in the pathway's own curated reactions in Neo4j but not in our network — measured several generator releases ago. Every large accuracy win on record was a network-construction fix.

## Phases (each produces a committed table before conclusions)

- **P1 — Lost routes, refreshed.** On the frozen record (build `20261026-1349_8b6527e`, solver `05dc867`): `failure_structure.py` → `case_oracle.py` (Reactome's own bipartite route per failed case, strict and lenient, same scope as the network) → `route_breaks.py` (the FIRST step of Reactome's route our network does not reach, typed). Deliverable: the refreshed count of severed-route errors and the break-class table, held-out and all.
- **P2 — Reaction-level structural diff.** Independent of errors: for every benchmark pathway, every curated ReactionLikeEvent in Neo4j (inputs, outputs, catalysts, regulators, within-pathway scope as the generator sees it) checked against the generated bundle — reactions absent entirely, participants dropped, roles mis-mapped, orientation flipped. Deliverable: per-defect-class counts with five traced examples each, distinguishing deliberate generator policy (cofactor handling, set expansion, variant caps — each matched to its recorded spec) from unexplained loss.
- **P3 — Pathway connection.** How sub-events are joined: `hasEvent` recursion into the pathway's own sub-pathways, `precedingEvent`/`hasEncapsulatedEvent` edges the generator uses or ignores at the pathway boundary, and whether P1's breaks concentrate at sub-pathway seams. Deliverable: the seam classes, counted, with traced cases.

## Rules

- Analysis only: no generator or solver change in this spec. Each defect class found becomes its own registered fix with a parity/benchmark gate (the specs/018 pattern).
- Trace cases, not aggregates: every claimed defect class carries at least five traced (pathway, reaction) examples, verified against Neo4j directly.
- Deliberate policies are not defects: anything explained by a recorded decision is listed separately with its spec reference.
- Neo4j queries bounded; the poetry venv python runs them.

## Success Criteria

- **SC-001**: A ranked table of generator defect classes by error mass (P1) and by structural mass (P2), each class traced, each separated from recorded policy.
- **SC-002**: A named, bounded candidate fix for the largest class, ready to register as its own spec.
- **SC-003**: Negative findings (classes that turn out to be policy or curation absence) recorded with the same completeness.
