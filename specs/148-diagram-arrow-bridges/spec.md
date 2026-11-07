# Feature Specification: Diagram-Arrow Bridges — the Curators' Drawn Connections, Served as the Fall-Back Class

**Feature Branch**: `specs/148-diagram-arrow-bridges`

**Created**: 2026-11-07

**Status**: Draft (pre-registered before any measurement)

**Input**: Adam, 2026-11-07, on specs/147's P4 finding ("feel free to do so"): the benchmark diagrams carry 115 curator-drawn `pathway→entity` arrows (release 97; 113 on 98) that exist in no event-graph relation — an upstream pathway drawn as a box feeding a NAMED entity of this pathway ("Signaling by EGFR → RAS GEFs"). These are the entity-anchored cross-pathway joins specs/132/133 had to approximate: 132's tracing showed the entity anchor separates right bridges from wrong ones, and 133's transitive rule passed its balanced-accuracy gate but was declined for silently mixing fall-back answers with in-pathway ones. specs/146's machinery removes that objection: fall-back answers can be served as a separate, labelled class. The error mass this attacks: two-thirds of validation readouts are untouched within the curated pathway, wrong 70% of the time against experiments that expected a change.

## The rule (fixed before measurement)

For each centre pathway P, each release-97 arrow (Q → e into P) with Q in the catalog adds Q as a whole, unreduced source (the specs/133 amendment-2 discipline, `--whole`, with `context_origin.csv`; **amendment 1, 2026-11-07**: the arm runs `DS_CONTEXT_BRIDGE_ORDER=first` — the registered `last` fails the bench's carrier-twin capacity guard, which counts possible twins composite-wide while `last` correctly drops cross-pathway rule edges, a census/solver mismatch no 133 arm happened to trip; under `first` the cross-pathway rule edges are kept, counted and reported, and the F-discipline scoring bounds their exposure to the untouched class), bridged **feed-in only** at the drawn entity: from Q's produced copies of e (or of e's Neo4j member/component closure — the drawn target is often an overview set the per-pathway build excludes, e.g. RAS GEFs) to P's root copies of the same closure. No identity joins anywhere else; no bridge leaves P; arrows whose Q is outside the catalog are counted and skipped (reported).

**Scoring rule (the served semantics, the specs/133 F discipline):** P's own answer stands wherever P alone routes the readout; the bridged answer is used ONLY for readouts P alone leaves untouched — exactly the class the product serves as `confidence_stratum: "untouched"` today. The arm never changes an in-pathway answer by construction.

## Hypotheses

- **H1 (the fall-back class carries information)**: on korkut scored readouts whose frozen-control fold is exactly 1.0 (the untouched class), the bridged answers beat always-no-change on balanced accuracy, by-drug-set bootstrap excluding 0.
- **H2 (the curator axes tolerate it)**: under the F discipline, curator held-out and experimental case accuracies are unchanged except through newly-answered untouched readouts, and those move net-positive on each axis (fixed > broken among cases whose control answer was an untouched NORMAL).

## Protocol

- Catalog: the frozen build `20261026-1349_8b6527e` as base; arrows from the release-97 census (regenerate with `specs/147-generator-fidelity/p4/census.py ~/reactome-diagrams/97`); dbIds mapped to stIds in Neo4j; member closures via the oracle's `fetch_composition` machinery.
- Build: `bench/context_build.py` gains an `--arrows` planner (entity-anchored, whole sources, feed-in only, origin table). BUILD.json records the arrow list and skip counts.
- Measure: the curator bench (both splits), experimental, and korkut (raw and balanced, both thresholds) on the arm vs the frozen control, all through `run_arm.sh`-grade provenance; the F discipline applied in scoring (control answer where control routes).
- Traces before claims: which arrows fire per fixed/broken case; the PIP3 switch arrows specifically (the specs/133 failure mode) — if PIP3's rails carry through its arrows, that is reported case-level.
- korkut is FIT; the curator held-out is SELECT: development, honestly measured. Adoption (serving the class in the product) is Adam's on the numbers.

## Success Criteria

- **SC-001**: H1 and H2 as registered, or the honest negatives with per-arrow traces.
- **SC-002**: The skip census (arrows whose source pathway the catalog lacks) recorded — it sizes the whole-Reactome version.
- **SC-003**: Nothing about in-pathway answers changes under the F discipline (checked case-for-case).

## Non-goals

- No event-graph bridges (132/133 measured those); no identity joins; no UI work until the measurement stands; no new generator run (release-97 consistency throughout); no change to the frozen headline numbers — this creates a new, separately-labelled answer class, it does not touch the existing one.
