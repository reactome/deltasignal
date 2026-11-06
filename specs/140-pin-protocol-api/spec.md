# Feature Specification: The API resolves a perturbation the way the benchmarks did

**Feature Branch**: `specs/140-pin-protocol-api`
**Created**: 2026-11-06
**Status**: Draft
**Input**: Adam, 2026-11-06: "continue trying to figure out how to make our system usable."

## The problem

Every accuracy figure this project can cite — 89.8% on the curator held-out
split, 643/849 on the MP-BioPath experimental axis, 10th of 75 on HPN-DREAM —
was measured under the **protocol of record** (Adam, 2026-09-25): perturb only
ROOT inputs that are, or contain, the gene; every downstream copy is
**computed**, never set. That protocol lives in one place:
`pin_occurrences` in `bench/benchmark_vs_mpbiopath.py`, Python the product
never runs.

The product runs something else. The Pathway Browser integration
(WebsiteAngular, on its main branch) pins **every node whose `reactome_id`
equals the clicked entity's stable id, wherever it sits in the network, at
confidence 1** (`matchingNodes` → `buildObservations` in
`projects/pathway-browser/src/app/deltasignal/`). That is the `all`-family
protocol retired by specs/023, whose numbers were recorded as "not
comparable": mid-pathway copies are *set* rather than computed. Concretely:

- Overexpressing a protein pins its produced, recycled and mid-pathway copies
  at the dose instead of letting the network answer whether they rise.
- A gene whose only starting form is a complex containing it (common — that
  is why the bench pins is-OR-CONTAINS) cannot be perturbed correctly by
  clicking the protein: none of its exact-id copies is a root, so the UI pins
  only downstream copies.
- A gene with no root form at all (a catalytic cycle regenerates it) needs the
  specs/024 `root_cycle` fallback; the UI has no notion of it.

So the answers users see are from a model configuration we never validated and
have specifically measured to be wrong. Making the system usable starts with
making the product give the measured model's answers.

## The rule

The engine gains a pin-resolution function and endpoint implementing the
protocol of record, so any client can ask: *"I want to perturb this entity —
which nodes do I pin?"*

**Resolution** (a port of `pin_occurrences`, scope `root_cycle`,
`DS_PIN_ROOT_IGNORE=depletion`, the defaults since specs/024 + specs/061):

1. **Occurrences.** A query is a list of stable ids (the clicked entity, or
   every form of a gene). Its occurrences are every node indexed by any of
   those ids, where a node is indexed by its `diagram_entity_id` and every
   stable id in its `member_leaves` (the generator's `nodes.csv`); `set_pool`
   nodes are never pin targets (specs/033).
2. **Roots.** Keep the occurrences with in-degree 0, where in-degree ignores
   `recycle` edges (specs/081), `context_bridge` edges (specs/056) and
   `depletion` edges (specs/058) — except that a node fed only by ignored
   depletion edges and with no out-edge is NOT a root (week review
   2026-10-15: a depletion-fed dead end pins nothing that moves).
3. **`root_cycle` fallback** (specs/024). If a query has occurrences but no
   root: candidate nodes are those with at least one out-edge whose every
   producer is downstream of them or reaches them only by a derived
   `dissociation`/`depletion` edge. Prefer candidates whose
   `diagram_entity_id` is one of the query's own ids; otherwise the
   candidates no other candidate reaches (or, if they all regenerate each
   other, all of them).
4. **Basis.** Each query reports how it resolved: `root`, `root_cycle`
   (the fallback fired), or `none` (no pinnable form — the honest answer the
   UI must surface instead of pinning mid-pathway copies).

**Endpoint** (additive; nothing existing changes):

`POST /api/pins` `{network_id, queries: [{id, stable_ids: [...]}, ...]}` →
for each query, `{id, pins: [{uuid, name}], basis, occurrences}` plus a
top-level `protocol` object naming the scope and root-ignore rule in force —
the run's protocol on the record, as every bench log prints it. A network
parsed without a node table (`nodes.csv`) cannot resolve the protocol and is a
`400` saying so, not a silent fall-back to exact-id matching.

## Requirements

- **FR-001**: Pure resolution functions in `src/` (no HTTP in them), taking
  the parsed network's edges and a node index; unit-tested with fixtures that
  compute (the fixture rule: the quantity a rule reasons about must be
  computed by the fixture, not pinned).
- **FR-002**: `POST /api/pins` as above. Unknown `network_id` → the existing
  re-parse `400`. Missing node table → `400` with a message naming
  `nodes.csv`. Malformed queries → `400`.
- **FR-003 (parity gate)**: On the frozen catalog
  (`20261026-1349_8b6527e`), for every benchmark (pathway, gene) the curator
  bench scores, the engine's pins EQUAL `pin_occurrences`'s pins (as sets)
  and the basis agrees (root / fallback / none), under default env. Measured
  by a dump-and-diff script committed with the feature; the diff must be
  empty. This is the gate: a port that is "close" is a second protocol, which
  is the disease being treated.
- **FR-004**: `docs/API.md` documents the endpoint (additive) and the client
  note steers the UI flow through it.
- **FR-005**: Tests live in an asserting suite listed in
  `test/asserting_suites.txt`; the CLAUDE.md test table row is added in the
  same PR (`scripts/check_doc_counts.py` enforces it).

## Non-goals

- No change to `/api/solve`, the solver, or any benchmark number.
- No readout-side protocol (entity-mean aggregation, the specs/134 abundance
  rule for total readouts) — that is the next feature of the usability
  programme, not this one.
- No gene-symbol resolution (the bench uses Neo4j for gene → stable ids; the
  API's clients already hold stable ids from the diagram).
- No UI change in this repo. The WebsiteAngular change (call `/api/pins`,
  show what was pinned and why, refuse honestly on `none`) is specified in
  that workspace once this lands.

## The usability programme (context, not scope)

This is the first of four gaps between what we validated and what a user
gets; recorded here so the sequence is on the record:

1. **Pins** (this spec): the product must perturb the way the benchmarks did.
2. **Readouts**: per-entity answers (specs/105 `mean_entity`) and the
   total-vs-forms distinction (specs/134) exist only bench-side; a UI user
   reading a single variant node's activity is reading a number we aggregate
   differently when we score it.
3. **Honesty surfacing**: `converged`/`residual` are served but not shown;
   railed readouts (specs/108: ≤0.01× or ≥99×) and bistable switches
   (specs/054 RAF) should be flagged in the UI, not presented as confident
   magnitudes.
4. **Scope statement**: what the model is validated FOR (direction within a
   curated pathway under a root perturbation; ranking readouts after a known
   inhibitor) and what it is NOT (unseen-drug direction, specs/121; per-cell
   effects, specs/121/130/139) belongs in the product's own words, not only
   in specs.
