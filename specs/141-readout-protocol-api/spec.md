# Feature Specification: The readout protocol in the engine — /api/readouts reads an entity the way the benchmarks did

**Feature Branch**: `specs/141-readout-protocol-api`

**Created**: 2026-11-06

**Status**: Draft

**Input**: Usability programme gap 2 (the programme is recorded in `specs/140-pin-protocol-api/spec.md`; Adam 2026-11-06: "continue trying to figure out how to make our system usable").

## The problem

Every citable DeltaSignal figure was produced by reading a readout in a specific way: resolve the curated entity to the network nodes that represent it (direct stable-id match, else the flux of its proxy reactions, else its set members), weigh each curated entity once rather than once per variant copy (entity-mean, specs/105, adopted at the freeze), and — for a total-protein question — cut the resolved nodes to the copies nothing makes from the protein itself (the abundance rule, specs/134, adopted on korkut/CPPA). That protocol lives only in bench Python (`bench/benchmark_vs_mpbiopath.py`, `bench/korkut_bench.py`).

The product reads a readout a way no figure was measured under: the Pathway Browser flat-means every node copy whose `reactome_id` equals the clicked entity, shows `change = activity − baseline` on the display scale rather than the fold the validated ±15% calls are made on, and shows **nothing at all** for the two cases the protocol handles by construction — a curated species expanded into variants (no node of its own; the bench reads its producing-reaction flux grouped by output entity) and a set readout (the bench aggregates member species, one group each).

This is the same shape of defect specs/140 fixed for pins, on the readout side. The fix follows the same pattern: the engine serves the protocol of record, a parity gate proves the port equals the bench on every benchmark case, and the UI then consumes the endpoint instead of reimplementing the rules.

## The protocol of record (what is being ported)

1. **Resolution** (bench `run_pathway`, specs/027/105): try the entity's stable id in the bundle's stable-id-to-uuid mapping — one group holding those uuids; if it resolves nothing, read the entity's proxy reactions from the bundle's entity–reaction proxy mapping, grouped by the sorted set of non-cofactor output entities each proxy produces (proxies making the same entities share a group; a proxy with no non-cofactor output is its own group); if that also resolves nothing and the entity is a set, resolve each member stable id direct-only, one group per member that resolves. The bench learns members from Neo4j; the engine has no Neo4j, so members arrive from the client (the Pathway Browser resolves them from the Reactome ContentService, the same source the bench's `neo4j_set_member_stids` reads).
2. **Aggregation** (specs/105 `entity_mean`): the mean over non-empty groups of each group's mean activity. A curated entity counts once, not once per variant copy.
3. **The value a user reads**: the fold over rest (`pred_ui` = aggregated internal activity × the ceiling, exactly 100 at the default baseline), and the validated direction call at 0.85 / 1.15 on that fold.
4. **Abundance mode** (specs/134, `abundance_a3`, as `korkut_bench.py` implements it): a total-protein readout resolves differently — the union of direct matches over **all** its stable ids (the korkut `readout_nodes` rule; no proxy or set-member steps), one flat group, order kept and deduplicated. The cut then keeps only copies nothing makes from the readout's own copies: first close the resolved set over the carry edges (`released_twin`, `set_member`, `specialization`, `variant_pool`, `variant_split` — a set or pool of the protein's copies IS the protein, amendment 2); then drop any copy with a made-from edge in from the carried set (`assembly`, `dissociation`, `released_twin`, `set_member`, `specialization`, `variant_pool`, `variant_split`; `recycle` excluded by amendment 3) or produced by a reaction that takes an input in the carried set. Aggregation is the korkut statistic: the mean over **unpinned** kept copies; a readout whose every kept copy is pinned leaves scoring (the endpoint exposes per-copy pinned flags and lets the client decide, per the Assumptions). This is per query — "what happens to total X" and "what happens to this form of X" are different readouts under different measured protocols.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Resolve and aggregate a readout the measured way (Priority: P1)

A client (the Pathway Browser, or any API consumer) solved a perturbation and wants the predicted level of a curated entity. It asks the engine, which resolves the entity exactly as the benchmarks did and returns the curated-entity groups, the entity-mean value over a supplied solve's activities, the fold, and the direction call at the validated cutoffs.

**Why this priority**: this is the gap — the product's displayed values are produced by a rule no figure was measured under, and expanded species and sets display nothing.

**Independent Test**: upload a fixture network with a direct two-copy readout, solve, query the endpoint, and check groups/value/fold/call against hand-computed entity-mean.

**Acceptance Scenarios**:

1. **Given** a network where entity E maps directly to two variant copies, **When** the client queries E with a solve's activities, **Then** the response has basis `direct`, one group of two uuids, and value = the mean of the two activities.
2. **Given** an entity with no node of its own but proxy reactions making two distinct non-cofactor output-entity sets, **When** queried, **Then** basis is `proxy` and there are two groups (one per output-entity set), and the value is the mean of the two group means — not the flat mean over proxies.
3. **Given** a set readout whose two members resolve to 1 and 3 copies, **When** queried with the member stable ids, **Then** basis is `set_members`, there are two groups, and the value weighs each member once (entity-mean ≠ flat mean distinguishes the fixture).
4. **Given** an entity that resolves nowhere, **Then** basis is `none` and no value is invented.

---

### User Story 2 - Total-protein (abundance) readouts (Priority: P2)

A client asks what happens to **total** protein X. The engine cuts X's resolved nodes to its abundance copies under the adopted `abundance_a3` rule before aggregating.

**Why this priority**: the korkut/CPPA figures for total readouts were measured only under this rule (specs/134: accuracy 0.674 → 0.796); a UI "total X" answer made from all forms would be the same class of unmeasured answer as the old pins.

**Independent Test**: a fixture where a complex is assembled from the protein's starting copy and a recycle edge returns a released copy; the abundance query must keep the starting copy (recycle is not "made from") and drop the complex.

**Acceptance Scenarios**:

1. **Given** protein P with a starting copy, a complex assembled from it, and a released copy recycled back to the start, **When** queried with mode `abundance`, **Then** only the starting copy remains in the group (the complex is made from P by `assembly`; the recycle edge does not make the starting copy "made"; and the value is the mean over the unpinned kept copies).
2. **Given** the same network, **When** queried with the default mode, **Then** all resolved copies remain (the modes differ on this fixture).
3. **Given** a readout whose every abundance copy is pinned by the perturbation itself, **Then** the response says so rather than scoring the pin back (the specs/134 own-target rule is the client's concern for scoring, but the response must expose which group members are pinned so a client can apply it) — resolved as: the response flags per-uuid pinned status when the request names pinned uuids.

---

### User Story 3 - Honest degradation and errors (Priority: P2)

A bundle missing the proxy table or node table cannot run the full protocol. The endpoint must say what it could not do, never approximate silently; malformed requests get the existing allow-listed error treatment.

**Acceptance Scenarios**:

1. **Given** a network uploaded without the proxy mapping, **When** an entity would need proxy resolution, **Then** resolution falls through to `set_members`/`none` and the response names the missing table.
2. **Given** an unknown network id, a query without `id`, `stable_ids` of the wrong type, or an unknown `mode`, **Then** the client gets a 400 whose message is on the user-facing allow-list.
3. **Given** `abundance` mode on a bundle without the tables the cut needs, **Then** the response refuses the mode for that query (basis reported, mode error stated) instead of returning a forms answer labelled abundance.

---

### User Story 4 - The parity gate (Priority: P1)

Before the endpoint is trusted, its resolution must be shown equal to the bench's on every benchmark readout, the way specs/140's gate did for pins.

**Independent Test**: run the gate; it fails on any difference.

**Acceptance Scenarios**:

1. **Given** the frozen build `20261026-1349_8b6527e`, **When** every (pathway, key_output) readout of the curator benchmark is resolved by the bench's own functions and by the engine, **Then** basis and groups agree on every case, and the comparison lists any difference and exits non-zero on one.
2. **Given** the korkut benchmark's total readouts on the frozen build, **When** the abundance cut is computed by `korkut_bench.abundance_copies` and by the engine, **Then** the kept uuid sets agree on every readout.
3. **Given** matching resolutions and a shared activities vector, **Then** engine entity-mean values equal the bench's `entity_mean` to floating-point equality on a sampled pathway's full solve.

### Edge Cases

- A readout with several candidate stable ids (the bench tries a resolved stId then the historical `R-HSA-{dbid}` guess): first candidate that resolves wins, across the whole resolution ladder per candidate (direct for all candidates, then proxy for all candidates — the bench's order) — the gate decides this question against the bench's actual behaviour, not this sentence.
- Empty groups (a member that resolves nothing) are skipped by entity-mean, as the bench skips them; a query where all groups are empty is basis `none`.
- A uuid in `activities` that is not in the network: ignored for that network's nodes; a resolved uuid absent from `activities` reads baseline, as the bench reads `activities.get(u, BASELINE)`.
- Duplicate stable ids across `stable_ids` and `members`: deduplicated the way the bench deduplicates member expansion (`seen_m`).
- `activities` omitted: resolution-only response (groups + basis), no value/fold/call fields.
- Cofactor outputs: a proxy whose outputs are all cofactors groups under its own uuid, as the bench's `"|".join(...) or u` does.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The engine MUST resolve a readout query by the three-step ladder of the protocol of record (direct → proxy grouped by non-cofactor output entities → client-supplied set members direct-only), returning the basis and the curated-entity groups.
- **FR-002**: The engine MUST aggregate with the specs/105 statistic exactly: mean over non-empty groups of group means; the fold and the 0.85/1.15 direction call MUST be computed from that value as the bench computes `pred_ui` and `classify`.
- **FR-003** (the gate): On the frozen build `20261026-1349_8b6527e`, the engine's resolution MUST equal the bench's on **every** (pathway, key_output) readout of the curator benchmark, with the bench side produced by the bench's own functions (`load_stid_to_uuids`, `load_entity_reaction_proxies`, `load_output_entities`, the `run_pathway` resolution block, `entity_mean`), never a reimplementation; `bench/analysis/141/` holds the dump/compare harness in the specs/140 pattern, any difference fails it, and the result is recorded in `research.md` before the endpoint is called done.
- **FR-004**: The abundance mode MUST reproduce the korkut total-readout protocol exactly: union-of-direct-matches resolution (`korkut_bench.readout_nodes`), the `abundance_a3` cut (`korkut_bench.ABUNDANCE_RULES` / `abundance_copies`), and mean-over-unpinned aggregation — gated on every korkut total readout of the frozen build with the bench side produced by those functions.
- **FR-005**: Everything is additive: no existing endpoint, response field, solver behaviour, or bench number changes; a bundle missing a table the protocol needs degrades honestly (the response names the missing table; abundance refuses rather than approximates).
- **FR-006**: Unit fixtures MUST compute the thing each rule reasons about (never pin the intermediate): the four US1 scenarios, the US2 abundance/forms contrast, missing-table degradation, and malformed-query 400s on the allow-list.

### Key Entities

- **Readout query**: an id, ordered candidate stable ids, optional member stable ids, a mode (`forms` default | `abundance`).
- **Resolution**: a basis (`direct` | `proxy` | `set_members` | `none`) plus groups — lists of node uuids, one group per curated entity.
- **Readout value**: entity-mean internal activity over a supplied solve's activities; its fold over rest; the direction call at the validated cutoffs.
- **Readout index**: the per-network tables resolution reads (stable-id map, proxy map, output entities, cofactors, abundance edge walk), built at parse time and cached beside the network, like the pin index.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The parity gates pass with zero differences on every curator-benchmark readout (resolution) and every korkut total readout (abundance) on the frozen build.
- **SC-002**: A client can obtain, for any curated entity in a parsed pathway, either the measured protocol's value/fold/call or an explicit statement of why not (`none`, missing table, refused mode) — no silent approximations exist in any response path.
- **SC-003**: The existing test suites pass unchanged; the new suite's fixtures distinguish entity-mean from flat mean and abundance from forms (each fixture fails under the wrong rule).
- **SC-004**: The two cases today's UI shows nothing for (expanded species, set readouts) return resolved groups through the endpoint on a real catalog pathway (smoke-tested on the dev API).

## Assumptions

- The client supplies set members (from the Reactome ContentService); the engine never queries Neo4j. The gate accounts for this by feeding the engine the same member lists the bench resolved from Neo4j, so member-list provenance is outside the gate.
- The GTPase readout remap (`DS_KO_REMAP`, bench-only recuration shim) and the release-skew exclusion bookkeeping are bench scoring concerns, not part of the protocol a product reads an entity by; they stay out of the endpoint. If the gate shows remapped readouts differing, those cases are excluded from the gate with the exclusion recorded, as specs/140's gate recorded `none` cases.
- The own-entity pin interaction of specs/134 (a drug's own fully-pinned target leaves scoring) is served as information (per-uuid pinned flags when the request names pins), with the leave-scoring decision made by the client, because "leaves scoring" is a benchmark concept while a UI must still show the user something about their own perturbation target.
- The solve's activities travel in the request rather than being stored server-side, keeping the API stateless beyond the existing parse cache.
- Non-goals: no change to how the solver computes node activities; no change to `/api/solve`; no UI work in this repo (WebsiteAngular consumption is its own spec once this is served); no new catalog build; no Neo4j dependency in the engine; no change to bench scoring.
