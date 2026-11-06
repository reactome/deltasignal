# Research: the readout protocol in the engine

Everything below was read from the bench code, not guessed. File references are
to the repo at the commit this lands on.

## D1 — The direct step is the pin index, verbatim

`bench/benchmark_vs_mpbiopath.py::load_stid_to_uuids` builds stable_id → uuids
from `nodes.csv` exactly as `src/api/pins.jl::build_pin_index` does: indexed by
`diagram_entity_id` and every `member_leaves` stable id, `set_pool` rows
skipped, `(key, uuid)` deduplicated, insertion order kept; legacy
`stid_to_uuid_mapping.csv` only when `nodes.csv` is absent (the engine refuses
that case for pins, and does the same here).

**Decision**: the readout's direct step reads `PinIndex.stid_to_uuids`. No
second index; the specs/140 parity gate already held this map equal to the
bench's on every benchmark occurrence.

**Alternative rejected**: a separate loader — a second copy of the same rule is
the disease specs/140 treats.

## D2 — The resolution ladder, as the bench actually runs it

From `run_pathway` (the `key_output_uuids` block):

1. **direct**: first candidate stable id with a non-empty `stid_to_uuids`
   entry; groups = `[uuids]` (one group).
2. **proxy**: first candidate with entries in
   `entity_reaction_proxy_mapping.csv` (`load_entity_reaction_proxies`: every
   row, producing AND consuming roles, order kept, duplicates kept). Grouped by
   `"|".join(sorted(output_entities[u])) or u` in first-occurrence order, where
   `output_entities` (`load_output_entities`) maps each reaction uuid to the
   non-cofactor `diagram_entity_id`s of its `output`-edge targets
   (`nodes.csv` entity map is ALL rows — set_pool included — matching
   `PinIndex.diagram_entity`; cofactors from `cofactors.csv` `stable_id`).
3. **set_members**: members resolved direct-only, one `groups.append` per
   member in order, deduplicated across candidates (`seen_m`); empty appends
   filtered at the end (`[g for g in groups if g]`). `uuids` empty ⇒ basis
   none. Members come from Neo4j in the bench; the API takes them from the
   client, and the gate feeds the engine the bench's own member lists.

## D3 — Aggregation and the call

`entity_mean`: mean over non-empty groups of group means, values
`activities.get(u, BASELINE)` with `BASELINE = 0.01` (internal scale — the
same scale `/api/solve` returns `node_activities` in, kept 0–1 at
`server.jl` line ~973). `pred_ui = agg * CEILING` with `CEILING = 100.0` at
the default. `classify`: `< 0.85 → down`, `>= 1.15 → up`, else `normal`
(strict `<`, inclusive `>=`).

**Decision**: the endpoint takes `activities` as `/api/solve` returned them
(optional; omitted ⇒ resolution-only response) and reports `value` (internal),
`fold` (`value × 100`) and `call`. Stateless beyond the existing parse cache.

**Alternative rejected**: a server-side solve store (a new statefulness class
for no protocol gain).

## D4 — The abundance mode is the korkut protocol, not the ladder

From `bench/korkut_bench.py`:

- **Resolution** (`readout_nodes`): union of direct matches over ALL the
  readout's stable ids, `dict.fromkeys` order-preserving dedup, one flat group.
  No proxy or set-member steps.
- **Cut** (`abundance_copies`, version `abundance_a3` = what `abundance`
  names): close the resolved set over carry edges
  (`released_twin, set_member, specialization, variant_pool, variant_split`);
  drop any node that is the target of a made-from edge
  (`assembly, dissociation, released_twin, set_member, specialization,
  variant_pool, variant_split` — `recycle` excluded by amendment 3) whose
  source is in the carried set, or that an `output` edge produces from a
  reaction taking an `input` in the carried set. Order kept.
- **Statistic** (the scoring loop, line ~284): `free = kept not in pinned`;
  no resolved nodes ⇒ `no_copy`; none free ⇒ `trivial` (leaves scoring);
  else `ui = 100 × mean(act.get(u, 0.01) for free)`, classified at the same
  cutoffs.

**Decision**: `mode: "abundance"` reproduces all three; the request may name
`pinned` uuids, and the response carries `status`
(`scored` / `trivial` / `no_copy`) plus which resolved copies were pinned, so
a client can state "your perturbation pins every copy of this readout" instead
of inventing a value.

## D5 — Server wiring

- `ReadoutIndex` (proxies + cofactor set + which tables existed) built in
  `parse_handler` beside the pin index, from
  `entity_reaction_proxy_mapping.csv` and `cofactors.csv` next to the logic
  network, cached in a `READOUT_INDEX_CACHE` under the same lock.
- Uploads: `PARSE_TABLE_FIELDS` gains `"proxies" =>
  "entity_reaction_proxy_mapping.csv"`; `cofactors` and `nodes` parts already
  exist.
- Edge orientation: `LogicNetworkEdge.parent_uuid` is the bench's `source_id`,
  `child_uuid` the `target_id` (as `pin_root_indegree` reads them; the 140
  gate validated the convention end to end).
- Errors: `READOUTS_QUERY_ERROR`, `READOUTS_NO_TABLE_ERROR`, both on the
  `user_facing_error` allow-list. Missing proxy/cofactor tables do not 400 a
  query that never needed them; a query that needed the proxy step on a bundle
  without the table falls through (D2 order) and the response names the
  missing tables.

## D6 — The parity gate (FR-003/FR-004)

`bench/analysis/141/` in the 140 pattern:

- `dump_readouts.py`: for every (pathway, key_output) of the curator benchmark
  on the frozen build `20261026-1349_8b6527e`, call the bench's own loaders and
  resolution block (imported, not copied) and dump `{basis, groups}` — plus,
  for every korkut total readout, `readout_nodes` + `abundance_copies` kept
  lists and the member lists fed from Neo4j.
- `readouts_141.jl`: the engine side over the same cases (members injected from
  the dump).
- `compare_141.py`: field-by-field equality; any difference listed; non-zero
  exit. Result recorded below before the endpoint is called done.

## Parity result

Run 2026-11-06 on the frozen build `20261026-1349_8b6527e`
(`bench/analysis/141/bench_readouts.json` / `engine_readouts.json`,
`compare_141.py`):

```
curator: 983 (pathway, key_output) readouts compared; agreeing: 983 (direct 889, none 8, proxy 81, set_members 5)
korkut:  272 (pathway, total readout) cases compared; agreeing: 272
no differences — the gate passes
```

Every (pathway, key_output) readout of the curator benchmark resolves to the
same basis and the same curated-entity groups in the engine as in the bench's
own code, and every korkut total readout resolves to the same occurrence list
and the same abundance_a3 kept copies. The bench side was produced by
`resolve_key_output` (the function `run_pathway` itself now calls — extracted,
not copied; the bench suite passes unchanged, 408 passed / 3 skipped) and by
`korkut_bench.readout_nodes` / `abundance_copies`. Set members were fed to the
engine from the bench's Neo4j answer, as registered.

What this does NOT establish: member-list provenance is outside the gate by
design (the client supplies members in production).

## Live smoke (dev API, same build)

Parse → `/api/pins` → `/api/solve` → `/api/readouts` on the restarted dev
API: a proxy-resolved readout (ko 1268261 in R-HSA-5673001, RAF/MAP kinase
cascade), a set readout (ko 1268261 in R-HSA-5654736) and their pathways'
direct readouts all return the bench's groups, and every engine `value`
equals the bench's `entity_mean` over the same solve's activities to within
1e-12 (observed exactly equal). Abundance smoke: BAK in R-HSA-109606
(Intrinsic Apoptosis) resolves 4 occurrences and keeps exactly the bench's 1
abundance copy. The two cases today's UI shows nothing for (expanded species,
sets) resolve through the endpoint on real catalog pathways — SC-004.
