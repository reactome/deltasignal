# DeltaSignal HTTP API

This is the stable contract between the DeltaSignal engine (this repo) and any
UI/client — notably the Angular frontend that lives in the WebsiteAngular
workspace. The seam between the two repos is this API: the engine exposes it,
the UI consumes it, and neither shares code with the other.

Server: `src/api/server.jl`. Run locally via `docker compose -f
docker-compose.dev.yml up julia-api` (binds `0.0.0.0:8080` in the container,
published on the host as `localhost:8080`).

## Conventions

### Activity scale — read this first
The scale is **asymmetric by design** (it matches what every benchmark client
already does; do not "fix" it without changing those clients too):

- **Input** (`observations` you send to `/api/solve`): **0–100 UI scale.**
  `0` = no activity, `1` = normal baseline (1×), `100` = 100× normal.
- **Output** (`node_activities` returned by `/api/solve`): **0–1 internal
  scale.** Multiply by `100` to display on the same 0–100 UI scale.

```
display_value = node_activities[uuid] * 100.0
```

### Error shape
Every error returns a sanitized JSON body (raw Julia exceptions are never
leaked):
```json
{ "status": "error", "message": "human-readable, safe message" }
```
Status codes: `400` (bad/missing/malformed input), `404` (unknown route),
`500` (internal). Server-side the original exception is logged with a backtrace.

### CORS
Responses carry CORS headers on every request; `OPTIONS` preflight returns `200`.
- `Access-Control-Allow-Origin`: `DS_CORS_ORIGIN` env var, default `*`.
- `Access-Control-Allow-Methods`: `GET, POST, PUT, DELETE, OPTIONS`
- `Access-Control-Allow-Headers`: `Content-Type, Authorization`

For local Angular dev, either proxy `/api` → `http://localhost:8080` from the
Angular dev server (recommended — no CORS involved), or rely on the default
`*` origin for direct browser calls.

## Endpoints

### `GET /api/health`
Liveness probe.
```json
{ "status": "ok", "timestamp": "1.7e9", "version": "0.1.0",
  "service": "deltasignal-api", "julia_version": "1.10.10" }
```

### `GET /api/pathways`
The pre-generated Reactome pathway catalog (from the logic-network-generator),
deduplicated and pretty-named. This is the primary "pick a pathway" list.
```json
[ { "id": "Cell_Cycle_Checkpoints_R-HSA-69620",
    "stable_id": "R-HSA-69620",
    "name": "Cell Cycle Checkpoints" }, ... ]
```
`id` is the opaque catalog key to pass to `/api/parse` as `pathway_id`.

### `POST /api/parse`
Parse a logic network into nodes + edges. Three input modes (checked in order):

1. **Catalog pathway** (primary) — JSON body:
   ```json
   { "pathway_id": "Cell_Cycle_Checkpoints_R-HSA-69620" }
   ```
2. **File upload** (escape hatch) — `multipart/form-data` with parts
   `logic_network` and `uuid_mapping` (required), `set_mappings` (optional),
   and `boundary_edges` (optional, specs/044): a generator bundle's
   `boundary_edges.csv`, the edges it derives at root inputs and terminal
   outputs. Send it with a bundle from a split generator, or the network is the
   curated one alone. A catalog `pathway_id` reads it automatically.
   The bundle's tables can come as optional parts too, named `drugs`,
   `pathogens`, `carrier_loops`, `pools`, `pool_transitions`, `pool_carriers`,
   `containment`, `cofactors` and `nodes` (each the bundle's `<name>.csv`).
   Without them the solve runs, but the rules that read them are inert and
   say so (`drug_rule: "inert: no drug table"` and so on): send them to get
   the same model a catalog `pathway_id` gets.
   A request that declares `multipart/form-data` but is missing either required
   part is a **`400`** — it is not silently treated as mode 3.
3. **No body** — falls back to the bundled sample network. This applies only to
   a genuinely empty body. A request that sends a body but no usable
   `pathway_id` gets a **`400`**, because it asked for something specific and
   must not receive unrelated demo data with a `200`.

Upload is actually checked before the body is read, so a multipart request
never reaches mode 1; the ordering above matters only when more than one mode
could apply, which the `400`s now prevent.

Response:
```json
{
  "status": "success",
  "message": "Network parsed successfully",
  "network_id": "pw:Cell_Cycle_Checkpoints_R-HSA-69620",
  "nodes": [
    { "uuid": "...", "name": "TP53", "reactome_id": "R-HSA-...",
      "entity_type": "protein", "baseline": 0.01, "set_id": null } ],
  "edges": [
    { "parent_uuid": "...", "child_uuid": "...", "is_and": false,
      "is_positive": true, "stoichiometry": 1.0, "edge_type": "output" } ],
  "pathways": [ { "id": "...", "name": "...", "members": ["uuid", ...] } ],
  "containment": { "<stable_id>": ["<contained stable_id>", ...] },
  "cofactor_stids": ["R-ALL-113592", ...]
}
```
`containment` (restricted to stable ids in this network) and `cofactor_stids`
are what the solver needs to run the default model; send them back with an
inline network. `drug_stids` (specs/032, additive) lists the bundle's
drug-derived stable ids, or is `null` when the bundle has no `drugs.csv`; send
it back too, or `DS_DRUG_MODE=inert` reports `"inert: no drug table"`.
`pathogen_stids` (specs/048, additive) is the same for the bundle's
`pathogens.csv`. Without it, the default `DS_PATHOGEN_MODE=inert` reports
`"inert: no pathogen table"` and holds nothing.
`pools` (specs/039, additive) is the interconversion-pool table, or `null`
when the bundle has no `pools.csv`. Send it back with an inline network, or the
default `DS_CYCLE_MODE=balance` reports `"balance: no pool table"` and loops
iterate instead (until 2026-10-02 it was not returned at all, so every inline
solve silently ran without it).
`carrier_loops` (specs/049, additive) is the carrier-loop table, a list of
`{free_uuid, loop_source_uuid, closed}`, or `null` when the bundle has no
`carrier_loops.csv`. Send it back with an inline network, or the default
`DS_CARRIER_LOOP_MODE=closed_twin` (and every other mode except `off`) reports
`"<mode>: no carrier table"` and solves as `off`.
`node_leaves` (specs/052, additive) maps each node uuid to `{kind, leaves}`
(the generator's `nodes.csv` `node_kind` and `member_leaves`). It is `null`
when the bundle has none, and also whenever the server runs the default
`DS_ZERO_SUPPLY_MODE=off`, because the table is large and only that rule reads it. Send it back with an inline network, or
`DS_ZERO_SUPPLY_MODE=hold` reports `"hold: no node table"` and solves as `off`.

Node display names are enriched from the Reactome ContentService when the
generator only provided stable ids (degrades gracefully if that service is
down — names fall back to the stable id).

**`network_id`** is the key to the parse-once/solve-many flow: hand it to
`/api/solve` instead of re-shipping the whole network on every solve. It is an
in-memory, per-process cache — it does **not** survive a server restart, so a
client must be prepared to re-`parse` on a cache miss. A solve naming a
`network_id` the server does not hold returns **`400`** with a message saying
to re-`parse`; it does **not** fall through to the inline `network` or the
sample, which previously returned a confident `200` for a different network.

### `POST /api/solve`
Solve the steady state under a set of perturbations. JSON body:
```json
{
  "network_id": "pw:Cell_Cycle_Checkpoints_R-HSA-69620",
  "observations": { "<node_uuid>": [80.0, 0.9] }
}
```
- **Network source** (in priority order): `network_id` (from a prior parse) →
  inline `network` (same shape as the parse response's nodes/edges/pathways,
  plus its `containment` and `cofactor_stids`; without `containment` the
  self-inhibitor rule is inert and the response says so) →
  bundled sample. A `network_id` that is present but unknown is a `400`, not a
  fall-through: naming a network the server does not have is an error, whereas
  naming none at all still selects the sample.
- **`observations`**: object keyed by node `uuid`; each value is
  `[activity, confidence]` where `activity` is **0–100** and `confidence` is
  `0–1`. Confidence `> 0` pins the node as a hard constraint; `0` ignores it.
  A uuid repeated in the object with **different** values is a `400`
  (two contradictory measurements are a conflict); an exact repeat is fine.

Response:
```json
{
  "status": "success",
  "message": "Steady-state solved successfully",
  "node_activities": { "<uuid>": 0.30008 },
  "influence_scores": { "<uuid>": 6.87 },
  "converged": true,
  "iterations": 13,
  "solve_time": 0.016
}
```
- **`node_activities`**: `uuid → activity` on the **0–1 internal scale**
  (`× 100` for display — see Conventions).
- **`self_inhibitors`** (additive, specs/022): the number of inhibitor slots FLAGGED as self-contained (the inhibitor contains one of its own reaction's inputs AND is reachable from it in the network). A flagged slot is damped only when the inhibitor moves with that input. On by default (`DS_SELF_INHIBITOR_WEIGHT=0.1`; `off` disables it).
- **`drug_rule`** / **`drugs_held`** (additive, specs/032): `"inert"`
  (default since 2026-10-13, specs/055), `"propagate"`, or `"inert: no drug table"`, and how many drug-derived
  nodes this solve held at baseline. Under `inert` a drug participates at fold
  1.0 and never carries a perturbation; an explicit observation still wins.
- **`pathogen_rule`** / **`pathogens_held`** (additive, specs/048): `"inert"`
  (default since 2026-10-10), `"propagate"`, or `"inert: no pathogen table"`,
  and how many pathogen-derived nodes this solve held at baseline. A node is
  held as a drug is under `inert`.
  The built-in sample network has no pathogen table, so a bodyless
  `/api/parse` → `/api/solve` reports `"inert: no pathogen table"`; nothing
  is wrong with that solve.
- **`zero_supply_rule`** / **`zero_supply_leaves`** / **`zero_supply_held`**
  (additive, specs/052): `"off"` (default), `"hold"` or `"hold: no node
  table"`; how many proteins pinned to 0 were left with no unpinned source;
  and how many nodes built from them this solve held at 0.
- **`self_inhibitor_rule`** (additive): `"on"`, `"off"`, or `"inert: no containment table"`. The last means the solve did NOT run the default model, because the network came without a containment table. `/api/parse` returns `containment`; send it back with a POSTed network to keep the rule active.
- **`pool_supply_reach`**, **`cycle_supply_regulated_kept`** (additive,
  specs/071): which edges the pool supply test follows (`"all"`, the default,
  or `"material"`), and how many pool supply producers the solve kept because
  the pool reaches them only through catalyst or regulator edges.
- **`cycle_carrier_release`**, **`cycle_carrier_release_added`** (additive,
  specs/072): which producers of a carrier's free form count as its release
  steps (`"listed"`, the default, or `"reaction"`: every copy of a listed
  release reaction), and how many copies `"reaction"` added.
- **`nonconv_fallback`**, **`nonconv_components`**, **`nonconv_improved`**,
  **`nonconv_rescued`**, **`nonconv_skipped_pools`** (additive, specs/066):
  the fallback for a cyclic component still above tolerance at the iteration
  cap (`"lm"`, the default, or `"off"`); how many such components the solve
  had, how many the fallback improved, how many it brought below tolerance,
  and how many it skipped because they hold pool nodes.
- **`nonconv_pools`**, **`nonconv_pool_components`**,
  **`nonconv_pool_improved`**, **`nonconv_pool_rescued`** (additive,
  specs/073): what the fallback does with such a component when it holds pool
  nodes (`"skip"`, the default, or `"joint"`: re-solved with the pool update
  as part of the least-squares problem); how many it re-solved, improved and
  brought below tolerance.
- **`cycle_rule`**, **`cycle_pools_solved`**, **`cycle_pools_multistep`**,
  **`cycle_carriers`**, **`cycle_carriers_held`** (additive, specs/039):
  `"balance"` (default), `"off"`, `"balance: no pool table"`, or
  `"balance: inert under DS_SCC_METHOD=…"` / `"…DS_SCC_SOLVE=0"` (a solver
  method that never solves pools, so they were not applied); how many
  interconversion pools this solve handled as conserved pools, how many of
  them multi-step; the carrier rule (`"on"`/`"off"`) and how many enzyme free
  forms it held. `balance` is the default (specs/039); `off` restores
  iterated loops.
- **`cycle_shared_nodes`** (additive, specs/049): intermediates and step
  reaction nodes that more than one pool passes through. Each is written once,
  from the flux of every pool through it. It is 0 unless the bundle was built
  with `LNG_POOL_SHARED=1`.
- **`residual`** (additive, 2026-10-15): the solve's final residual. A solve
  that stopped at the iteration cap (`converged: false`) is still returned; this
  says how far from steady it stopped.
- **`depletion_own_product`** / **`own_source_depletions`** (additive,
  specs/061): the own-product depletion rule that ran
  (`"suppress_only_source"` by default since 2026-10-15; `"full"`,
  `"suppress_only"`, `"suppress_only_entity"`) and how many depletion edges
  were gated on their source copies. **`recycled_twin_rule`** /
  **`recycled_twin_edges`** (specs/059): `"on"` by default since 2026-10-15.
- **`carrier_loop_rule`**, **`carrier_loop_forms`**,
  **`carrier_loop_edges_dropped`**, **`carrier_loop_twin_edges`** (additive,
  specs/049): `"closed_twin"` (default since amendment 6), `"off"`,
  `"closed"`, `"closed_supply"`, `"all"`, or `"<mode>: no carrier table"`
  (which solves as `"off"`); how many carrier-loop free forms lost an in-edge
  (from their own loop, or under `closed_supply`/`closed_twin` a depletion
  edge), and how many edges that was. Influence scores are computed on the
  network the solve ran, with these edges dropped and the twin edges added. `closed_supply` also drops depletion in-edges into
  closed free forms, and counts them. `closed_twin` also gives each closed
  free form an in-edge from every root copy of the same entity, so it
  follows its own gene when that gene is perturbed; `carrier_loop_twin_edges`
  counts them. Under `closed`, the free form of a loop that nothing outside it
  produces then reads baseline (unless a silo bridge or boundary edge still
  feeds it). Under `all`, every listed free form reads only
  its producers outside the loop. An observation still wins.
- **`scc_order`**, **`scc_order_ties`** (additive, specs/013): the Gauss-Seidel
  sweep order inside cyclic components (`flow` by default; `structure`, or
  `label` for the old uuid order), and how many adjacent reactions it could not
  tell apart (they keep uuid order).
- **`scc`** (additive, specs/017): how the cyclic components were resolved —
  `method` (the `DS_SCC_METHOD` in force), `pooled`, `iterated`,
  `fallback_negative`, `fallback_inconsistent`, `pooled_nodes`. Under the
  default `fixed_point`, `pooled` is 0 and `iterated` counts the cyclic
  components. A benchmark arm that expects pooling checks `pooled > 0` here
  instead of trusting the environment. With `DS_SCC_BREAK_ROLES` set
  (specs/018) it also carries `break_roles`, `closures_catalyst` /
  `closures_assembly` / `closures_depletion` (edges read at their entry value)
  and `cyclic_before` / `cyclic_after` / `largest_after` (component census
  before and after recomputation).
- **`influence_scores`**: `uuid → Σ|∂F/∂input|` at the solution, computed with
  the same propagator as the solve (so inhibitors carry real influence). Use for
  a "what drove this?" view; higher = more influential.

## Notes for the WebsiteAngular client
- Treat this contract as the integration boundary. Keep all Angular code in the
  WebsiteAngular workspace; call these endpoints over HTTP.
- Typical flow: `GET /api/pathways` → user picks one → `POST /api/parse
  {pathway_id}` (keep the `network_id`) → render nodes/edges → user perturbs →
  `POST /api/solve {network_id, observations}` → overlay `node_activities × 100`
  and rank drivers by `influence_scores`.
- On a `network_id` miss after a backend restart, re-`parse` and retry.
