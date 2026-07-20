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
   `logic_network` and `uuid_mapping` (required), `set_mappings` (optional).
3. **No body** — falls back to the bundled sample network.

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
  "pathways": [ { "id": "...", "name": "...", "members": ["uuid", ...] } ]
}
```
Node display names are enriched from the Reactome ContentService when the
generator only provided stable ids (degrades gracefully if that service is
down — names fall back to the stable id).

**`network_id`** is the key to the parse-once/solve-many flow: hand it to
`/api/solve` instead of re-shipping the whole network on every solve. It is an
in-memory, per-process cache — it does **not** survive a server restart, so a
client must be prepared to re-`parse` on a `404`-style cache miss (a solve with
an unknown `network_id` falls through to the inline `network`, then the sample).

### `POST /api/solve`
Solve the steady state under a set of perturbations. JSON body:
```json
{
  "network_id": "pw:Cell_Cycle_Checkpoints_R-HSA-69620",
  "observations": { "<node_uuid>": [80.0, 0.9] }
}
```
- **Network source** (in priority order): `network_id` (from a prior parse) →
  inline `network` (same shape as the parse response's nodes/edges/pathways) →
  bundled sample.
- **`observations`**: object keyed by node `uuid`; each value is
  `[activity, confidence]` where `activity` is **0–100** and `confidence` is
  `0–1`. Confidence `> 0` pins the node as a hard constraint; `0` ignores it.

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
