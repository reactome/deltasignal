# Contract: SCC resolution diagnostics (additive)

`POST /api/solve` response gains one optional object; every existing field is
unchanged.

```json
"scc": {
  "method": "pool_parity",
  "pooled": 12,
  "iterated": 3,
  "fallback_negative": 0,
  "fallback_inconsistent": 3,
  "pooled_nodes": 1840
}
```

- `method` — the resolved `DS_SCC_METHOD`.
- `pooled` / `iterated` — cyclic components resolved by the pool rule / by the
  damped fixed point (or minimiser). Under `fixed_point` and `minimize`,
  `pooled` is 0.
- `fallback_negative` — components sent to iteration because they contain an
  internal negative edge (`pool` only).
- `fallback_inconsistent` — components sent to iteration because their internal
  signs are inconsistent (odd negative cycle).
- `pooled_nodes` — non-pinned member nodes written by the pool rule.

`SolverResult.diagnostics` carries the same keys prefixed `scc_`. The CLI
`solve_provenance` already records `scc_method`.
