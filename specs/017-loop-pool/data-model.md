# Data model: Loop as a conserved pool

No persistent data. Solve-time structures only.

| entity | representation | rules |
|---|---|---|
| Component `c` | `comp_id::Vector{Int}` (node → component), `comp_rxns[c]` (member reaction indices) — existing | cyclic iff `length(comp_rxns[c]) > 1` or a self-loop (existing test) |
| Member node set | `{i : comp_id[i] == c}` | pinned members (`i in obs_set`) are never written |
| External fold `e_r` | `compute_reaction_output_vec(xpool, r) / baseline_vec[r.target_idx]` | `xpool` = `x` with non-pinned member nodes at baseline |
| Entry | member reaction with `e_r != 1` (tolerance 0, exact) | its target node is the BFS root for parity |
| Pool fold | product of entry folds (sorted by value), then `clamp(·, 0, 100)` | 0 absorbs; cap equals hill_sat's 100x |
| Parity sign `s_e[m]` | ±1 per (entry, member), BFS over internal signed edges | conflict ⇒ component inconsistent ⇒ fall back |
| Member value | `clamp(baseline[m] * pool_m, 0, 1)` where `pool_m = prod_e e_fold^(s_e[m])` (`pool`: all `s = +1`) | written once; no iteration |
| Mode | `DS_SCC_METHOD ∈ {fixed_point, minimize, pool, pool_parity}` | validated at `solve_scc_ordered!` entry via the existing check |
| Fallback | `pool`: any internal negative edge; both: inconsistent parity; both: component with no member reactions | falls back to the existing damped fixed point for that component only |
| Report | `(pooled, iterated, fallback_negative, fallback_inconsistent, pooled_nodes)` | into `SolverResult.diagnostics` and the API `scc` object |
