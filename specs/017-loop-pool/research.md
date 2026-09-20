# Research: Loop as a conserved pool

All five unknowns resolved from the code in `src/solvers/steady_state.jl`
(`solve_scc_ordered!`, lines 569–820 on this branch) and
`src/core/reaction_model.jl`.

## R1. Evaluating a member reaction with in-component inputs at baseline

**Decision**: build `xpool = copy(x)`; for every member node `m` of component
`c` (`comp_id[m] == c`) that is not in `obs_set`, set `xpool[m] =
baseline_vec[m]`; evaluate `compute_reaction_output_vec(xpool, r; config)` for
each member reaction `r` and take `e_r = out / baseline_vec[r.target_idx]`.

**Rationale**: `compute_reaction_output_vec` reads activators
(`x[rxn.activator_indices[k]]`), inhibitors (`x[rxn.inhibitor_indices[k]]`,
line 1714) and depleters (`x[rxn.depletion_indices[k]]`, line 1735) from the
same vector, so one overlaid copy covers every input class. External inputs
keep their solved upstream values (the component is visited in topological
order, so they are final); internal ones read baseline, i.e. fold 1, which is
exactly "no internal gain". Pinned members keep their pinned value and
therefore act as entries.

**Alternatives**: the `supply` overlay used by `DS_SCC_BREAK_CATALYST` freezes
only catalyst back-edges at the component's *entry state* (which is baseline
for interior nodes anyway — specs/016 §5); it is a per-edge mechanism and
cannot express "all internal inputs". Rejected.

## R2. Detecting internal negative edges per component

**Decision**: for member reaction `r` (target `t`), an internal negative edge
exists iff any `i in r.inhibitor_indices` or `i in r.depletion_indices` has
`comp_id[i] == comp_id[t]`. Computed once per component from
`comp_rxns[c]` (built at line 690).

**Rationale**: these are the only two negative input classes the propagator
has (`pos_neg == "neg"` edges land in `inhibitors` or, for `edge_type ==
"depletion"`, in `depletions` — `reaction_model.jl:205–210`). Activators are
positive regardless of edge_type. `comp_neg_frac[c]` (line 697) already counts
negative in-component edges for the `DS_SCC_NEG_MODE` layer and can be reused
as the gate: `comp_neg_frac[c] > 0`.

## R3. Parity with multiple entries and multiple paths

**Decision**: build the component's internal signed node graph once:
`adj[src] = [(tgt, sign)]` for every member reaction `r` and every input index
`i` with `comp_id[i] == c`, `sign = +1` for activators, `-1` for inhibitors and
depleters. For each entry node `e` (target of a reaction with `e_r != 1`), BFS
from `e` assigning `s[e] = +1` and `s[v] = s[u] * sign(u→v)`; a node reached
with two different signs from the same entry marks the component
*inconsistent* → fall back to the fixed point. Member value:
`baseline[m] * prod_e e_fold^(s_e[m])`.

**Rationale**: in an SCC every node is reachable from every entry, so every
member gets a sign; inconsistency is exactly an odd negative cycle reachable
from that entry, which is genuine negative feedback and is the case the spec
leaves to iteration. Multiple entries combine as a product per Adam's rule,
each with its own parity.

**Determinism**: fold multiplication is not associative in floating point, so
the product is taken over entry folds sorted by value (ties by exponent sign),
never by node index — bit-identical under relabelling (FR-008).

## R4. Surfacing pooled / iterated counts

**Decision**: `solve_scc_ordered!` returns a third value, a `NamedTuple`
`(pooled, iterated, fallback_negative, fallback_inconsistent, pooled_nodes)`;
the call site (line 925) stores them in `SolverResult.diagnostics` under
`"scc_method"`, `"scc_pooled"`, `"scc_iterated"`, `"scc_fallback_negative"`,
`"scc_fallback_inconsistent"`, `"scc_pooled_nodes"`. The API solve response
gains an additive `"scc"` object carrying the same keys (documented in
`docs/API.md`; nothing existing changes). The CLI provenance already records
`scc_method`.

**Rationale**: `diagnostics::Dict{String,Any}` exists for this; the API does
not currently echo diagnostics at all (`server.jl:762–769`), so an additive
object is the smallest contract change. A benchmark arm verifies the mode from
the server env *and* from a probe solve's `scc` counts, so an arm cannot run
the old solver by accident (constitution V).

## R5. Reusable fixtures

- `test/test_loop_elasticity.jl::posloop()` — U → A, A → r1 → B → r2 → A: the
  positive two-node loop with one external drive (US1 scenarios 1–2).
- `test/test_solver_determinism.jl::relabel`, `RELABELLINGS`,
  `solve_by_original` — the relabelling harness (FR-008, SC-002); its
  `cyclic_net()` has a single root and one component.
- New fixtures needed: four-entry AND loop (scenario 3), OR-entry loop
  (scenario 4), the X → M, M ⊣ T parity fixture (US3), an odd-negative-cycle
  fixture (US3 scenario 3), a three-component network for the count report
  (US2 scenario 2). All go in a new `test/test_loop_pool.jl` (the repo's rule:
  assertions in a file that has them).

## Constitution check

- I (processing is the solver's job): no generator change. ✓
- II (measure, then claim): both variants measured on the wide set, both axes,
  held-out, p; predictions committed first. ✓
- III (negative results are results): the decision rule in SC-004 makes a
  negative outcome a recorded result. ✓
- IV (default-OFF): `DS_SCC_METHOD` default `fixed_point` unchanged;
  differential check against the pre-feature tree in the tasks. ✓
- V (honest reporting): fallback counts in diagnostics and API. ✓
- VI (API boundary): additive field only. ✓
