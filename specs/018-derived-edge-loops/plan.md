# Implementation Plan: Derived-edge recycling closures

**Branch**: `feat/018-derived-edge-loops` | **Date**: 2026-09-19 | **Spec**: [spec.md](spec.md)

## Summary

`DS_SCC_BREAK_ROLES=<comma list of catalyst|assembly|depletion>`. In
`index_reactions`: first Tarjan pass as today; every edge of a listed role
whose source and target share a component is marked a closure; a second
Tarjan pass runs on the adjacency without those edges and its components drive
the topological order and the per-component iteration. Closure activators join
the existing `activator_break` (read through `supply`, the component-entry
state); closure depletions get a new `depletion_break` read the same way; the
acyclic single-pass evaluation passes `supply = x` so a closure into an acyclic
node reads the current (entry) state. Counts per role and the component census
before/after go into `SolverResult.diagnostics` and the API `scc` object.

## Technical Context

Julia 1.10; no new dependencies. Files: `src/core/reaction_model.jl`
(`_break_roles_env`, `DS_BREAK_ROLES`, two-pass census, `depletion_break`
field, supply read in the depletion term), `src/solvers/steady_state.jl`
(`use_supply` extended, singleton branch passes `supply = x`, stats sink,
diagnostics), `src/api/server.jl` + `docs/API.md` (additive keys),
`test/test_scc_break_roles.jl` (26 assertions).

## Constitution Check

I ✓ solver-side only · II ✓ measured against the current-tree control, both axes, held-out, p, pre-registered · III ✓ SC-004 names the negative outcome · IV ✓ default unset = bit-identical (asserted); legacy `DS_SCC_BREAK_CATALYST` unchanged · V ✓ closures and census reported · VI ✓ additive API keys.

## Design notes (from research)

- Closure marking uses the FIRST-pass components (the welded ones); recomputation may leave some closure edges inside a second-pass component when another cycle survives through them — they still read `supply`, so no feedback, and the surviving cycle is iterated as before.
- `fwd_adj` is a `Set` per node, so an s→t pair with several edges is dropped only if every edge on the pair is a closure (rebuilt from the surviving edges).
- Determinism: closure marking depends only on component membership, which is label-free; tested.
- `DS_SCC_BREAK_CATALYST` (legacy, no recomputation) is left byte-identical; `catalyst` in the role list is the recomputing form and is measured separately.
