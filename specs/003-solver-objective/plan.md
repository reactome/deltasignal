# Implementation Plan: Solve the objective the model specifies

**Branch**: `003-solver-objective` | **Date**: 2026-09-10 | **Spec**: [spec.md](./spec.md)

## Summary

Replace the damped fixed-point iteration inside cyclic components with a
minimisation of the objective the design specifies. Keep everything else:
the SCC decomposition, the topological ordering, and the exact single-pass
evaluation of acyclic nodes are correct and stay untouched.

## Technical Context

**Language**: Julia 1.10
**Objective**: `Σ ωᵢ(xᵢ − yᵢ)² + μ‖x − F(x)‖² + γ‖x − x₀‖²`
**Dependencies**: none new. `Optim` and `ForwardDiff` are already imported in
`src/DeltaSignal.jl`, and `compute_reaction_output_vec` is declared
`where {T<:Real}` over `AbstractVector{T}` — it is already autodiff-ready.
That scaffolding was laid deliberately and never used.

**Problem size** (measured, not assumed — this decides feasibility):

| pathway | nodes | components | cyclic | largest cycle | nodes in cycles |
|---|---|---|---|---|---|
| TP53 R-HSA-3700989 | 2,307 | 1,472 | 1 | **836** | 836 |
| Cell Cycle Checkpoints | 1,490 | 717 | 4 | 643 | 777 |
| PIP3 | 1,187 | 1,042 | 2 | 97 | 147 |

The optimisation is over **836 variables at worst**, not 2,307 and not the
16k of a whole catalog. That is comfortably inside what LBFGS with
ForwardDiff gradients handles interactively, which is what makes SC7
achievable.

**NEEDS CLARIFICATION**: none.

## Constitution Check

| Principle | Status | Note |
|---|---|---|
| I. Processing is DeltaSignal's job | **PASS — the clearest case yet** | The networks are fine; the unperturbed TP53 network reaches a consistent state with shortfall 5.2e-18. The failure is entirely in how we process it. |
| II. Measure, then claim | **PASS** | Comparisons share one catalog build (002 R1). Every arm reports per-pathway attribution. |
| III. Negative results are results | PASS | If the minimisation does not beat 407, that is recorded with its numbers, not buried. |
| IV. New behaviour ships default-OFF | **JUSTIFIED DEPARTURE** | See Complexity Tracking. |
| V. Honest solver reporting | **PASS — this feature is that principle** | It replaces a pass/fail verdict, currently false for 122 cases, with a reported shortfall that is always meaningful. |
| VI. API boundary | N/A | No contract change; the response gains a shortfall figure where it had a boolean. |

## Approach

### Stage 1 — Objective and gradient, tested in isolation

Build the objective as a pure function of the component's node vector, with
the rest of the network's already-final values held fixed. Verify
`ForwardDiff` differentiates it cleanly through the propagator — including
the `min`/`max` branches in the aggregators, which are the parts most likely
to produce a zero or undefined derivative.

**Risk to check first, before any integration**: `min`, `max` and `clamp` are
not differentiable at their switching points, and the propagator uses all
three. If the gradient is unusable, the fallback is a derivative-free method
over the same objective — slower, but the objective is what matters, not the
optimiser.

### Stage 2 — Replace only the cyclic inner loop

`solve_scc_ordered!` keeps its structure. For a singleton acyclic component
the exact evaluation stays. For a cyclic component, minimise the objective
over that component's nodes, **warm-started from the current damped
iteration's result** — a good starting point costs nothing and the existing
iteration produces one.

Observations become weighted terms rather than skipped indices: the
`t in obs_set && continue` skip is what makes an observation absolute, and it
is what FR4 requires removing.

### Stage 3 — Wire the weights

`μ`, `γ` and per-observation `ω` become live. They are already carried in
`SteadyStateParams` and observation tuples, so this is connection, not new
plumbing. Starting values follow MP-BioPath's relative scale, where
measurements dominate model consistency by roughly 10³ and the baseline prior
is weakest.

### Stage 4 — Measure

Against MP-BioPath and the structural baseline on identical cases, sharing
one catalog build. Report per-pathway attribution, and separately for the 442
currently-stable and 122 currently-broken cases — SC5 requires showing the
stable set did not regress.

## Non-goals, and why

- **No JuMP/Ipopt.** MP-BioPath's route, but it needs the model expressed as
  constraints, which means re-encoding aggregation logic that already exists
  and is verified in Julia. Optim over the existing propagator reuses it.
- **No binarisation.** MP-BioPath splits n-ary AND into pairwise
  `_PSEUDONODE_` chains for solver ergonomics. Multiplication is associative
  so it is exact for them, but our saturation is applied to the summed
  log-fold, so pairwise splitting would apply the tanh once per pair and
  attenuate harder — the very failure being fixed in 002. FR7 forbids it.
- **No change to the aggregation rules.** Out of scope; 002 owns the clamp.

## Project Structure

```
specs/003-solver-objective/
├── spec.md, plan.md, checklists/    # written
├── research.md                       # Phase 0: autodiff viability, weights
└── quickstart.md                     # Phase 0: reproducing every figure

src/solvers/steady_state.jl          # the only change site
```

No `data-model.md` or `contracts/`: no new entities, no interface change.

## Complexity Tracking

| Departure | Why needed | Simpler alternative rejected because |
|---|---|---|
| Principle IV: new behaviour ships default-OFF. A solver that only minimises when asked leaves the default path broken. | The convention protects against modelling changes of unknown direction. Here the current default is measurably indefensible — 122 cases return a value determined by when we stopped computing. Leaving it on by default is the risk. | Shipping behind an opt-in flag means the benchmark, the API and the UI all keep the broken solve. Mitigation: keep the fixed-point path reachable by flag for A/B, and gate the switch of default on SC5 (the stable subset must not regress). |

## Open risks

1. **Non-differentiable branches** in the aggregators. Checked in Stage 1
   before anything is built on it.
2. **Local minima.** The objective is non-convex, so the result depends on the
   start. Warm-starting from the feed-forward pass is both the natural choice
   and a reproducible one, which SC's edge case on determinism requires.
3. **Runtime.** 836 variables is tractable, but the objective evaluates the
   whole component's propagator per iteration. If SC7 is threatened, the
   fallback is to minimise only over components that fail to converge — a
   hybrid keeping the cheap path where it already works.
4. **It may not be enough.** The clamp (002) accounts for 21 cases and this
   for an unknown number; together they may still not reach 407. The spec
   says exceed, and if the gap survives both fixes, the residual becomes the
   next feature rather than a disappointment.
