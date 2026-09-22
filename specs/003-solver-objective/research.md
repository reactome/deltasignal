# Phase 0 research: autodiff viability, and why the γ term is the whole point

**Created**: 2026-09-18
**Feature**: 003-solver-objective
**Status**: Phase 0 complete — the gating risk is cleared, proceed to Stage 2

## 1. The gating risk is cleared: ForwardDiff differentiates the propagator

plan.md Stage 1 named this as the thing to check before building anything:
`min`, `max` and `clamp` are not differentiable at their switching points and
the propagator uses all three.

Measured on a two-reaction cycle, gradient of `‖F(x) − x‖²`:

| x | gradient finite | nonzero entries | ‖g‖ |
|---|---|---|---|
| 0.0 | yes | **0 / 7** | 0.0 |
| 1e-9 | yes | 6 / 7 | 2.00e-9 |
| 0.01 (baseline) | yes | 6 / 7 | 9.32e-17 |
| 0.5 | yes | 6 / 7 | 1.414 |
| 1.0 | yes | 4 / 7 | 1.58e-9 |

Finite everywhere, informative in the interior. **No fallback to a
derivative-free method is needed.** The `where {T<:Real}` scaffolding on
`compute_reaction_output_vec` works as intended.

## 2. The all-zero root is a global minimum when γ = 0

The `x = 0.0` row above is not an artifact. It is the observation, stated
exactly: *"if we had a loop and the numbers should all be 1 or higher, in the
way we are doing the optimization if all the nodes of the loop are zero the
math works even though all the inputs to the loop are either one or greater."*

With `L = μ‖F(x) − x‖² + γ‖x − x₀‖²`, evaluated at the all-zero state against
the correct baseline state:

| γ | ‖∇L‖ at all-zero | stationary? | L(all-zero) | L(correct) |
|---|---|---|---|---|
| **0.0** (what ships) | **0.0** | **yes** | **0.0** | 4.35e-33 |
| **0.1** (spec default) | 5.29e-3 | no | 7.0e-5 | 4.35e-33 |

At γ = 0 the all-zero state is a **stationary point and a global minimum**,
tied with the right answer. Nothing in the objective prefers the correct
state. At γ = 0.1 the stationarity is broken and the correct state is ~10²⁸
times better.

**The design anticipated this failure mode and the implementation dropped the
one term that prevents it.** FR5 is therefore not a truthfulness requirement
about provenance; γ is load-bearing.

## 3. FR3 is violated harder than the spec recorded

spec.md FR3: *"Repeating a solve with a larger computational budget MUST NOT
materially change the answer."* Re-measured 2026-09-18 on the 92-pathway
catalog, comparing two **labellings of one network** (verified isomorphic —
identical stable-id edge multiset for 92 of 92 pathways):

| `DS_MAX_ITERS` | predictions differing between labellings | of which held-out |
|---|---|---|
| 500 (default) | 14 — all TP53 | 0 |
| **10,000** | **25** — TP53 14 + RUNX2 11 | **11** |

A 20x budget made it **worse**, and pushed the instability into the held-out
half, which had been clean. Accuracy barely moved (83.36→83.40, 83.42→83.42),
so this is not a convergence trade-off: the extra iterations let two orderings
of the same graph drift further apart. Iterating longer does not find the root
— there is no unique root to find, which is section 2.

## 4. The answer currently depends on node NAMES

Same experiment, five labellings of one network:

| labelling | accuracy | vs original | net | held-out net |
|---|---|---|---|---|
| original | 83.36% | — | — | — |
| perm | 83.42% | 14 | +14 | +0 |
| p101 | 83.41% | 53 | +13 | +5 |
| **p202** | **83.01%** | **112** | **−84** | **+0** |
| p303 | 83.47% | 28 | +28 | +0 |

**Spread 0.47pp from node names alone.** Cause: node indices come from
`collect(keys(network.nodes))` (Julia `Dict`, i.e. UUID hash order), and that
order drives the Gauss-Seidel sweep inside a cyclic component, which selects
which root it lands in. Zero acyclic pathways churn; all 10 churning pathways
contain a cycle.

A minimiser does not have this problem: the answer is defined by the objective,
not by the visit order. `DS_SCC_SWEEP=jacobi` (implemented, see
specs/013) removes the label-dependence but still does not give the loop a
unique answer — it is a symptom fix and should not be merged as if it were
more.

## 5. The forward model also drifted from the design

Recorded here because it changes what the existing A/B results mean, not
because this feature changes it.

| design (Sections 2.2–2.3) | what ships |
|---|---|
| activators weighted geomean `A = exp(Σwᵢ log(x̃ᵢ+ε))` | **`DS_AND_MODE=hill_sat`** — product of fold-changes, capped at 100 |
| inhibitors `H = ∏ 1/(1+βⱼxⱼ^mⱼ)`, β=0 | **`DS_INHIBITION_MODE=divide`** — `H = (baseline+ε)/(x+ε)` |
| output `y_r = s_r^h/(s_r^h + K_r^h)` | **absent** — `reaction_model.jl:1512` returns `clamp(A·H·H_dep·L, 0, 1)` |
| `x₀ = 0.01`, UI 0–100 | matches |
| sensitivity `α(x)` | present, `s=0` default ⇒ identity |

Each substitution has its own spec folder and measured A/B. But **every one of
those A/Bs was run inside the fixed-point solve**, so they establish which
forward equation works best *when iterated to a fixed point* — not which works
best inside the specified minimisation. Once this feature lands, the AND mode,
the inhibition form and the missing output Hill should be re-compared under it.
Until then the current defaults are validated against the solver we have, not
against the design.

## 6. γ has a floor AND a ceiling, and the spec default is above the ceiling

Stage 2 is implemented (`DS_SCC_METHOD=minimize`). Two measurements fix γ.

**Ceiling — γ distorts well-determined answers.** Two-reaction loop, minimiser
vs the fixed point it warm-starts from:

| γ | `U1` knocked out → A (UI) | `U1` at 80x → readout (UI) |
|---|---|---|
| 0.0 | 2.3e-10 | 100.0 |
| 1e-9 | 2.3e-10 | 100.0 |
| 1e-6 | 3.9e-6 | 99.999 |
| 1e-4 | 4.0e-4 | 99.91 |
| 1e-2 | 0.037 | 91.1 |
| **0.1 (spec default)** | **0.23** | **51.5** |

At γ ≤ 1e-9 the minimiser reproduces the fixed point exactly, which validates
the implementation — it finds the same root, it is not doing something else.
At the spec's γ = 0.1 a 100x perturbation reads 51x at the readout. **γ = 0.1
is not usable in this formulation**; the per-node prior competes directly with
per-reaction consistency.

**Floor — below ~5e-7 the prior cannot escape a collapsed loop.** Warm-started
at the near-zero state the fixed point actually reaches, supply pinned at
baseline, so the correct answer is 1.0 UI:

| γ | escapes? | A (UI) | residual |
|---|---|---|---|
| 0.0 | **no** | 1e-28 | 5.0e-10 |
| 1e-9 | no | 1e-28 | 5.0e-10 |
| 1e-7 | no | 1e-28 | 5.0e-10 |
| **1e-6** | **YES** | **1.0** | **3.3e-17** |
| 1e-4 | YES | 1.0 | 7.8e-17 |

Identical from starts of 1e-30, 1e-10 and 1e-7. The threshold matches the
predicted floor `g_tol / (2·x₀)` = 1e-8 / 0.02 = 5e-7: below it the escape
gradient `2γx₀` is smaller than the optimiser's gradient tolerance, so LBFGS
stops where it started.

**Note the collapsed state's residual is 5.0e-10 — below the 1e-8 tolerance.**
The fixed-point solver therefore reports `converged = true` on it. The loop
returns a wrong answer and says it succeeded, which is Principle V's case
exactly.

**A caveat on an earlier probe, recorded so it is not repeated.** Starting the
minimiser at *exactly* zero escapes at every γ including 0 — but only because
`Fminbox` refuses to begin on the box boundary and nudges the point interior,
from where consistency alone suffices. That is an artifact of the optimiser,
not evidence about γ. The table above starts strictly inside the box for that
reason.

**Chosen: γ = 1e-6**, three orders above the floor and five below the spec
default. Consistent with the project's epsilon-sizing rule (baseline is 0.01,
so a tie-breaking term belongs orders of magnitude below it). To be confirmed
at catalog scale in Stage 4.

## Decisions for Stage 2

- **Optimiser**: LBFGS via `Optim.jl` with `ForwardDiff` gradients, box
  constrained to [0,1] (`Fminbox`). Justified by section 1.
- **γ ships non-zero.** Section 2 makes γ = 0 indefensible: it makes the
  degenerate state a global minimum. Spec default γ = 0.1 is the starting
  value, swept in Stage 4.
- **Warm start** from the existing feed-forward/damped pass, per plan.md
  Stage 2 — a reproducible starting point, which the spec's determinism edge
  case requires.
- **Acyclic components keep their exact single-pass evaluation** (FR6). Only
  cyclic components change.

## 7. NEGATIVE RESULT — the minimiser loses to the fixed point at every γ

Stage 4, run on the 92-pathway catalog (`cat_os`), scored against the wide
curator set with the held-out split. LM optimiser, observations still hard
pinned (FR4/Stage 3 not implemented).

| arm | accuracy | macro-F1 |
|---|---|---|
| **fixed point (baseline)** | **83.36%** | **0.7998** |
| minimise, γ=1e-6 | 82.85% | 0.7952 |
| minimise, γ=1e-4 | 82.61% | 0.7917 |
| minimise, γ=1e-2 | 82.48% | 0.7855 |

Both metrics decline monotonically in γ. Held-out, paired:

| γ | net | fixed | broke | p | pathways moved | readouts |
|---|---|---|---|---|---|---|
| 1e-6 | −71 | 1 | 72 | <0.0001 | 3/70 | 20 |
| 1e-4 | −87 | 3 | 90 | <0.0001 | 6/70 | 30 |
| 1e-2 | −148 | 125 | 273 | <0.0001 | 24/70 | 107 |

Distributed across many pathways and readouts, so this is not the
single-readout concentration artifact that invalidated LNG #89's p-value.

**The γ=1e-6 arm is the informative one.** There the prior is nearly inert, so
the objective reduces to minimising `‖F(x) − x‖²` — of which the fixed point is
a *global* minimum, residual ~1e-9. The minimiser still lands 1-fixed/72-broke
worse. It is therefore finding **different zero-residual roots**, and the
damped Gauss-Seidel iteration's choice among them is the better predictor.

**How it fails.** 136 of 137 broken cases had truth = NORMAL and the fixed
point got them right; the minimiser moved them off baseline (103 → DOWN,
33 → UP). Only 4% are full collapses — most are small drifts (`1.0 → 0.816`)
that cross the NORMAL band edge. That is **false change**, already the largest
error class at 48.2%. A weak prior lets unperturbed nodes wander; a strong one
(γ=1e-2) flattens genuine changes instead. Neither end wins.

**A prediction that failed, recorded so it is not re-made.** From the
false-change diagnosis at γ=1e-6 I predicted larger γ would help, because a
stronger baseline prior is exactly what holds unperturbed nodes at baseline.
It made things monotonically worse. The false-change mechanism was real and the
inference from it was wrong.

### What survives, and what does not

Survives — all independently measured:
- Loops have multiple roots; the solve is not a function of the graph
  (0.47pp from UUID relabelling alone, specs/013).
- At γ=0 the all-zero state is a stationary point AND a global minimum, tied
  with the correct answer (section 2). that mechanism was right.
- A collapsed loop reports `converged = true` at residual 5e-10.
- More iterations makes divergence worse, not better (section 3).
- LM makes the minimisation affordable: 814s catalog wall, same as the fixed
  point, against LBFGS's 29s on a single 840-node component.

Does not survive:
- **"Solving the specified objective will improve accuracy."** It does not, at
  any γ tested. The objective's optimum is a worse predictor than the
  iteration's arbitrary root.

### Not tested

Observations remain hard constraints. The design's `Σ ωᵢ(xᵢ − yᵢ)²` term
(FR4) is unimplemented, so the full L_SS has not had its turn. It is unlikely
to reverse −71 but it is a real untested part of the design.

### Recommendation

Do not make `minimize` the default. Keep it reachable (`DS_SCC_METHOD=minimize`)
because it is needed for the parameter-learning direction — learning requires a
stable `d(prediction)/dθ`, which the label-dependent fixed point does not
provide — and because the correctness findings above stand on their own. But
the accuracy case for it is refuted on current evidence.
