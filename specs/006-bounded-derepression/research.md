# Research: Bound de-repression on purpose

Phase 0. Four questions, all resolved by reading the code and the arithmetic.

## R1 — The pattern already exists, one branch over

`DS_DEPLETION_H_MAX` (default 10) is an explicit, named, documented bound on
de-repression for **depletion** edges, with a comment saying exactly that:
"Capped at DS_DEPLETION_H_MAX (default 10) to bound runaway de-repression."

The **inhibitor** branch has no equivalent. It has a hard-coded
`clamp(result, 0, 10.0)` with the comment "Cap H to prevent catastrophic
de-repression from multiple knockouts" — so the problem was known — but no
configuration control, no per-inhibitor bound, and the per-inhibitor maximum
is set implicitly by the division guard.

**Decision**: give inhibitor edges the same treatment depletion edges already
have. This is a consistency fix with a precedent in the same function, not a
new concept.

## R2 — Exactly what the guard does, and why it cannot simply be lowered

`h = (bl + ε)/(x + ε)`, `bl = 0.01`.

| ε | max de-repression (x=0) | at x = bl/2 | max suppression (x=1) |
|---|---|---|---|
| 1e-9 | **10,000,000×** | 2.000 | 0.0100 |
| 1e-4 | 101× | 1.980 | 0.0101 |
| **1e-3 (current)** | **11×** | 1.833 | 0.0110 |
| 1e-2 | 2× | 1.400 | 0.0198 |

Three separate effects move together:

1. **The de-repression ceiling** (column 2) — the defect.
2. **The curve shape at intermediate values** (column 3) — ε compresses the
   whole response, not only its end.
3. **The suppression floor** (column 4) — *smaller* ε means *stronger*
   suppression.

So the probe that raised ε from 1e-3 to 1e-2 weakened de-repression (11× →
2×), flattened the intermediate response (1.833 → 1.400) **and** weakened
suppression (0.0110 → 0.0198). Three changes, one number, one measured gain.

**Decision**: lowering ε alone is not viable — at 1e-9 the ceiling becomes
ten million. ε can only be reduced *together with* an explicit ceiling, which
is why FR-001 and FR-002 are one piece of work rather than two.

## R3 — The attribution experiment is a 2×2, and it discriminates

With the two controls separated, the arms isolate each effect:

| arm | ε | ceiling | isolates |
|---|---|---|---|
| **A** | 1e-3 | 11 | current behaviour exactly |
| **B** | 1e-9 | 11 | the curve-shape and suppression effects of ε, ceiling held |
| **C** | 1e-9 | 2 | the de-repression ceiling, ε held |
| **D** | 1e-2 | none | the original probe, for continuity |

The design discriminates because the effects push opposite ways: going from
A to B *strengthens* suppression (0.0110 → 0.0100). If the probe's gain came
from weaker suppression, **B scores worse than A and C does not recover it**.
If it came from the ceiling, B ≈ A and C captures the gain. A result where
both move is also readable, as a proportion.

**Alternatives considered**: sweeping ε more finely — rejected, it cannot
separate effects no matter how many points it has, which is the whole problem.
Holding ε and only adding the ceiling — insufficient, because it leaves the
curve-shape and suppression effects unmeasured and FR-002 unmet.

## R4 — Per-inhibitor or per-reaction, and the compounding question

Current behaviour bounds **per reaction** (`clamp(result, 0, 10.0)`) after
multiplying the individual factors, and applies **no per-inhibitor bound**
beyond what ε implies. Nothing bounds compounding **across** reactions along
a path, which is where the observed 25× and 100× predictions come from: each
reaction may legitimately return up to 10×, and a path multiplies them.

**Decision**: bound **per inhibitor**, and keep a per-reaction bound as well.
Per-inhibitor is the level at which the modelling claim lives — "removing
*this* inhibitor can raise its target at most k-fold" — and it makes the
claim checkable in isolation. An end-to-end bound is rejected: it would make
a reaction's output depend on path history, which the propagator has no
representation for and which would break the node-local evaluation the SCC
solver relies on.

Compounding across reactions therefore remains possible and is a stated
limitation, not an oversight: a chain of genuinely de-repressed steps should
compound. What changes is the per-step magnitude, from 11× to whatever the
evidence supports.

## Open questions carried into implementation

1. Does the measured gain attribute to the ceiling, to suppression, or both?
   That is US1 and it is unrun.
2. What ceiling value does the 92-pathway catalog support? The 2× candidate
   comes from a sweep on the smaller set and may not survive.
3. Does the `DS_DEPLETION_H_MAX` default of 10 deserve the same scrutiny?
   Out of scope here, but the same argument applies to it and it should be
   recorded as a follow-on rather than quietly changed alongside.
