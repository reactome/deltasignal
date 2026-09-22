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

---

# Phase 3–4 results

## R5 — The attribution succeeded, and the answer is two-part

The gate first: **arm A is behaviour-preserving**. Every prediction,
`predicted_ui`, convergence flag and iteration count is identical to the
pre-change run across all 847 rows, so the arms below are interpretable.

| arm | ε | ceiling | correct | macro-F1 | vs A |
|---|---|---|---|---|---|
| **A** | 1e-3 | 11 | 491/742 | 0.5771 | — (current) |
| **B** | 1e-9 | 11 | **483** | 0.5687 | **−8** |
| **C** | 1e-9 | 2 | 497 | 0.5833 | **+6** |
| **D** | 1e-2 | none | **506** | 0.5986 | **+15** |

Reading the design:

- **A → B is −8.** Tightening ε with the ceiling held *hurts*. That isolates
  ε's non-ceiling effects — a sharper curve at intermediate values and
  stronger suppression — and they are **negative**.
- **B → C is +14.** Adding the ceiling with ε held is the largest single
  effect measured, and it is unambiguously the de-repression bound.
- **C → D is +9.** Same effective 2× ceiling, looser ε. That is ε's
  non-ceiling effects again, in the *helpful* direction.

**So both halves contribute and they are separable: the ceiling is worth
about +14, and a looser ε is worth about +8–9 independently of it.** The
original +15 probe was getting both at once, which is why it could not be
attributed. FR-004 is satisfied and the answer is not the single-cause story
R15 implied.

**This partially corrects R15/R16 in `specs/004-loop-handling`**, as FR-010
requires. Those entries attributed the whole gain to de-repression. The
de-repression ceiling is real and is the larger effect, but roughly 60% of
the probe's improvement came from it and the rest from ε's effect on the
curve and on suppression strength — a separate mechanism that was not
identified at all.

### The targeted subgroup — arm C is the cleanest result

| arm | KO→UP calls | correct | accuracy | wrong removed | right lost |
|---|---|---|---|---|---|
| A | 78 | 28 | 0.359 | — | — |
| B | 81 | 29 | 0.358 | 3 | 0 |
| **C** | **70** | **29** | **0.414** | **10** | **0** |
| D | 68 | 28 | 0.412 | 12 | **1** |

**Arm C removes ten wrong de-repression calls and loses none.** That is
precisely the behaviour the bound was designed for, and it is a better
outcome than arm D, which scores higher overall but pays one correct call for
its twelve. The subgroup view is the only place that distinction is visible —
FR-009 earning its place.

## R6 — FR-002 is met, and it exposes what ε actually was

With the ceiling fixed at 2 and ε swept 1e-6 → 1e-9: **zero differing
predictions** across all four arms. The guard is a guard over that range.

But the internal values differ on ~340 of 742 cases even there, and at the
**current default** the effect is large enough to move 8 predictions (A→B).
The reason is arithmetic:

| ε | ε as a share of baseline | h at x = bl/2 |
|---|---|---|
| **1e-3 (current default)** | **10%** | 1.833 |
| 1e-4 | 1% | 1.980 |
| 1e-6 | 0.01% | 1.9998 |
| 1e-9 | 0.00001% | 2.0000 |

**At the current default, ε is ten percent of baseline.** It was never a
guard at that value; it was shaping the entire interior of the response
curve, and the ceiling was only its most visible side effect. Below about
1e-6 it is genuinely only preventing division by zero.

That reframes the defect. It is not merely "a guard accidentally set the
ceiling" — it is "a guard was sized like a model parameter", and the ceiling
was one of three things it was silently deciding.

## Open, carried to US3

1. Which combination to adopt. Arm D scores best on this set but costs a
   correct call and leaves the ceiling implicit again. Arm C is the clean
   mechanism. A fifth arm — loose-ish ε with an explicit ceiling — is the
   obvious candidate and has not been run.
2. Whether any of it survives the 92-pathway catalog. Unrun, and per FR-007
   nothing is adopted until it is.
3. Whether the gain is threshold-adjacent (SC-006). Unrun.

## R7 — The best solution is no epsilon at all, and this is the same bug twice

On being shown the attribution, the design intent was restated: epsilon was only supposed
to avoid divide by zero errors."* It was. It wasn't doing that.

A divide-by-zero guard has to be small enough to be invisible. Internal
activities live in [0,1] and Float64 handles 1e-300, so a real guard could be
1e-12. At **1e-3 it is ten percent of baseline** — a model parameter wearing
a guard's name, silently deciding three things: the de-repression ceiling,
the shape of the whole interior curve, and the maximum suppression.

### This is the same defect as feature 002's, in the sibling parameter

`specs/002-upregulation-propagation/research.md` R6 said of the *other*
epsilon:

> `DS_HILL_SAT_EPS` defaults to 0.001, which is LARGER than the internal
> values where knockouts live (~0.0006), so it acts as a floor and lifts the
> whole network upward.

Same value, same diagnosis, different parameter. Commit 17bc6be changed
`DS_HILL_SAT_EPS` from `0.001` to `1e-5` and **never checked whether 0.001
appeared anywhere else doing the same damage.** It did.

| parameter | default | share of baseline | status |
|---|---|---|---|
| `DS_HILL_SAT_EPS` | 1e-5 | 0.1% | fixed in feature 002 |
| `DS_INHIBITOR_EPS` | 1e-3 | **10%** | the defect here |

The general lesson, worth more than either fix: **a constant whose name
claims numerical safety needs a stated relationship to the scale it guards.**
Both of these were 1e-3 against a baseline of 0.01, almost certainly by
copying, and both quietly became model parameters.

### The fix: remove epsilon from the interior entirely

    h = min(ceiling, baseline / max(x, baseline / ceiling))

The **ceiling is the divide-by-zero guard** — `max(x, baseline/ceiling)` can
never be zero — so no epsilon appears anywhere. Above `baseline/ceiling` the
value is exactly `baseline/x` with no distortion; below it, the bound binds.

| x | exact `bl/x` | old (ε=1e-3) | new (ceiling 2) |
|---|---|---|---|
| 0 | ∞ | 11.0000 | 2.0000 |
| 0.005 (bl/2) | 2.0000 | 1.8333 | **2.0000** |
| 0.01 (bl) | 1.0000 | 1.0000 | 1.0000 |
| 1.0 | 0.0100 | 0.0110 | **0.0100** |

One parameter, one stated meaning, no side effects. Rejected alternatives:
a smaller epsilon — it shrinks the distortion without removing it, and leaves
a number in the code whose size is load-bearing for reasons its name denies;
and clamping only the final result — that is the existing per-reaction cap,
which does not stop a single edge asserting an eleven-fold rise.

**Implementation verified against the measured arm C** (ε=1e-9 + ceiling 2,
which is the ε→0 limit of the old form): 3 differing predictions of 847,
**0 of the 3 converged in both arms**, so the difference is the known
stopping-point artifact and the new form is the same model.

### Consequence for the default

The default ceiling of 11.0 **no longer reproduces old behaviour exactly**,
because the interior of the curve is now undistorted. Arm B measured that
change in isolation at −8 cases. So this is a deliberate behaviour change
requiring validation before merge, not a no-op, and the contract's
"default reproduces today" guarantee applies only to the parameter split, not
to the formula replacement. Nothing merges until the wider validation lands.

---

# Phase 5 — the wider validation (US3)

## R8 — The 742-case set was wrong about BOTH changes, in opposite directions

89 pathways, **23,788 curator cases**, same catalog (`cat92`):

| arm | correct | accuracy | macro-F1 | false_positive_change | propagator_missed |
|---|---|---|---|---|---|
| old formula (ε=1e-3, implicit 11×) | 19,431 | 81.68% | 0.7814 | 1,825 | 709 |
| **epsilon-free, ceiling 11** | **19,456** | **81.79%** | **0.7834** | **1,819** | **690** |
| epsilon-free, ceiling 2 | 19,417 | 81.63% | 0.7807 | 1,826 | 722 |

Two separate verdicts:

**The formula change is POSITIVE at scale: +25 cases, +0.0020 macro-F1**, and
it has the lowest count in *both* relevant failure categories. On the
742-case set the same change measured **−8**. The small set had the sign
wrong.

**Tightening the ceiling to 2 is NEGATIVE at scale: −39 cases against the
matched baseline, −0.0027 macro-F1**, and the failure bucket it was designed
to fix — `false_positive_change` — gets marginally *worse* (1,819 → 1,826)
while `propagator_missed` rises by 32. On the 742-case set it measured
**+6**. The small set had that sign wrong too.

**So the ten-pathway set mispredicted both decisions, in opposite
directions.** This is the clearest evidence yet for R13's warning, and it
retroactively justifies FR-007 having been written into the spec before any
of these numbers existed.

**Decision: the ceiling stays at 11.0. Not adopted, per SC-004.** The
de-repression weakness is real — 78 knockout→UP calls at 0.359 against
MP-BioPath's 0.769 — but a single global ceiling is not its fix. Whatever
the answer is, it is not one number applied to every inhibitor edge.

## R9 — The formula change: the axes disagree, and it is a judgement call

| axis | cases | result |
|---|---|---|
| curator (`cat92`) | 23,788 | **+25 cases, macro-F1 +0.0020** |
| experimental (`cat7`) | 742 | **−5 cases** among both-arms-converged |

The experimental figure is not noise and should not be dismissed as such: of
16 changed predictions, 5 converged in both arms and **all 5 went the wrong
way**. The headline −8 rests on 11 non-converged cases, but the converged
core is a real −5.

SC-004 requires improvement on both, so **this is raised rather than
resolved**. The asymmetry to weigh: the curator evidence is 32× the sample
and spans 89 pathways; the experimental set is 742 cases of which 60% come
from two pathways, and carries a known ~19% subset with no directed path that
no propagator change can move.

**Recommendation — adopt the formula change, reject the ceiling change.** The
epsilon removal is a correctness fix that stands independent of score: a
constant documented as a divide-by-zero guard was sized at 10% of the scale
it guarded and was silently setting the de-repression ceiling, compressing
the interior of the response curve, and shifting maximum suppression. That is
wrong whatever it scores, and at scale it also happens to score better. But
the −5 is the lead's call to accept, not mine to absorb.


---

# R10 — The original suggestion was right: the fix is one line, and my rewrite earned nothing

On the proposed formula change, the design intent was restated: it was needed to avoid
divide by zero and I thought the solution would just be to make it really
really small."*

That is correct, and I should have tested it before building anything. The
blow-up as x → 0 was **never unhandled** — `clamp(result, 0, 10.0)` already
sat at the end of the branch. So a tiny epsilon gives an exact `bl/x` curve
in the interior, and the existing clamp handles the division.

Measured, 89 pathways / 23,788 curator cases:

| arm | correct | accuracy | macro-F1 | false_positive_change |
|---|---|---|---|---|
| current (ε=1e-3) | 19,431 | 81.68% | 0.7814 | 1,825 |
| **ε=1e-12, formula unchanged** | **19,459** | **81.80%** | **0.7835** | **1,814** |
| rewrite, ceiling 11 | 19,456 | 81.79% | 0.7834 | 1,819 |
| rewrite, ceiling 2 | 19,417 | 81.63% | 0.7807 | 1,826 |

**Shrinking epsilon captures the entire gain and then some** — +28 cases,
three better than the rewrite, with the lowest false-positive count of any
arm. The formula restructuring and `DS_DEREPRESSION_MAX` contributed
*nothing*, so both are dropped. A parameter that earns nothing is worse than
no parameter.

## What the feature becomes

One line — `DS_INHIBITOR_EPS` from `0.001` to `1e-12` — plus two comments:
one saying why an epsilon here must stay orders of magnitude below baseline,
and one on `clamp(result, 0, 10.0)` naming it as the de-repression ceiling.
That comment is the knowledge that was actually missing; its absence is why
nobody noticed the identical defect in `DS_HILL_SAT_EPS`.

## The process failure, recorded because it is the useful part

I built a 2×2 attribution experiment **specifically** to avoid mis-attributing
a gain, and then did not run the simplest arm until it was asked for. The
2×2 tested my proposal against itself; it never tested it against doing less.

Every arm in it was a variant of "restructure the formula and add a
parameter". The null hypothesis — "change the constant and touch nothing
else" — was not among them, because by the time I designed the experiment I
had already decided what the fix was. Finding a real defect is not the same
as knowing the smallest thing that fixes it, and an attribution design
inherits the blind spots of whoever chose its arms.

## Still standing

The de-repression weakness is unchanged by any of this: 78 knockout→UP calls
at 0.359 accuracy against MP-BioPath's 0.769 on the same cases. A global
ceiling is falsified as the fix (2× cost 39 cases). The `10.0` clamp remains
an assumption nobody has justified from data — now at least it is labelled.
