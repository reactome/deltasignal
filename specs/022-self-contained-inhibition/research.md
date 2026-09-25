# 022 — Self-contained inhibition: the input sets the direction, the inhibitor only damps it

## Motivation

specs/012 found that an inhibitor which *contains* its own reaction's input
(a sequestering complex such as WIF1:WNT on *WNT binds to FZD and LRP5/6*)
double-counts that input under `divide` inhibition:
- **one** such inhibitor cancels the input exactly (x · 1/x = 1);
- **two** invert it, giving x · min(x⁻², 10) (specs/021, WNT5A trace).

Deleting those edges was measured at −61 held-out, but that arm was run on a
hash-seed-dependent flag set and is not reproducible (specs/012 addendum).

Adam, 2026-09-25: *"when we have a direct inhibition from an input I feel like
the inhibition should be a lot less than the input. so if the input is up the
output of the reaction should just be a little less up and if it is down
regulated it should be down regulated but just a little less."* And: *"later on
we could learn from data the actual weight it should be."*

## The rule (`DS_SELF_INHIBITOR_WEIGHT=w`, off by default)

An inhibitor I of a reaction is **self-contained** when Reactome's containment
table says I contains the stable id of one of that reaction's activators (the
shared inputs S). Only `divide`-mode inhibitor edges are affected; depletion
edges are not.

I's fold f_I is split into two parts:

- **the part explained by the shared inputs**, f_S = ∏ fold(s ∈ S). This part
  is kept only at a small power w / n, where n is the number of self-contained
  inhibitors on the reaction, so that the reaction's total damping exponent is
  w however many there are;
- **the independent part**, f_I / f_S. This keeps full strength, so a WIF1
  overexpression still inhibits exactly as today.

    effective fold of I = (f_I / f_S) · f_S^(w/n)          (f_S → 0: independent part = 1)

With every self-contained inhibitor tracking its input exactly, the reaction
reads x^(1−w): 2x in gives 2^(0.9) = 1.87x out, and 0.5x in gives 0.54x out.
At baseline the rule changes nothing. w is a single structural parameter. It is
**fixed at 0.1 before the run**, not tuned on the evaluation set, and is the
first candidate for learning from data later.

## Pre-registration (committed before any arm runs)

**Arms.** Both run on build `20260925-1039_d4f4f64`, code defaults otherwise.
- **Control:** `DS_SELF_INHIBITOR_WEIGHT` unset. It must reproduce the
  production scoring (`ae84de9`) byte-for-byte, or the arm is not interpreted.
- **Arm:** w = 0.1.

**Decides:**
- curator **held-out** net cases (fixed − broken) and macro-F1, with an exact
  McNemar p;
- how many pathways moved out of how many were scored, and how many distinct
  readouts the discordant cases span.

**Reported alongside:** tuning split, experimental axis, per-pathway net, and
the WNT5A readout's dose ladder.

**Reading rule, fixed now:**
- **Adopt as the default** only if held-out net > 0 with p < 0.05, the gain is
  not concentrated in one pathway or one readout, and the experimental axis is
  not worse by more than the noise floor (15 cases).
- **Record as neutral** if held-out |net| is within the regeneration noise
  floor (15).
- **Record as negative** if it is below −15.

**Prediction:** a small effect. Only 26 distinct reactions (174 variant nodes)
carry two or more self-contained inhibitors, and 628 pairs have at least one.
Most of them are in PIP3, which is a tuning pathway.
