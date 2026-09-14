# Contract: inhibition bounds

Two numeric controls, previously one.

## Before

| control | stated purpose | actual effect |
|---|---|---|
| `DS_INHIBITOR_EPS` (1e-3) | prevent division by zero | that, **plus** the de-repression ceiling (11×), **plus** the curve shape at intermediate values, **plus** the suppression floor |
| `DS_INHIBITOR_FLOOR` (0.0) | bound suppression | as stated, scoped by `DS_INHIBITOR_FLOOR_SCOPE` |
| hard-coded `10.0` | bound per-reaction de-repression | as stated, not configurable |

## After

| control | purpose | default |
|---|---|---|
| `DS_INHIBITOR_EPS` | prevent division by zero, and nothing else | unchanged at 1e-3 until US3 decides |
| **`DS_DEREPRESSION_MAX`** | **the most a single inhibitor's removal may raise its target** | **11.0 — reproduces current behaviour** |
| `DS_INHIBITOR_FLOOR` | unchanged | unchanged |

Named after the depletion path's existing `DS_DEPLETION_H_MAX`, which is the
same idea already present in the same function.

## Guarantees

1. **Default reproduces today.** With `DS_DEREPRESSION_MAX=11.0` and
   `DS_INHIBITOR_EPS=1e-3`, predictions are unchanged.
2. **The guard is only a guard.** With the ceiling fixed, varying
   `DS_INHIBITOR_EPS` across three orders of magnitude changes no prediction
   except where a division by zero would otherwise occur (FR-002).
3. **The bound is per inhibitor**, applied before inhibitors are combined, so
   the claim it encodes is about one edge and is checkable in isolation. The
   existing per-reaction cap stays.
4. **Compounding across reactions is still possible**, and that is a stated
   limitation rather than an oversight — a chain of genuinely de-repressed
   steps should compound. Only the per-step magnitude is bounded.
5. **Validation is loud.** A non-numeric or negative value is a configuration
   error, not a silent default.
6. **Provenance.** The value in force is recorded in the solve output, so a
   result can be traced to the assumption that produced it.

## What this does NOT change

- The inhibition mode (`divide`) — settled, out of scope.
- `DS_DEPLETION_H_MAX` — the same assumption on the depletion path, left
  alone so it cannot confound this feature's attribution.
- The suppression direction, except as a measured consequence of ε that
  User Story 1 exists to quantify.
