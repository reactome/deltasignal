# Flag cleanup: a flag is a measurement instrument with an expiry

**Decided with Adam, 2026-09-21.** The generator carries 19 `LNG_*` variables
and the solver 47 `DS_*`; most are answered questions. Constitution IV says new
behaviour arrives behind a flag "**until it is measured**" — the second half was
never honoured, so measurement scaffolding accumulated as permanent surface.

## Policy

A flag exists only while a question is open. On resolution it expires:

- **Adopted** → delete the flag; the behaviour becomes the code; the spec keeps
  the evidence.
- **Refuted** → delete the flag **and the code**; the spec keeps the numbers so
  the question is not re-opened. Git holds the history; a spec that says
  "measured, −185 held-out, here is the mechanism" is a better deterrent than
  dead code behind a switch.

Permanent knobs are only: genuine model parameters a user might legitimately
tune (epsilons, Hill exponents, damping, caps) and operational config (paths,
threads, debug, determinism).

## Queue (one reviewable PR per repo, after the specs/020 measurement lands)

| flag | status | action |
|---|---|---|
| `LNG_SHARE_VARIANT_NODES` | specs/020 (running): better design, size −35%, reachability-neutral | **becomes the behaviour**, flag removed |
| `LNG_BOUNDARY_LEAF_REUSE=any` | specs/018: fix adopted default-on | drop the escape hatch |
| `LNG_COMPOSITION_EDGES` | specs/016 and 018: refuted twice (−79, −185 held-out) | delete code + flag |
| `LNG_SINK_BRIDGES`, `LNG_SINK_BRIDGE_MAX_FANOUT` | specs/019: refuted (−11 held-out, false change +195) | delete code + flag |
| `DS_SCC_METHOD=pool\|pool_parity\|pool_all` | specs/017: refuted on both axes | delete; `DS_SCC_METHOD` returns to two values |
| `DS_SCC_METHOD=minimize` | specs/003: loses at every γ | delete (with `ComponentProblem`, `lm_minimize!`, `OverlayVector`, `DS_MU`/`DS_GAMMA*`/`DS_LM_ITERS`/`DS_SCC_OPTIMIZER`) |
| `DS_COMPOSITION_MODE` | specs/016: refuted | delete |
| `DS_DEPLETION_OWN_PRODUCT` | specs/016: inert | delete |
| `DS_SCC_BREAK_ROLES` | specs/018: superseded by the generator fix (≈0 on the fixed catalog) | delete |

Kept deliberately: `DS_LOOP_ELASTICITY*` (one measurement, +29/−85, still an
open question), `DS_SCC_SWEEP=jacobi` (the label-independence tool, and
specs/019 C3 shows why we still need it), `DS_SKIP_EDGE_TYPES` (the A/B
instrument itself), and every genuine parameter.

## Conditions on the cleanup PRs

1. Deleting code is where a silent behaviour change hides most easily, so each
   PR must show a **differential run**: production defaults before and after,
   bit-identical on the catalog.
2. Every assertion file must still pass, and any test that only existed to
   cover a deleted flag is deleted with it — tests for *kept* behaviour stay.
3. Adversarial review of each PR before merge, pointed specifically at
   behaviour change hidden in deletion.
4. `CLAUDE.md`'s knob list and the specs' "current state" sections updated in
   the same PR.
