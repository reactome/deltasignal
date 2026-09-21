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

## Queue — REVISED after adversarial review (2026-09-21)

The first draft of this queue was overconfident and in two places misdescribed
this repo's own specs. Verified corrections:

- **Three rows contradicted the spec they cited.** `specs/003:251` says
  "**Keep it reachable** (`DS_SCC_METHOD=minimize`)" because parameter learning
  needs a stable d(prediction)/dθ, and CLAUDE.md still marks 003 **Open**.
  `specs/018:292` says `DS_SCC_BREAK_ROLES` "**stays available**, off by
  default, as the measured record". `specs/019:165` says sink bridges are
  "**kept** as an off-by-default generator flag". A policy may override a
  spec, but the queue presented the specs as agreeing with it. They do not.
- **"Composition: refuted twice" is wrong.** The two numbers are different
  interventions (default `max`, −79; `limit_novel`, −185 on v1), and between
  them `specs/016:1132-1135` records `limit_novel` at **+12 held-out alone and
  +40 with elasticity (p 0.013)**. 016 declined it on *concentration*, not
  harm, and asked for a re-measure on the share-deduplicated catalog that was
  never run.
- **Sink bridges are a null, not a refutation** — p 0.59 and p 0.87. Interferon
  α/β **+100 fixed / 0 broken** is real and unclaimed.
- **Five of nine rows rest on `cat_os` or the v1 fix catalog**, both superseded.
  A mechanism refuted on the welded catalog may have been compensating for the
  weld, and the welds are gone. That bites hardest on the pool family and on
  composition, whose stated failure mechanism is "they re-weld what the fix
  unwelded" — a claim about a cycle census that has since changed twice.

### Executing now (safe, verified)

| action | why safe |
|---|---|
| **Close PR #92** (`LNG_SINK_BRIDGES`) | `git diff main..feat/sink-bridges` is **293 insertions, 0 deletions** — nothing exists on `main` to delete, no shared-code fix to salvage, and the solver has no `sink_bridge` handling. Everything is recorded in specs/019 on `main`. `019:165` amended to match. |
| **Drop `LNG_BOUNDARY_LEAF_REUSE=any`** | the only row with current-catalog, significant, spec-endorsed evidence (specs/018: held-out +173, p<1e-4). No shared helper dies; delete exactly 2 of 6 tests. |

### Blocked, with the specific blocker

| row | blocker |
|---|---|
| `LNG_SHARE_VARIANT_NODES` | **no measurement yet** — specs/020's P2–P4 have not run. The queue stated the outcome as fact. |
| `DS_SCC_METHOD=minimize` | specs/003 says keep it reachable and is Open. Also far larger than described: `SteadyStateParams.mu/.gamma` are **positional fields constructed 38 times across 19 files**, CLI `--mu/--gamma` are **documented in CLAUDE.md's quick-start**, and the provenance keys are asserted in `test_cli_observations.jl` — a *kept* assertion file. |
| `DS_SCC_METHOD=pool` / `pool_parity` | measured on the **welded** catalog, where most SCCs were our own welds (TP53 836→56). Not among the ten arms specs/018 re-ran. Other arms reversed by that much (composition +12→−185, assembly closures +214→−5). One cheap same-catalog arm settles it. |
| `LNG_COMPOSITION_EDGES` / `DS_COMPOSITION_MODE` | see the correction above; and `016:489` required a re-measure on fixed+share that has never been run. |
| `DS_DEPLETION_OWN_PRODUCT` | no p-value in any of its three measurements, and 016 itself says it "was isolated on the wrong base" — the weld that made it untestable is exactly what v2 fixed. |
| `DS_SCC_BREAK_ROLES` | see the plumbing hazard below. |

### Hazards that any future deletion PR must handle

1. **`reaction_model.jl:194` is NOT flag-gated.** `push!(activator_is_assembly, edge.edge_type in ("assembly","composition"))` is the network→reaction mapping, pinned by `test_propagator_invariants.jl:178`. Deleting it as "composition machinery" silently turns every composition edge from a ≤1 cap into a fold multiplier. Catalogs carrying such edges exist on disk (31,159 and 5,625), and `/api/parse` upload accepts any `edge_type`. **Keep it regardless.**
2. **`DS_SCC_BREAK_ROLES` shares plumbing with the kept legacy knob, and yesterday's `d83cfaf` changed the *legacy* knob's behaviour** (a closure activator with no `supply` read the target's baseline; it now reads the live state on the flat path, in the final residual and in influence scores). Deleting by reverting those commits silently reverts the legacy knob too. Delete by **editing the disjuncts** at `steady_state.jl:594/1076/1237`, keep the `elseif stats !== nothing` census (it feeds `cyclic_before/after/largest_after` on the non-roles path), and **port `test_scc_break_roles.jl:105-126` into a new `DS_SCC_BREAK_CATALYST` testset first** — `grep DS_SCC_BREAK_CATALYST test/` returns **nothing**, so that file is the kept knob's only coverage.
3. **`DS_COMPOSITION_GROUP` is in neither list** and reads `activator_group`, which only composition edges populate. Deleting `activator_group` breaks it and ~20 assertions; keeping it leaves a flag whose input nothing can emit. Also `test_propagator_invariants.jl:366-374` tests the **kept** `DS_DEDUP_ACTIVATORS` guarantee *using* `DS_COMPOSITION_MODE=limit` as its probe — rewrite, don't delete.
4. **`test_connectivity_ladder.py:318` inverts.** Its tier-4 branch is conditional on composition edges existing; with the emitter gone it permanently takes the dead-defect path, asserts the inverse, and xfails with a message naming a variable that no longer exists.
5. **Condition 1 as written cannot fail.** Every DS flag being deleted is off by default, so a default-vs-default differential is bit-identical by construction, and the mounted catalog contains none of the affected edge types. Replace with: isomorphism under stable-id relabelling, **plus** an explicit before/after A/B of each *kept* non-default knob (`DS_SCC_BREAK_CATALYST=1`, `DS_COMPOSITION_GROUP=1`, `DS_DEDUP_ACTIVATORS=1`), **plus** a run on a catalog that actually contains the desupported edge types.
6. **Deleted flags become silently ignored.** There is no unknown-`DS_*` guard, so every reproduction command in specs/016–019 would quietly measure the default instead — which contradicts "the spec keeps the evidence". Add the guard in the same PR, or state that those commands are no longer runnable.
7. **Deleting a Julia test file requires editing `.github/workflows/test.yml:55`**, which names each file.
8. **The evidence lives in `/tmp`.** Every catalog backing specs/016–020 is in a session scratchpad; one is already down to 4 of 93 pathways while specs/019 quotes full-catalog numbers from it. "The spec keeps the numbers" needs the numbers to be reproducible.

### Sequencing

1. Close #92; amend `019:165`. 2. Drop `=any`. 3. Finish specs/020 and record P2–P4 before deciding the share row. 4. Re-measure `pool` and `DS_DEPLETION_OWN_PRODUCT` on `cat_fix2` — one cheap same-catalog arm each. 5. Composition: re-measure on fixed+share as `016:489` required, or write down explicitly that the arm is foreclosed — and rewrite the ladder's tier 4 either way. 6. Split `DS_SCC_BREAK_ROLES` from the specs/018 supply fix per hazard 2. 7. Get the catalogs out of `/tmp`; add the flags to `_FINGERPRINTED_ENV`. 8. Add the unknown-`DS_*` guard.
