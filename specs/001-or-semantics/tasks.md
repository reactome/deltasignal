# Tasks: Faithful OR semantics for set-valued catalysts

**Spec:** [spec.md](./spec.md) · **Plan:** [plan.md](./plan.md)

## Phase 1 — Representation (LNG) — DONE

- [X] **T001** Confirm the defect: `and_or` derived from sign, zero `pos/or`
      regulatory edges in the 92-pathway catalog. *(catalyst/pos/and 33,531;
      regulator/pos/and 3,849; regulator/neg/or 2,409)*
- [X] **T002** Establish scope: verify set-derived catalyst members are
      **co-located on one VR** (not split). *(117 groups, all ≥2 members,
      2–26 each, zero single-member groups)*
- [X] **T003** Add `LNG_SET_MEMBERS_OR` in `append_regulators`, positive edges
      only, EntitySet/DefinedSet/CandidateSet labels only. *(commit 2fdfbda)*
- [X] **T004** Add the flag to the cache fingerprint so it invalidates stale
      caches. *(`_FINGERPRINTED_ENV`)*

## Phase 2 — Processing (DeltaSignal) — DONE

- [X] **T005** Forward `inhibitor_is_and` into `IndexedReaction` and honor it
      (`DS_INHIBITOR_OR`). *(was parsed then dropped)*
- [X] **T006** Add `DS_OR_MODE=capacity` fold-space aggregator with
      `DS_OR_REDUNDANCY` weight.
- [X] **T007** Characterise the aggregator (KO 1-of-N table); identify that
      `w=1.0` leaves N≥8 above the DOWN cutoff.
- [X] **T008** Diagnose that `max(and, or)` discards OR-cluster loss entirely
      (fold 1.0 even at `w=0`); add `DS_OR_COMBINE=gate`.

## Phase 3 — Validation — DONE

- [X] **T009** Regenerate the 3 eligible pathways at Release97, control vs
      treatment; confirm control reproduces the shipped catalog exactly.
- [X] **T010** Benchmark OR + `max` sweep → catastrophic (−0.19 to −0.21).
- [X] **T011** Benchmark OR + `gate` sweep → parity at every `w`.
- [X] **T012** McNemar test every arm → **no arm significantly differs from
      control**; record that `w=0.5` (+2 cases, p=0.688) is not a win.
- [X] **T013** Keep every knob default OFF; no shipped behaviour change.

## Phase 4 — Remaining / blocked

- [ ] **T014** Decide whether to adopt `gate` + `LNG_SET_MEMBERS_OR` as
      defaults. Trade: faithful curator semantics + a logically correct
      combination rule, at statistically neutral accuracy. **Needs Adam's call.**
- [ ] **T015** Re-measure T014 on a wider slice than the 3 MP-BioPath-eligible
      pathways (e.g. the curator ground truth, 72 pathways) before adopting —
      223 cases cannot resolve a ±2-case effect.
- [ ] **T016** Expression-aware set membership: weight alternatives by measured
      expression in the relevant cell line instead of the global `w` fudge.
      Data may already exist (`mp-biopath-pathways/expression/`,
      `gene_to_activity_map.tsv`).
- [ ] **T017** Nested complex-containing-a-set boundary assembly edges still AND
      their set alternatives — needs per-member provenance out of
      `_decompose_regulator_entity`.
- [ ] **T018** `DS_OR_COMBINE=gate` changes semantics for **every** reaction with
      mixed AND/OR activator clusters, not just set-derived ones. Audit that
      blast radius before adopting as default.
- [ ] **T019** Consider whether `max` should be retired entirely as a
      combination rule (it treats a catalyst as a substitute for a substrate).

## Phase 5 — Self-review of this branch (adversarial diff review)

A high-effort adversarial review of **my own diff** raised 15 findings. Fixed:

- [X] **T020** `free_residual` filtered non-finite values out of the max, so a
      NaN/Inf state reported `converged: true` with `null` activities. The code
      it replaced got this right only incidentally (`NaN < tol == false`). Now
      tracked explicitly; non-finite ⇒ never converged.
- [X] **T021** **λ scale mismatch — corrected a wrong claim.** The inner SCC loop
      broke on the *damped* step `λ·|F−x| < tol`, while the new verdict tested
      `|F−x| < tol`. With the default `λ=0.5` a component could legitimately stop
      at up to 2× tolerance and be reported non-converged. This is what produced
      the "122/223 non-converged" figure. With the inner loop measuring the
      undamped residual, it is **0/223**, and the residual is a live measurement
      (2.4e-17 unperturbed, 6.8e-7 for a KO, all under the 1e-6 tolerance).
      **The earlier "55% of solves don't converge" claim was my own bug, not a
      property of the solver.**
- [X] **T022** Dedup kept `is_and` from the first occurrence, making AND/OR
      cluster membership depend on edge order (the repo asserts edge-order
      invariance, and `gate`/`capacity` make membership load-bearing). Now
      AND-wins for `is_and`, union for catalyst, and assembly requires *all*
      duplicates to be assembly so a mixed input+assembly pair isn't
      reclassified into the `assembly_min` rule. Verified order-independent.
- [X] **T023** `isa Real` accepted JSON booleans (`Bool <: Integer <: Real`), so
      `[false, true]` pinned a full knockout at confidence 1.0 and returned 200.
- [X] **T024** For an array-shaped `observations`, `String(node_uuid)` in the
      `@warn` threw before the intended `ArgumentError`, so the client got the
      generic message instead of the specific shape one.
- [X] **T025** Client-supplied `baseline` was unvalidated while being a divisor
      throughout the propagator (`baseline: 0` + `gate` amplifies ~1/ε; negative
      ⇒ NaN). Now required finite and in (0, 1].
- [X] **T026** `DS_OR_REDUNDANCY` was unclamped; `w=2` inverts a knockout into a
      ~2× increase. Clamped to [0,1] at resolve time.
- [X] **T027** Mapping `ErrorException`/`DomainError` to 400 reported genuine
      server faults as client errors (`error(...)` is this codebase's internal
      invariant failure, e.g. `aggregators.jl` parameter-length checks). Narrowed
      to `MethodError`/`InexactError`, and the TSV schema checks now throw
      `ArgumentError` at the read site so they still surface as 400.

Deferred, with reasons:

- [ ] **T028** The final consistency evaluation omits `supply`, so under
      `DS_SCC_BREAK_CATALYST` (default off) `F(x)` is a different operator than
      the one solved and `converged` would always be false. Needs
      `solve_scc_ordered!` to report its own residual rather than a global
      recompute.
- [ ] **T029** Self-loop singletons are invisible to the negative-edge census
      (`comp_size[c] > 1 || continue`), so `DS_SCC_NEG_MODE=transient` (default
      `converge`) treats `A ⊣ A` differently from the equivalent `A ⊣ B → A`.
      Also `DS_SCC_NEG_FRAC=0` would mark every self-loop singleton negative.
- [ ] **T030** `comp_has_self_loop` ignores `activator_break`, so under
      `DS_SCC_BREAK_CATALYST` a deliberately-broken self-catalyst is routed into
      the damped branch and reads its entry value instead of baseline — a
      numeric change smuggled into a reporting fix.
- [ ] **T031** Dedup covers activators only; duplicate inhibitor/depletion edges
      still square their suppression. Low impact in practice — only 3 of 11,901
      catalog duplicate groups are negative.
- [ ] **T032** An observation for a uuid absent from the network is still
      silently dropped (200, unperturbed solve). Pre-existing; the node set is
      available at validation time, so it could be rejected or counted.
- [ ] **T033** A fully-pinned solve (no free nodes) reports `converged: true`
      vacuously. Kept as-is deliberately — "nothing to converge" reads as
      converged, and the empty-network case would otherwise report false.

## Out of scope

- Input/output edge operators — VR splitting already encodes their OR-ness.
- The duplicate `input`+`catalyst` squaring (`DS_DEDUP_ACTIVATORS`) and the
  "bugs partially cancel" hypothesis — separately measured, disproved; see
  the adversarial-review notes in DS PR #12.
