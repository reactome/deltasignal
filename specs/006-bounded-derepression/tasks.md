# Tasks: Bound de-repression on purpose, not by accident

**Feature**: `specs/006-bounded-derepression` | **Branch**: `006-bounded-derepression`
**Plan**: [plan.md](./plan.md) | **Spec**: [spec.md](./spec.md) | **Contract**: [contracts/inhibition-bounds.md](./contracts/inhibition-bounds.md)

**Tests**: requested. T005's byte-identity check is the one that gates
everything — if the default does not reproduce current behaviour, no
comparison below means anything.

---

## Phase 1: Setup

- [x] T001 Record the pre-change baselines in
      `specs/006-bounded-derepression/research.md` so every arm has something
      to move against: experimental 491/742 macro-F1 0.5771 (`$S/n_sets`),
      curator 3016/4346 macro-F1 0.6711 (`$S/q_cur_base`), and the targeted
      subgroup — knockouts predicted UP, **28 correct of 78**.
- [x] T002 Confirm `$S/cat7` and `$S/cat92` are present and carry
      `node_resolution.csv`, since every arm runs with
      `--resolve-set-readouts` and a missing table changes the denominator
      rather than failing.

---

## Phase 2: Foundational — split the parameter (BLOCKING)

- [x] T003 Add `DS_DEREPRESSION_MAX` to the config in
      `src/core/reaction_model.jl` via `_float_env` with default **11.0**,
      chosen because `(0.01 + 0.001) / 0.001 = 11` is the current implicit
      ceiling, so the default reproduces today's behaviour. A non-numeric or
      negative value must fail loudly, not silently default.
- [x] T004 Apply it in the `"divide"` branch as `min(ceiling, h_k)` **per
      inhibitor, before inhibitors are combined** (before the AND product and
      the OR minimum), per contract guarantee 3. Keep the existing
      `clamp(result, zero(T), T(10.0))` per-reaction cap — it bounds a
      different thing and removing it is not in scope.
- [x] T005 **The gate**: verify the default is behaviour-preserving. Run the
      experimental benchmark with no env overrides on `$S/cat7` and confirm
      the case table is **byte-identical** to `$S/n_sets`. If it is not, stop
      — every later arm is uninterpretable until it is.
- [x] T006 Add assertions to `test/test_config_validation.jl`: the default is
      11.0; a non-numeric value raises; a negative value raises; and the
      value appears in `resolve_reaction_eval_config()`'s result.
- [x] T007 Add an arithmetic test to `test/test_and_curves.jl` (or a sibling
      with real assertions) pinning the bound: with ceiling 2, a knocked-out
      inhibitor raises its target by at most 2× **at one reaction**; with
      ceiling 11 and eps 1e-3 the value matches the current formula to
      floating-point equality.
- [x] T008 Record the value in force in the solve provenance so a result can
      be traced to the assumption that produced it (constitution V).

**Checkpoint**: the two concepts are separable and the default still
reproduces today.

---

## Phase 3: US1 — attribute the measured gain (P1) 🎯 MVP

**Goal**: say which of the three effects ε was moving actually produced the
+15 experimental / +59 curator gain.

**Independent Test**: each arm benchmarked with the other control held fixed,
and the reported cause names the effect responsible.

- [x] T009 [US1] Arm **A** — `DS_INHIBITOR_EPS=1e-3 DS_DEREPRESSION_MAX=11`
      into `$S/r_A`. This is current behaviour and must match T005.
- [x] T010 [US1] Arm **B** — `DS_INHIBITOR_EPS=1e-9 DS_DEREPRESSION_MAX=11`
      into `$S/r_B`. Isolates ε's curve-shape and suppression effects with
      the ceiling held. Note this **strengthens** suppression (0.0110 →
      0.0100), so it moves the opposite way from the original probe.
- [x] T011 [US1] Arm **C** — `DS_INHIBITOR_EPS=1e-9 DS_DEREPRESSION_MAX=2`
      into `$S/r_C`. Isolates the ceiling with ε held.
- [x] T012 [US1] Arm **D** — `DS_INHIBITOR_EPS=1e-2` with the ceiling
      effectively inactive, into `$S/r_D`. Reproduces the original probe for
      continuity; expect ≈506/742.
- [x] T013 [US1] Repeat T009–T012 against `--ground-truth curator`. One arm
      at a time on a fresh port, killing the Julia server between arms.
- [x] T014 [US1] Report all arms through
      `bench/analysis/compare_set_rules.py` — macro-F1, coverage delta,
      per-pathway net, and the both-arms-converged count. **Read coverage
      first**; read the both-converged count before believing any small
      delta, since only 9 of 29 changed predictions converged in both arms in
      the original probe.
- [x] T015 [US1] Report the targeted subgroup separately for every arm
      (FR-009): knockouts predicted UP, showing **both** wrong calls removed
      and correct ones lost, against the baseline of 28 correct of 78. An arm
      that removes all 78 would look fine on the total and be wrong.
- [x] T016 [US1] State the attribution in `research.md`. If B is worse than A
      and C recovers it, the gain is the **ceiling**. If B is better than A,
      part of the gain is **suppression** and this feature's framing is wrong
      — say so plainly per FR-010 rather than reframing. If both move, give
      the proportion and do not pick the flattering half.

**Checkpoint**: SC-001 met. The reported cause matches the actual cause.

---

## Phase 4: US2 — the guard is only a guard (P1)

- [x] T017 [US2] With `DS_DEREPRESSION_MAX=2` fixed, run
      `DS_INHIBITOR_EPS` at 1e-6, 1e-7, 1e-8 and 1e-9 into `$S/r_eps<v>`.
- [x] T018 [US2] Assert the four case tables are **identical** (FR-002). Any
      difference means ε is still doing modelling work and the separation is
      nominal rather than real.
- [x] T019 [US2] Document in `contracts/inhibition-bounds.md` the ε value
      below which only division-by-zero protection remains, so the next
      reader knows what the guard is for.

**Checkpoint**: SC-002 met.

---

## Phase 5: US3 — decide the default on the wider sample (P2)

- [x] T020 [US3] Run the surviving candidate ceilings on **`$S/cat92`**
      against both ground truths. This is the primary evidence; the 742-case
      set is secondary (FR-007), because 60% of it is two pathways and the
      candidate came from sweeping it.
- [ ] T021 [US3] Measure the bound against the classification cutoffs
      (SC-006): re-score each candidate at 0.80/1.20, 0.85/1.15 and
      0.90/1.10. A gain that disappears under a small cutoff shift is a
      threshold artifact, not a modelling win.
- [x] T022 [US3] Adopt a default only if it improves macro-F1 on the
      92-pathway catalog against **both** ground truths. A candidate that
      improves the small set and not the large one is **not** adopted and the
      discrepancy is recorded.
- [ ] T023 [US3] If a default changes, update `CLAUDE.md`'s Solver
      Configuration — current state only, one line, pointing here for the
      rationale.

---

## Phase 6: Polish & Cross-Cutting

- [x] T024 [P] Record every arm in `research.md` including the neutral and
      negative ones, per constitution III.
- [ ] T025 [P] Update `quickstart.md` with the commands and the numbers
      actually observed.
- [ ] T026 File the follow-on rather than absorbing it: `DS_DEPLETION_H_MAX`
      (default 10) carries the identical unexamined assumption on the
      depletion path. It was left alone deliberately so it could not confound
      this attribution; it deserves the same scrutiny separately.
- [ ] T027 Re-check the accuracy decomposition in
      `specs/004-loop-handling/research.md` R12/R15 against whatever ships,
      since those numbers were measured under the 11× ceiling.

---

## Dependencies

```
Phase 1 (T001-T002)
      ↓
Phase 2 (T003-T008)   ← BLOCKING; T005 gates everything
      ↓
   ┌──┴──────────────┐
   ↓                 ↓
Phase 3 US1       Phase 4 US2     (same edit, different evidence)
(T009-T016)       (T017-T019)
   └──┬──────────────┘
      ↓
Phase 5 US3 (T020-T023)
      ↓
Phase 6 (T024-T027)
```

- **T005 gates every arm.** A default that does not reproduce current
  behaviour makes all four arms uninterpretable.
- **T016 gates Phase 5.** Choosing a ceiling before knowing whether the
  ceiling is what mattered would be choosing blind.
- T017–T019 need only T003–T004, so Phase 4 can run alongside Phase 3.

## Parallel opportunities

- Phases 3 and 4 are independent bodies of evidence on the same edit.
- T024 and T025 are different files.
- **Benchmark arms are never parallel.** One at a time, `pkill` the server
  between them; concurrent servers have caused an out-of-memory crash here.

## Implementation strategy

**MVP is Phase 2 + Phase 3.** The parameter split is worth landing on its own
— it turns an accidental modelling assertion into a stated one even if the
default never changes — and Phase 3 answers the question the feature exists
to ask.

**Phase 5 may conclude "change nothing".** That is a legitimate outcome: if
the 92-pathway catalog does not support a smaller ceiling, the right result
is an explicit 11× with the evidence recorded, which is still strictly better
than an implicit 11× nobody chose.

## Notes on discipline for this feature

- **The attractive wrong move is to skip Phase 3** and ship the ceiling
  because the probe looked good. The probe moved three things at once.
- **The ceiling is a free parameter** on a project with a documented history
  of fitting them to this case set. FR-007 puts adoption on the larger
  catalog; do not sweep finely on the 742.
- **A result that overturns the framing is still a result.** If the gain was
  suppression, R15 and R16 in `specs/004-loop-handling/research.md` need
  correcting, and that correction is the deliverable.
