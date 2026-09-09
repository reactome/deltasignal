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

## Out of scope

- Input/output edge operators — VR splitting already encodes their OR-ness.
- The duplicate `input`+`catalyst` squaring (`DS_DEDUP_ACTIVATORS`) and the
  "bugs partially cancel" hypothesis — separately measured, disproved; see
  the adversarial-review notes in DS PR #12.
