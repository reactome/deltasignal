# Tasks: Every node maps, both ways, including on the diagram

**Feature**: `specs/005-node-identity-mapping` | **Branch**: `005-node-identity-mapping`
**Plan**: [plan.md](./plan.md) | **Spec**: [spec.md](./spec.md) | **Contract**: [contracts/node-resolution.md](./contracts/node-resolution.md)

**Tests**: requested, and one of them is load-bearing. T010's negative
control is not polish — a completeness check that cannot fail is the exact
defect this project has shipped (`_decomposed_ids` passed 11 of 11 while
masking 18 dropped catalysts).

**Repos**: tasks touch both. `logic-network-generator` (LNG) produces the
mapping; `deltasignal` (DS) consumes it.

---

## Phase 1: Setup

- [x] T001 Record the pre-change baseline in
      `specs/005-node-identity-mapping/research.md` so every later claim has
      something to move against: orphan rate per role from
      `node_reaction_context.csv` (input 3,641/0, output 3,491/0, catalyst
      2,228/2,228, regulator 1,018/1,018), the glyph join rate (114 of 156,
      73.1%), and the benchmark's discarded-case count (204 of 847).
- [x] T002 Confirm `LNG_DIAGRAM_DIR` resolves to the Release97 diagrams and
      that `R-HSA-1257604.json` and `.graph.json` are both present — the
      glyph work silently degrades to nothing if only one is there.

---

## Phase 2: Foundational — LNG #67 (BLOCKING)

**⚠️ Nothing downstream measures correctly until this lands.** Catalysts and
regulators are 31.3% of the context export, they are where the causality
is, and 25 of the 42 glyph-join misses are this same bug.

- [x] T003 In `logic-network-generator/src/logic_network_generator.py`, make
      `append_regulators` record every `(member_uuid, reaction_id, role)`
      triple it actually emits, at the point it appends the edge (~line
      1590), where `member_uuid` is the decomposed terminal member and role
      is `catalyst` or `regulator` per `edge_type`. Return it alongside the
      existing outputs rather than reconstructing it later.
- [x] T004 Change `export_node_reaction_context` in the same file to consume
      those recorded triples for the catalyst/regulator branch instead of
      `r.get("uuid")` from `catalyst_map` / `negative_regulator_map` /
      `positive_regulator_map`. Those fetch rows carry the **undecomposed
      parent** entity's uuid, which is why 100% of both roles orphan.
- [x] T005 Add an assertion to `export_node_reaction_context` that every
      emitted `context_node` appears in the logic network, and fail the
      export rather than writing an orphaned row. Per the contract, "never
      emit a row whose uuid is absent from `logic_network.csv`".
- [x] T006 Add a regression test in `logic-network-generator/tests/`
      asserting the orphan rate is 0 for **all four** roles. It must fail
      against the current code — run it before T003 to prove it can.
- [x] T007 Regenerate the ten-pathway catalog into a fresh directory and
      re-measure the orphan rate. **`rm -rf <dir>/*/cache` first**:
      regenerating into a populated directory adopts pre-fingerprint caches
      and silently reproduces the old output while reporting success.

**Checkpoint — DONE.** Measured before and after on the ten-pathway
catalog with `bench/analysis/resolution_audit.py`:

| role | rows before | orphaned before | rows after | orphaned after |
|---|---|---|---|---|
| input | 3,641 | 0 | 3,641 | 0 |
| output | 3,491 | 0 | 3,491 | 0 |
| catalyst | 2,228 | **2,228 (100%)** | 571 | **0** |
| regulator | 1,018 | **1,018 (100%)** | 219 | **0** |
| TOTAL | 10,378 | 3,246 (31.3%) | 7,922 | **0 (0.0%)** |

**Zero orphans is trivially achievable by dropping rows, so coverage was
checked separately**: all **790 of 790** distinct resolvable
`(node, reaction, role)` triples now appear in the export, against **0 of
790** before. The row-count drop is not loss — 3,758 catalyst edges collapse
to 790 distinct triples because the same member catalyses the same reaction
through many parallel virtual-reaction instances, and the old 2,228 was
counting a different and wrong population (parent × reaction).

**The fix is export-only**: all ten regenerated networks are structurally
identical to the originals (edge counts, node counts and the
edge_type/sign/and_or/stoichiometry signature all match 10/10), so no
benchmark number moves because of it.

**The bug was held in place by a test.** `test_export_node_reaction_context`
asserted the parent uuid appears in the export — the exact wrong behaviour —
so the suite was green while a third of the export was broken. Inverted and
kept as the record. Negative control run: reintroducing the parent uuid
makes the new test fail; restoring the fix makes it pass. 167 tests pass in
the no-database tier.

Landed on `logic-network-generator` branch `fix/67-context-export-orphans`.

**Still open from this phase**: catalyst and regulator edges carry no
`edge_reaction_id` at all (3,758 and 1,124 of them). The export now recovers
the reaction through the target VR node instead, so nothing is blocked, but
the column is empty where it should not be. Tracked, not fixed here.

---

## Phase 3: US1 + US2 — the mapping and its invariant (P1) 🎯 MVP

**Goal**: every entity resolves to nodes and back, with set readouts
resolving as a consequence.

**Independent Test**: a completeness report over the ten pathways shows zero
unexplained absences in either direction, and a set readout scores.

- [x] T008 [US2] Create `logic-network-generator/src/set_resolution.py`
      resolving an EntitySet to its leaf members through
      `hasMember|hasCandidate`, recursing through nested sets. Carry `depth`
      (0 for self, incrementing per hop; **observed maximum is 5**). Keep a
      visited-set guard — **note honestly in the docstring that Release97
      has no membership cycles** (0 self-membership, 0 at lengths 2/3/4
      across 5,440 nested sets) and the guard is insurance against future
      curation error, not a response to an observed cycle.
> **T008 result (measured against all 20 real blocked readouts, not
> fixtures):** every one resolves — 18 fully, 2 partially, **0
> unresolvable**, against 3 unresolvable at one hop. Depth distribution 15×1,
> 4×2, 1×3, so recursion is load-bearing for 5 of 20.
>
> **New defect found, and it changes T010's expectation.** The 2 partial sets
> (ERBB2 `R-HSA-1963585`, `R-HSA-1963588`, 10 cases) are partial because some
> of their leaf Complexes have no node. Those leaves are **not** direct
> reaction inputs or outputs — each reaches its reactions only *via an
> EntitySet participant* (9–12 reactions each) — and structurally identical
> siblings differ: `R-HSA-1963593` and `R-HSA-1248703` got nodes while
> `R-HSA-1963583` and `R-HSA-1250316` did not. So **set expansion drops some
> members**, which is a generation bug distinct from this mapping.
>
> This does **not** falsify research.md R5, which measured *direct* input and
> output participants and found every gap to be an EntitySet. It adds a layer
> R5 did not look at. `node_exclusions.csv` will therefore **not** be empty,
> and T010 should treat a non-empty list as the finding it is rather than
> forcing it to zero.

- [ ] T009 [US2] Emit `node_resolution.csv` per pathway from
      `logic_network_generator.py` with exactly the columns in
      `data-model.md`: `stable_id, uuid, relation, depth, role,
      reaction_stid, glyph_id, diagram_stid, release`. `relation` is one of
      `self`, `set_member`, `complex_component`, `variant`, `reaction`,
      `boundary_collapsed`, `dissociation_sink` — **there is deliberately no
      catch-all value**; a node whose relation cannot be determined is a
      defect to surface, not an `other` bucket. `stable_id` is primary; a
      database id may be a separate column but never the key.
- [ ] T010 [US2] Emit `node_exclusions.csv` with `stable_id, reason,
      release`, `reason` required and non-empty. **Expect it empty**: every
      unresolved entity measured is an EntitySet (29 of 29 in PIP3, 9 of 9
      in Cell Cycle Checkpoints), so a non-empty list is a finding to
      investigate rather than a category to accept.
- [ ] T011 [US2] Add the bidirectional completeness check to
      `logic-network-generator/scripts/validate_logic_network.py`: forward
      (every entity participating as input or output of any reaction appears
      in the resolution table or the exclusion list with a reason), reverse
      (every uuid in `logic_network.csv` appears in the resolution table),
      and referential (every uuid in the resolution table exists in the
      network).
- [ ] T012 [US2] Add `logic-network-generator/tests/test_resolution_negative_control.py`
      building a fixture with rows deliberately removed and asserting the
      T011 checks **FAIL**. Run it against a correct mapping too, to prove
      it passes when it should. This is the task that stops T011 becoming
      another check that cannot fail.
- [ ] T013 [US1] In `deltasignal/bench/benchmark_mpbiopath_cases.py`,
      resolve a set-valued readout through `node_resolution.csv` to its
      `set_member` leaves instead of returning `absent_from_network` or
      `proxy_available_not_enabled`.
- [ ] T014 [US1] Implement the two-level combine: collapse the positional
      uuids of the **same member first**, then combine across members.
      Without the first level a member split into three positions outweighs
      one split into two purely by decomposition (`R-HSA-202074`'s members
      map to 3, 2 and 2 uuids).
- [ ] T015 [US1] Report a partially resolved set as **partial**, naming the
      missing members, and do not combine over the members that happened to
      resolve (FR-009).
- [ ] T016 [US1] Refuse to use a resolution table whose `release` differs
      from the networks' without an explicit override (FR-005). Version skew
      has produced a false finding on this project already.
- [ ] T017 [US1] Re-run the benchmark and report the discarded-case count.
      Expect **204 → ≤27**. Report the count alongside any accuracy figure
      (FR-012), because it sets the denominator.
- [ ] T018 [US2] Populate or remove `source_sets`, `chosen_members` and
      `compartment` in `nodes.csv` (FR-011). All three are empty in every
      row of every pathway and `compartment` is hardcoded `None`; a declared
      column that is always empty misleads its reader.

**Checkpoint**: SC-001, SC-002, SC-003, SC-006, SC-008 met.

> **The denominator moved.** ~741 scored cases, not 564. No accuracy figure
> from here is comparable with 365/564 or with anything published earlier
> unless both denominators are stated.

---

## Phase 4: US3 — glyph identity (P2)

**Goal**: clicking a glyph resolves to the right nodes, and a solved node
resolves back to the glyph to highlight.

- [ ] T019 [US3] In `logic-network-generator/src/diagram_connectivity.py`,
      stop discarding glyph ids: expose `(reaction_stid, entity_stable_id,
      role) -> glyph_id` from the layout, alongside the producer/consumer
      pairs it already computes.
- [ ] T020 [US3] Populate `glyph_id` and `diagram_stid` in
      `node_resolution.csv` by joining on that triple. **Never write one
      without the other** — a glyph id is unique only within its diagram.
- [ ] T021 [US3] Add a validator check that every `glyph_id` carries a
      `diagram_stid` and resolves within that diagram.
- [ ] T022 [US3] Create `deltasignal/bench/analysis/glyph_join_report.py`
      reporting join coverage: joinable triples, diagram-only (with role
      breakdown), and LNG-only. **LNG-only is expected** — generation
      descends below the diagram, 411 triples on PIP3 — and must be reported
      as expected rather than as failures.
- [ ] T023 [US3] Assert the duplicate case explicitly: ATP has 9 glyphs in
      R-HSA-1257604 and each must resolve distinctly. **0 of 156 diagram
      triples resolved to more than one glyph**, so any ambiguity is a
      regression, not a tolerance.
- [ ] T024 [US3] Where several glyphs genuinely collapse to one node (ADP: 9
      glyphs → 1 node via the boundary cache), record the collapse as
      `boundary_collapsed` rather than letting the mapping look lossy.

**Checkpoint**: SC-004, SC-005 met.

---

## Phase 5: US4 — how a set's members combine (P2)

- [ ] T025 [US4] Implement `mean`, `max` and `mean_reachable` as selectable
      level-2 rules in `deltasignal/bench/benchmark_mpbiopath_cases.py`,
      **defaulting to current behaviour** until measured (FR-010).
- [ ] T026 [US4] Benchmark each arm separately on the same catalog, one at a
      time, killing the Julia server between arms. Report macro-F1,
      per-pathway net change, both-arms-converged count and coverage delta.
- [ ] T027 [US4] Record the arithmetic for whichever rule is adopted, not
      only its score: with uniform baseline x₀, `Σxᵢ/(n·x₀) = mean(xᵢ)/x₀`,
      so `mean` is the sum-of-abundances reading and matches Adam's stated
      rule that OR configurations average. Record the losing arms too
      (constitution III).
- [ ] T028 [US4] Check the dilution failure explicitly: report how often
      `mean` and `mean_reachable` differ, and whether the difference is
      concentrated in TP53, where signal dilution is a recorded failure.

---

## Phase 6: Polish & Cross-Cutting

- [ ] T029 [P] Update `specs/005-node-identity-mapping/research.md` with
      every measured outcome including the negative ones.
- [ ] T030 [P] Update `specs/005-node-identity-mapping/quickstart.md` with
      the commands and numbers actually observed.
- [ ] T031 Re-baseline feature 004's remaining loop arms against the
      corrected case set. Its control result (403 / 405 / 364) stands, being
      measured on a case set scored in all three arms, but the intervention
      arms must not be compared to 365/564.
- [ ] T032 Add the resolution artifacts to `logic-network-generator/CLAUDE.md`
      as current state, one line, pointing here for rationale. Close LNG #67
      with the measured before/after.

---

## Dependencies

```
Phase 1 (T001-T002)
      ↓
Phase 2 #67 (T003-T007)   ← BLOCKING: everything downstream mismeasures until this lands
      ↓
Phase 3 US1+US2 (T008-T018)   ← MVP
      ↓
   ┌──┴──────────────┐
   ↓                 ↓
Phase 4 US3      Phase 5 US4
(T019-T024)      (T025-T028)
   └──┬──────────────┘
      ↓
Phase 6 (T029-T032)
```

- **T006 before T003.** The regression test must be shown to fail first.
- **T012 gates T011 being believed.** A completeness check without a
  negative control is not evidence.
- **T008 gates T009**; the resolution table cannot be written before set
  closure exists.
- **T019 gates T020.**
- Phases 4 and 5 are independent of each other and both depend on Phase 3.

## Parallel opportunities

- T009 and T010 are the same emitter and are **not** parallel.
- T019–T021 (LNG) and T022 (DS) are different repos once the column exists.
- T029 and T030 are different files.
- **Benchmark arms are never parallel.** One at a time, kill the server
  between them; concurrent servers have caused an out-of-memory crash here.

## Implementation strategy

**MVP is Phase 2 + Phase 3.** #67 has to land regardless — it is a
correctness bug in shipped output — and Phase 3 delivers Adam's actual
request plus 177 recovered cases. Phases 4 and 5 are the pathway-browser
deliverable and the accuracy question respectively, and both are better done
against a mapping that already exists and is verified.

**Phase 2 is worth landing on its own** even if the rest slips: a 31.3%
orphaned export is wrong output that other work will keep tripping over.

## Notes on discipline for this feature

- **The denominator change will masquerade as an accuracy change.** Going
  from 564 to ~741 scored cases moves every rate. State both denominators or
  neither, every time.
- **An empty exclusion list is the target, not a formality.** R5 measured
  every gap as an EntitySet. If something else shows up in the other eight
  pathways, that is the interesting finding of this feature.
- **The cycle claim was wrong once already.** research.md R1 retracts it.
  Do not reintroduce "Reactome sets can be cyclic" as a justification
  anywhere; the honest reason for the guard is future-proofing.
