# Tasks: Handle loops properly, and compare honestly

**Feature**: `specs/004-loop-handling` | **Branch**: `004-loop-handling`
**Plan**: [plan.md](./plan.md) | **Spec**: [spec.md](./spec.md)

**Tests**: assertion-bearing tests requested only where a classification rule
ships (T009). This feature is measurement; most tasks produce numbers, and a
number is validated by the A/B protocol, not by a unit test.

---

## Phase 1: Setup

- [x] T001 Confirm the shared catalog at `$S/cat` is intact and unchanged
      (10 pathway directories, LNG main 4ff0408, Release97) and record its
      per-pathway `logic_network.csv` line counts in
      `specs/004-loop-handling/research.md` so later arms can prove they read
      the same build.
- [x] T002 Confirm the baseline experimental arm `$S/i_defaults` and curator
      arm `$S/j_cur_new` are present and still score 365/564 and 2628/3914;
      every arm in this feature is compared against those two directories.

---

## Phase 2: Foundational (blocking)

**⚠️ Nothing in Phase 3+ may be reported before T003 and T004 exist**, because
every arm's headline requires the cyclicity numbers and the coverage guard.

- [x] T003 Create `bench/analysis/cycle_structure.py` implementing the
      contract in `specs/004-loop-handling/contracts/cycle-structure-report.md`:
      `--catalog` or `--mpbiopath` (exactly one; both or neither is an error,
      not a default), optional `--classify`, `--ratio-threshold` defaulting to
      15, optional `--json`. Per pathway report node count, edge count,
      self-loop count, components larger than one node, total cycle-resident
      nodes and largest component. Read-only; exit non-zero on a malformed
      network rather than reporting partial results as complete.
- [x] T004 Add to `bench/analysis/cycle_structure.py` the `--classify` output:
      per component the intra-edge count, negative intra-edge count, distinct
      `edge_reaction_id` count, `nodes_per_reaction`, and the class
      (`recycling_artifact` when `nodes_per_reaction >= threshold`, else
      `candidate_feedback`). **The threshold in force must be printed even
      when it is the default** — a classification whose parameter is implicit
      is not reproducible.
- [ ] T005 Add a coverage-delta helper to `bench/analysis/` that, given two
      benchmark output directories, reports cases scoreable in the baseline
      but not the arm and the reverse, plus macro-F1, per-class F1,
      per-pathway net change and the count of changed predictions where both
      arms converged. Every later task reports through this one helper so no
      arm can be quoted with a missing field.

**Checkpoint**: cyclicity is measurable and no arm can be scored without its
coverage delta.

---

## Phase 3: User Story 1 — Attribute the deficit (Priority: P1) 🎯 MVP

**Goal**: split the 42-case deficit into a propagator share and a network
share by running the DeltaSignal propagator on MP-BioPath's own acyclic
networks.

**Independent Test**: a third column appears next to DeltaSignal-on-LNG (365)
and MP-BioPath's published predictions (407), on one common case set.

- [x] T006 [US1] Create `bench/analysis/mpbiopath_network_adapter.py`
      converting each `~/gitroot/mp-biopath-pathways/pathways/*.tsv` into a
      DeltaSignal catalog directory. Columns are parent dbid, child dbid,
      polarity (`1` = pos, `-1` = neg), conjunction (**`0` = AND, `1` = OR**
      per research.md R4). Emit `logic_network.csv` with columns
      `source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id`
      and `stid_to_uuid_mapping.csv` with `uuid,stable_id` — those two files
      are the server's stated minimum (`src/api/server.jl:423-425`);
      `nodes.csv` and `entity_reaction_proxy_mapping.csv` are optional and
      should be omitted rather than faked.
- [x] T007 [US1] In the adapter, use the MP-BioPath dbid verbatim as the node
      uuid and emit `stable_id` as `R-HSA-<dbid>`, so that
      `load_dbid_to_uuids`'s `stable_id.rsplit("-", 1)[-1]` recovers exactly
      the dbid the benchmark cases key on. No translation layer is needed
      (research.md R2); do not add one.
- [x] T008 [US1] Name each output directory `<PathwayName>_R-HSA-<id>` to
      match `find_pathway_dir`'s suffix matching, taking the id from the
      corresponding `$S/cat` directory name. Report any MP-BioPath pathway
      with no catalog counterpart, and any catalog pathway with no
      MP-BioPath file, rather than silently dropping either.
- [x] T009 [US1] Add assertions for the adapter to a test file that already
      contains assertions or a new one, covering the conjunction mapping
      (`0`→`and`, `1`→`or`), the polarity mapping (`1`→`pos`, `-1`→`neg`),
      and the dbid round-trip through `rsplit("-", 1)[-1]`. A silent
      mis-mapping here would be indistinguishable from a propagator finding.
- [x] T010 [US1] **Validate the adapter before trusting the control**: run
      `bench/analysis/cycle_structure.py --catalog` over the adapted catalog
      and confirm it reproduces the MP-BioPath acyclicity figures from
      research.md R1 (zero self-loops; largest component 11 nodes). A
      conversion that introduced cycles would measure the adapter.
- [x] T011 [US1] Run the benchmark on the adapted catalog into `$S/k_control`
      on a fresh port with default solver config, then kill the server.
- [x] T012 [US1] Report DeltaSignal-on-MPB-networks, DeltaSignal-on-LNG
      (`$S/i_defaults`) and `mpbiopath_prediction` on the **common case set** —
      cases mappable in all three arms — excluding any case unmappable in any
      arm from all three, per FR-002.
- [x] T013 [US1] State the attribution explicitly in `research.md`: how much
      of the 42-case deficit is propagator and how much is network. If the
      control lands near 365, say plainly that the networks are exonerated
      and the propagator is the problem, and hand that to feature 003 rather
      than absorbing it here.

**Checkpoint**: SC-001 met. The answer redirects Phase 5 either way.

---

## Phase 4: User Story 2 — Report honestly (Priority: P1)

**Goal**: no DeltaSignal-versus-MP-BioPath figure is stated without the
acyclicity asymmetry beside it. Independent of every experiment.

**Independent Test**: the case table carries a cyclicity classification and
the split accuracies are reported.

- [ ] T014 [P] [US2] Add to `bench/benchmark_mpbiopath_cases.py` a per-case
      `readout_in_scc` and `perturbation_in_scc` column, computed from the
      solved network's strongly connected components, per data-model.md.
- [ ] T015 [P] [US2] Add `largest_scc_on_path` — the largest component
      containing any node on a shortest signed path from perturbation to
      readout, `0` when there is no path. A case counts as **cyclic** when
      `readout_in_scc` or `largest_scc_on_path > 0`; pathway-level cyclicity
      is not a substitute, because a pathway can be 36% cyclic while a given
      case never touches a cycle.
- [ ] T016 [US2] Report DeltaSignal and MP-BioPath accuracy and macro-F1
      separately on cyclic and acyclic cases (FR-003), on both ground truths.
- [ ] T017 [US2] Correct the record: amend
      `specs/002-upregulation-propagation/research.md` R8 to state that its
      365-versus-407 and 365-versus-393 comparisons are against
      hand-acyclicised networks, with the largest-SCC figures for both sets
      (FR-004). Constitution principle II requires retracting loudly; this is
      a qualification of a published number, so it goes in the spec that
      published it, not only here.
- [ ] T018 [US2] Record the second asymmetry alongside acyclicity: only
      **449 of 25,196 (1.8%)** MP-BioPath edges are negative. Our networks
      encode far more inhibition, and any claim about inhibition handling
      that compares the two must say so.

**Checkpoint**: SC-002 and SC-003 met regardless of what Phase 5 finds.

---

## Phase 5: User Story 3 — Test the loop hypothesis (Priority: P2)

**Goal**: intervene, in increasing order of how much curated structure each
intervention destroys, stopping at the first that both works and is
defensible under constitution principle I.

**Independent Test**: each arm reported with per-pathway breakdown and
both-arms-converged count, including the arms that lose.

- [ ] T019 [US3] Create `bench/analysis/loop_interventions.py` producing a
      derived catalog per arm. It must **never mutate the baseline catalog**;
      it writes a new directory and records `edges_removed`,
      `cycle_nodes_before`, `cycle_nodes_after` and
      `destroys_curated_causality` (true if it removes `catalyst`,
      `regulator`, `depletion`, `input` or `output` edges) into a manifest.
- [ ] T020 [US3] Arm `drop_diagram_bridge`: remove all 5,358 `diagram_bridge`
      edges. This is the only arm that could ship — it is generated
      connectivity, not a curator assertion. Expect cycle-resident nodes
      4,561 → 2,002 and no change to TP53.
- [ ] T021 [US3] Benchmark `drop_diagram_bridge` against both ground truths
      and report through the T005 helper. **Read the coverage delta before
      the score**: this arm removes edges and can disconnect readouts.
- [ ] T022 [US3] Arm `break_recycling`: for each component classified
      `recycling_artifact` (T004), remove the intra-component edges that
      close the cycle, preferring the orientation evidence in Reactome's
      `precedingEvent` where available and falling back to the ratio
      heuristic only where it is not. Record which mechanism decided each
      removal.
- [ ] T023 [US3] Benchmark `break_recycling`, report through the T005 helper,
      and state explicitly whether it destroys curated causality. If it does
      and it improves the score, it is a diagnostic result and **must not be
      adopted as a default** (FR-009).
- [ ] T024 [US3] Arm `drop_diagram_bridge + break_recycling` combined;
      benchmark and report the same way. Check whether the two compose or
      one masks the other.
- [ ] T025 [US3] Record every arm's numbers in `research.md` including the
      neutral and negative ones, per constitution principle III and FR-008.
      An arm that removes cycles without improving macro-F1 is a result and
      is the single most useful thing to write down.

**Checkpoint**: SC-004 and SC-006 evaluable.

---

## Phase 6: User Story 4 — Bound the prize (Priority: P3)

**Goal**: know what perfect acyclicity could buy, so the cost of principled
loop-breaking can be judged against it.

- [ ] T026 [US4] Arm `dagify` in `bench/analysis/loop_interventions.py`: an
      approximate minimum-feedback-arc-set removal making every network
      acyclic. Label it unshippable in the manifest and report the deleted
      edge count with the score.
- [ ] T027 [US4] Benchmark `dagify` and report the macro-F1 delta against the
      current default as **an upper bound, never as a result**. If the bound
      is small, the loop hypothesis is dead and Phase 5's principled arms are
      not worth further investment; say so.

---

## Phase 7: Polish & Cross-Cutting

- [ ] T028 [P] Write the answer to the feature's question in `research.md`:
      are loops the cause, and if so how much of the gap do they own.
- [ ] T029 [P] Update `specs/004-loop-handling/quickstart.md` with the exact
      commands and observed numbers for every arm actually run.
- [ ] T030 Update `CLAUDE.md` only if a default changed — current state only,
      one line, pointing here for the rationale. If nothing shipped, change
      nothing.
- [ ] T031 File the findings that fall outside this feature as tracked items
      rather than inline scope creep: the misleading `load_dbid_to_uuids`
      naming (correct but confusing, research.md R2), and any propagator
      finding surfaced by T013.

---

## Dependencies

```
Phase 1 (T001-T002)
      ↓
Phase 2 (T003-T005)  ← blocking; nothing is reportable before this
      ↓
   ┌──┴───────────────┐
   ↓                  ↓
Phase 3 US1        Phase 4 US2      (independent of each other)
(T006-T013)        (T014-T018)
   ↓
Phase 5 US3 (T019-T025)   ← interpretation depends on US1's answer
   ↓
Phase 6 US4 (T026-T027)
   ↓
Phase 7 (T028-T031)
```

- **T010 gates T011.** An adapter that introduced cycles would make the
  control measure the adapter.
- **T004 gates T022.** The recycling classification must exist before an arm
  can act on it.
- **T021, T023, T024, T027 each require T005.** No arm is quotable without
  its coverage delta and both-arms-converged count.
- Phase 4 has no dependency on Phase 3 and should not be sequenced behind it;
  the honest-reporting obligation stands whatever the experiments show.

## Parallel opportunities

- T014 and T015 are the same file and are **not** parallel; sequence them.
- Phase 3 and Phase 4 are genuinely parallel bodies of work.
- T028 and T029 are different files.
- **Benchmark arms are NOT parallel.** One at a time, `pkill` the Julia
  server between arms. Concurrent servers have caused an out-of-memory crash
  in this workspace before.

## Implementation strategy

**MVP is Phase 3 alone.** The control answers the question that redirects
everything else, and it is deliverable on its own: a third number next to
365 and 407, with the attribution stated. If the control shows the propagator
is at fault, Phase 5's arms are expected to be neutral and the honest thing
is to run them, record the neutral results, and hand the propagator finding
to feature 003.

**Phase 4 ships regardless.** It is a correctness obligation on numbers
already reported, and it does not depend on any experiment succeeding.

## Notes on discipline for this feature

- **The attractive wrong conclusion is "removing edges improved the score".**
  It will be true for some arm. Constitution principle I is the guard:
  processing a faithful representation is DeltaSignal's job, and deleting
  curated causality to make the graph easy is not a fix. T023 exists to force
  the question to be asked out loud when the number appears.
- **Coverage before score, every time.** Breaking a cycle can disconnect a
  readout, converting a wrong answer into an unscoreable one and inflating
  accuracy. That is a loss of coverage, not a gain.
- **179 of 564 cases do not converge** and their values depend on sweep
  order. A few changed cases concentrated in TP53 with no both-converged
  support is the known artifact signature, seen twice before — not a result.
