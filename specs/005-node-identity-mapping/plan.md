# Implementation Plan: Every node maps, both ways, including on the diagram

**Branch**: `005-node-identity-mapping` | **Date**: 2026-09-10 | **Spec**: [spec.md](./spec.md)

## Summary

The request was for a perfect two-way mapping between Reactome entities and
generated nodes, and specifically for knowing which *diagram glyph* a node
came from, since the same entity can be drawn in two places.

Phase 0 found the work is smaller and better-founded than the request
implies, in three ways. The positional provenance already exists — nodes are
minted from a registry keyed `(entity_dbId, reaction_uuid, role)` and
partially exported already. That same triple **uniquely identifies a
glyph**: 0 of 156 diagram triples resolve to more than one glyph, because an
entity drawn nine times is drawn once per reaction. And the gap is entirely
one category — every unresolved entity measured is an EntitySet, 29 of 29
and 9 of 9, so the exclusion list should come out empty rather than merely
short.

It also found a blocker and a retraction. The existing export is **31.3%
orphaned and it is the causal 31%**: 100% of catalyst and 100% of regulator
rows point at nodes absent from the network (LNG #67). And the spec's claim
that Reactome's set graph contains cycles is **false** — zero cycles, max
nesting depth 5. The spec has been amended.

## Technical Context

**Language/Version**: Python 3.11 (generator, benchmark); Julia 1.10 only if
the API later serves resolution, which is out of scope here.

**Primary Dependencies**: existing only — py2neo for set membership,
the diagram layout JSON already read by `src/diagram_connectivity.py`.

**Storage**: two new per-pathway CSVs beside the existing network files,
plus the fix to one existing CSV.

**Testing**: pytest in the generator, tiered (`not database and not
integration` by default). The negative control is the load-bearing test.

**Target Platform**: Linux; Neo4j Release97; diagrams under
`LNG_DIAGRAM_DIR`.

**Project Type**: data-generation pipeline plus its offline consumer.

**Performance Goals**: none. Set closure is depth ≤ 5 over 5,440 nested
sets — trivial. Do not add a cache; the stale-cache trap has already
invalidated an A/B on this project.

**Constraints**: stable ids primary, never database ids as keys. Every
persisted row carries its Reactome release. Regeneration into a populated
directory is a silent no-op unless caches are removed first.

**Scale/Scope**: ten pathways, 16,843 nodes, 10,378 existing context rows of
which 3,246 are orphaned, 157 and 101 participating entities in the two
pathways measured.

## Constitution Check

*GATE: checked before Phase 0, re-checked after Phase 1.*

| Principle | Status | Note |
|---|---|---|
| I. Processing is DeltaSignal's job | **PASS, and it decides the location** | The mapping is *representation* — which node stands for which curated entity, and which glyph drew it. That is the generator's knowledge, so it is emitted at generation time rather than reconstructed downstream. Reconstruction in DeltaSignal was considered and rejected: the information is lossy by the time it arrives, and a consumer querying Neo4j at read time would silently mix releases. |
| II. Measure, then claim | **PASS — and it caught me** | R1 retracts a cycle claim I put in the spec on plausibility alone. Amended in `spec.md`. US4 measures the combining rule rather than asserting mean. |
| III. Negative results are results | PASS | SC-007 requires the losing combining rules recorded. R1 is itself a recorded negative. |
| IV. New behaviour ships default-OFF | PASS | The combining rule defaults to current behaviour until measured. Emitting a new CSV changes no prediction. |
| V. Honest solver reporting | **PASS, extended** | The same principle applied to data: FR-006 forbids silent absence, FR-009 forbids reporting a partial set as whole, FR-011 forbids a declared-but-empty column. A partial result presented as complete is the same defect as a fake convergence report. |
| VI. API boundary is a trust boundary | N/A | No API change here. Serving resolution over HTTP is a later feature. |

**Post-Phase-1 re-check**: unchanged. One thing to watch — the contract puts
release-checking on the *consumer*, which is weaker than enforcing it at the
boundary. Accepted because there is no boundary yet; when the API serves
resolution, principle VI applies and it must be enforced there.

**Justified deviation**: none.

## Project Structure

```text
logic-network-generator/
├── src/
│   ├── logic_network_generator.py     # fix #67; emit node_resolution.csv
│   ├── set_resolution.py              # NEW — recursive closure to leaves
│   └── diagram_connectivity.py        # stop discarding glyph ids
├── scripts/validate_logic_network.py  # bidirectional completeness check
└── tests/
    └── test_resolution_negative_control.py   # NEW — must fail on corruption

deltasignal/
└── bench/
    ├── benchmark_mpbiopath_cases.py   # resolve set readouts; report denominator
    └── analysis/glyph_join_report.py  # NEW
```

**Structure decision**: the mapping is generated upstream and consumed
downstream, per principle I. DeltaSignal gains no ability to invent a
mapping — only to read one and to say when it is absent.

## Phase sequencing

**Prerequisite: fix LNG #67.** Not a related cleanup. Catalysts and
regulators are 31.3% of the export, they are where causality lives, and R4
shows they are also 25 of the 42 glyph-join misses. Everything downstream
measures wrong until this is fixed.

**MVP: US1 + US2 together.** US1 (set readouts) is a *consequence* of US2's
mapping, not separate work — once resolution recurses to leaves and is
recorded, the 177 cases resolve. Splitting them would build the same table
twice.

**US3 after.** Glyph identity needs US2's table to hang off, and the
benchmark does not need it. It is the pathway-browser deliverable.

**US4 last.** The combining rule only matters once sets resolve at all, and
it is the one question Phase 0 deliberately left unmeasured.

## Risks

- **The denominator moves and hides everything else.** Recovering ~177 cases
  takes the scored set from 564 to roughly 741. Any accuracy comparison
  across that boundary is meaningless. Every arm after this must state both
  denominators, and feature 004's remaining loop arms must be **re-baselined
  against the corrected case set** rather than compared to 365/564.
- **Sequencing against feature 004.** Both change the case set. 005 lands
  first because it fixes a measurement defect, whereas 004's remaining arms
  are hypothesis tests that are better run against a correct denominator.
  004's control result (403 vs 405 vs 364) stands regardless — it was
  measured on a common case set scored in all three arms.
- **Glyph ids may not survive a release bump.** Unresolved; only Release97
  diagrams are on disk. FR-005 mitigates by recording the release, but if
  they do move, anything persisted downstream must be regenerated rather
  than migrated.
- **"Every unresolved entity is a set" is measured on two pathways.** PIP3
  and Cell Cycle Checkpoints. The other eight may hold a category R5 did not
  see, which is precisely why the exclusion list exists and why a non-empty
  one is a finding rather than a failure.
- **A completeness check that cannot fail.** The known failure mode here:
  `_decomposed_ids` passed 11 of 11 while masking 18 dropped catalysts. The
  negative control is not optional polish.
