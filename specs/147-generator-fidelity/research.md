# specs/147 research: the generator vs Neo4j, measured

The question (Adam, 2026-11-07): "I feel like something has been going on with
how the logic networks are being generated that is not correct." The audit's
answer, in three measured parts: **at every layer this audit can reach, the
generation is sound** — the error mass lives in pathway *scope* and in the
curated knowledge itself, not in the generator's faithfulness.

## P2 — material fidelity: ZERO findings

Every ReactionLikeEvent Neo4j reaches by `hasEvent*` in all 92 benchmark
pathways — 4,975 reactions, 19,535 curated participants — exists in the bundle
and touches its reaction in some carried form, with every regulator sign
correct (`p2/structural_diff.py`, `p2/summary.json`). The naive diff without
composition resolution flags ~10% of participants; all traced ones are
EntitySets whose own id never survives the recorded set-expansion policy
(14-3-3 dimer, the BH3-only set, the p-AKT CandidateSet) — policy, not defect.
The generator drops nothing material.

## P1 — lost routes: 392 → 130, and the remainder decomposes

On the frozen record (`05dc867`), re-running the specs/080 chain
(`failure_structure` → `case_oracle` → `route_breaks`):

- **130** failed cases where Reactome's own within-pathway bipartite graph has
  a route and our network does not — down from 392 at the specs/080
  measurement: the generator work since (recycle pools, release pools,
  specialization bridges, balance choices, summation supplement) recovered
  two-thirds of the class.
- First unreached step, typed:
  - **94: `component → complex`, both nodes exist, edge missing** — the
    composition hop into a complex no reaction forms (specs/031's recorded
    curation-review class).
  - **28: the route starts from a form of the gene the protocol does not
    pin** (26 of them PTK6) — the pinning protocol, known since specs/024.
  - 6: the complex node absent; 2: member → set.
- **The guard test on the 94+6**: adding the missing composition edge only
  where it closes no cycle (the specs/018-safe subset) would serve just **6
  cases** (3 pairs: TP53 ×1, RUNX1 ×2). **58 cases sit on pairs where the
  edge welds a cycle** — the exact mechanism specs/016 measured at −185 and
  refused — and 36 need nodes the network does not have. The biggest break
  class is therefore *correctly* refused, not negligently missing:
  Reactome curates these complexes as appearing with their components inside,
  and any edge faithful to that containment also fabricates feedback Reactome
  does not contain.

## P3 — connection

Within-pathway: the oracle's scope is the same `hasEvent*` closure the
generator builds from, and P2's zero plus P1's 130-only say sub-pathway
seams are not losing routes. Between pathways: the benchmark networks are
deliberately unconnected (each pathway solved alone); the connection
machinery that exists (context builds, specs/056/132/133) was measured —
identity joins refuted, curated-connection bridges passing balanced accuracy
(+0.0141, T2+F) but declined because fall-back answers silently mixed with
in-pathway ones and dragged korkut raw down through PIP3's rails.

## What this says to do (and not do)

1. **Do not hunt the generator further for the benchmark pathways.** Two
   independent methods (participant diff, route oracle) put its remaining
   within-pathway defect mass at ~6 safely-fixable cases out of 23,511.
2. **The error mass is scope**: two-thirds of validation readouts are
   untouched by the perturbation *within the curated pathway* (specs/146),
   and against experiments those read wrong 70% of the time. The measured,
   gate-passing instrument for that is specs/133's curated-connection
   fall-back — now revivable with the specs/146 machinery: serve
   connected-pathway answers as a **separate, labelled, separately-calibrated
   confidence class** instead of silently mixing them (the exact failure that
   blocked adoption). That is the registered next candidate.
3. The 26 PTK6 pinned-form cases are a pinning-protocol question (one
   pathway), and the 6 safe composition pairs are a small guarded generator
   arm if ever worth the run; the 58 welding pairs and the no-curated-route
   mass stay on the curation-review channel (issue #237's class).

## Provenance

Frozen record `20261026-1349_8b6527e / 05dc867`; oracle chain run 2026-11-07
(`p1/curator_failure_structure.tsv`, `p1/curator_case_oracle.tsv`,
`p1/route_break_steps.tsv`); structural diff `p2/structural_diff.py` with
`fetch_reaction_rows`/`fetch_composition` imported from the specs/016 oracle.

## P4 (addendum, 2026-11-07, Adam's question): the diagram-only connection arrows

Adam: "the diagram files provide edges that are not in Neo4j but should be in
our networks." Checked, and he is right — with a precise shape:

- The generator **does** consume diagrams for reaction connectivity
  (`diagram_connectivity.py`: shared entity glyphs stand in for the
  under-annotated `precedingEvent`; LNG #39) — within a pathway, drawn
  connectivity is already used.
- What nothing consumes are the diagrams' **`Interaction` links**: across the
  90 benchmark diagrams, **245 of them — 115 `pathway→entity` arrows** (an
  upstream pathway drawn as a box feeding a NAMED entity, e.g. "Signaling by
  EGFR → RAS GEFs"), **117 `entity→pathway`**, 13 pathway→pathway. They carry
  **no reactomeId and no dbId** — purely curator-drawn, in neither the event
  graph nor our networks nor the oracle this audit used. The census is
  `p4/diagram_interactions.tsv`.
- The set-relation links (`EntitySetAndMemberLink` 270, set–set 53) duplicate
  Neo4j `hasMember` — no gap there.

**Why this matters**: these arrows are exactly the curated, entity-anchored
cross-pathway joins that specs/132/133 had to approximate from event-graph
connections — and the entity anchor is the part 132's tracing said decides
right from wrong bridges (the p-AKT handoff was right; root-set changes were
wrong). The registered next candidate therefore sharpens: **context bridges
drawn FROM the diagram Interaction arrows** (joined at their named entities),
served as the separate, labelled, separately-calibrated confidence class the
specs/146 machinery provides. The within-pathway generator remains sound
(P1/P2); this is the connection half of Adam's question, answered: the edges
exist, they are diagram-only, and we now have their census.
