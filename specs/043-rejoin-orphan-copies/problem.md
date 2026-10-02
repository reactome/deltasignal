# specs/043: rejoin copies of one entity that the generator split without cause

## Verified defect (specs/041 step 3b; specs/042 research)

The generator gives each occurrence of an entity its own node (uuid). Phase 2
then joins a producer's output to a consumer's input **only along curated
`precedingEvent` pairs** (`_get_or_create_entity_uuid`). That was deliberate.
- Merging every occurrence by stId welded loops Reactome does not contain:
  specs/018, and MP-BioPath's hand-cut loops (specs/034 §11).
- Per-position copies resolve about 86% of those.

**But where a producing reaction has no `precedingEvent` link to a consumer of
the same entity, the copies are never joined**, even when Reactome curates
them as one entity flowing from producer to consumer.

| case (verified in Neo4j) | Reactome | our build |
|---|---|---|
| WNT, CTNNB1 [cytosol] R-HSA-448839 | one entity: bound by the destruction complex (195304), released by it (201685), imported to the nucleus (201669) | 3 copies: a root that only binds, a release product that only goes to the nucleus, a dissociation dead end. The curated cycle is severed. |
| WNT, AMER1 [cytosol] | the same pattern | 3 copies |
| TP53, TP53 Tetramer R-HSA-3209194 | one entity: 3 producers, 21 consumers | 4 copies. The p14ARF (6804996) and USP7 (3215310) products have **no following events**, so they dead-end in dissociations. |

These copies account for 18 traced experimental errors, and an estimated
499 root-split stIds in 50 pathways (census, specs/042 derivation-fable §0).

## The fix: connect copies only where a Reactome diagram connects them (Adam, 2026-10-02)

**Criterion.** Copies P (an output of reaction r1) and C (an input of
reaction r2) of the same entity are joined **only if some Reactome diagram
draws r1's output and r2's input on the same glyph.** The diagram is the
curators' own statement of which occurrence connects to which.

- **Today the generator already merges on shared glyphs, but reads only one
  diagram** per pathway: the pathway's own, else its covering one
  (`diagram_connectivity.covering_diagram_stid`).
- **The fix:** consider every diagram (Reactome 97 layout + graph JSON in
  `~/reactome-diagrams/97`) that draws both r1 and r2 within the pathway's
  events.
- Behind `LNG_DIAGRAM_ALL=1` (default off) until measured.

**The verified cases (diagram check, 2026-10-02):**
- **TP53 Tetramer R-HSA-3209194: connected by the diagram.**
  - In R-HSA-6806003 (and 9723905), glyph 6 is the output of 6804996
    (p14ARF sequesters MDM2), 3215310 (USP7 deubiquitinates TP53) and
    "TP53 forms homotetramers", and the input of "MDM2 binds TP53" and
    "TP53 binds the CCNG1 gene".
  - Our 4 copies are a **generator artefact**: the connecting diagram is not
    the one read. This covers the 10 CDKN2A cases.
- **WNT CTNNB1 [cytosol] R-HSA-448839: NOT connected by any diagram.**
  - Association with the destruction complex (195304) draws it as glyph 2838,
    in R-HSA-195253 and its copies.
  - Release (201685) and nuclear import (201669) draw it as glyph 4423, in
    R-HSA-201681 and its copies.
  - No diagram shares a glyph between the release output and the association
    input. **Our split is faithful to the drawings.**
  - The APC KO sign inversion therefore comes from Reactome's representation
    (free cytosolic β-catenin is curated only as "released from the
    destruction complex"). It is a solver/semantics question (specs/042), or a
    curation question for Adam. It is **not** a rejoin.
  - The specs/042 research entry calling WNT a generator artefact is
    **corrected by this check**.
- **Open curation question (Adam):** when the same compartment entity appears
  in *different sub-pathway diagrams*, is it the same pool? Glyph identity
  answers the question within a diagram, not across diagrams. Until Adam
  decides, cross-diagram occurrences stay separate.

**Guard (kept):** a join that would create a strongly connected component
absent from Reactome's own reaction graph is refused and counted.

## What it does NOT do

It adds no solver rule. A rejoined binding/release cycle with no supply (WNT)
may still read wrong until a conservation or turnover rule handles it (specs/039
pools; specs/042 R1). This spec measures the rejoin alone first.

## Census before any arm

Count:
- joins per pathway;
- new strongly connected components, and their sizes;
- rejected joins (by condition 1, 2 or 3);
- cases whose route changes.

List every join in WNT and TP53 by name, for review.

## Gates

The same as specs/040, both axes:
- curator held-out > +15, p < 0.05;
- experimental ≥ 0, not significantly negative;
- no pathway loses > 10.

Plus:
- **zero new strongly connected components absent from Reactome;**
- relabel churn no worse than the control.

The rejoin is a faithfulness fix: it makes the network match Reactome. So
**neutral results on both axes are acceptable for adoption, with Adam's
agreement**, as variant sharing was in specs/020.
