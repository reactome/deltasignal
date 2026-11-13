# specs/151: where the Release-98 drift lives — tracing the re-curated pathways

**Type: analysis** (like specs/079/092/094): no rule, no gate, no adoption. The
questions are registered before looking, the answers are traces of named cases.

## Why

specs/150 R1: the frozen rulebook replicates on Release 98 (held-out −0.32pp,
inside the ±2.0pp band), but the drift is concentrated, not diffuse: on
identical case keys, WNT −25, IL-4/13 −24, RHO-PKNs −16, and ERBB2 +18. A
replication that passes on the aggregate can still be hiding a mechanism in
the concentration. The trace rule (Adam): cases, not aggregates.

## Questions (registered before measurement)

- Q1. For each of the four pathways, which exact cases changed class
  (fixed / broke) between `20261026-1349_8b6527e` (97, solver 05dc867) and
  `20261113-0316_8b6527e_r98` (98, solver e967a1b), on identical
  (pathway, gene, direction, key_output) keys?
- Q2. Are the broken cases concentrated in a few perturbations or readouts
  (the specs/148-style concentration test), and is the prediction change a
  class flip at a threshold, a route change (untouched <-> routed), or a
  rail/switch flip (specs/108's artefact classes)?
- Q3. For at least one representative broken case per pathway: what did the
  curators change? Named in Reactome terms (new/removed/rewired reactions or
  participants between releases), verified in the 98 Neo4j against the 97
  network's route for that case.
- Q4. Is any loss our machinery's fault (a rule misreading the new curation —
  actionable) rather than the curation legitimately moving the label's route
  (not actionable)? The specs/108 lens: was the 97 answer right through an
  artefact?

## Non-goals

No rule changes, no re-scoring, no gate. If Q4 finds an actionable defect it
is filed as its own follow-up, not fixed here. RAF is not in scope (set
aside; its drift is zero anyway).

## Data

The two case dumps named in Q1 (both on disk under
`~/deltasignal-catalogs/builds/<id>/results/<sha>/`), the two catalog bundles
for the four pathways, and the Release-98 Neo4j (`reactome-neo4j`).

## Amendment 1 (2026-11-13, during the trace): the drift is a build defect, not curation

Q3/Q4 answered before Q2 finished, because the first trace found a root cause
that supersedes the questions' framing: **Release 98 removed every HGNC
`referenceGene` link** (17,714 ReferenceGeneProducts with one in 97; **0 of
106,316 in 98**; nothing at all links to the HGNC ReferenceDatabase node in
98). The frozen generator's `get_reference_entity_id` keys on
`rd.displayName = 'HGNC'`, so on Release 98 it returned None for every
protein, silently disabling everything protein-identity-based in the build:
specs/101 balancing (`balance_dropped: 0` across the whole 98 build log),
output-slot binding (`output_slots_bound` absent, all slots free), and the
set-member signatures behind same-entity joins. The traced PKN case: 97
enumerates PDPK1-phosphorylates-PKN as 12 balanced copies (PKN_i → p-PKN_i);
98 as the full 36-copy cross product, so a PKN2 over-expression raises
p-T774-PKN1 — mechanism the curation does not contain. WNT's +49% node growth
on an identical reaction set is the same defect. The 98 BUILD.json also
records an empty environment (no PYTHONHASHSEED=0, no LNG_EMIT_ONE_SIDED=1) —
both default correctly in the frozen code, but the seed matters for
reproduction and should have been set.

Consequence: specs/150 R1 measured the frozen rules on a mis-built catalog.
The band still passed, but "the replication" was not run on the construction
the freeze intends. The serving switch is being reverted; the fix (a
release-robust reference lookup, gated on reproducing the frozen 97 catalog
identically) and the true 98 replication are follow-up work.
