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
