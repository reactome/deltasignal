# specs/070: are the networks what the curators drew?

**Status:** audit done 2026-10-18; arm `dsub` pre-registered below, not yet run.

Adam, 2026-10-18: "now that you understand how to check the networks vs what
was intended by the curator perhaps you could check to make sure the networks
are correct relative to what was intended ... perhaps a lot of our problems is
from creating links that are not faithful to what curators intended."

## The audit (method)

Scripts in `bench/analysis/070/`, inputs exported from Release 97 Neo4j
(`reactions.cypher`, `entities.cypher`, `subpathways.cypher`) and the
Release 97 diagram layouts (`~/reactome-diagrams/97`, 1,382 diagrams), run on
canonical build `20261017-1246_5e663d9`.

1. **Each reaction's wiring** (`participant_audit.py`): for every reaction of
   every catalog pathway (4,975), the entities the network joins to any copy
   of it as input, output, catalyst, positive or negative regulator, compared
   with Neo4j role by role (a set's member counts as the set). Checked with a
   planted fake participant, which it reports.
2. **Hand-offs the network makes** (`handoff_audit.py`): every place a
   reaction A's output node is read by a reaction B (directly, or through up
   to three generator-derived edges: assembly, variant_pool, variant_split,
   specialization, set_member), at stable-id level, cofactors excluded. Each is
   checked against Neo4j participants, precedingEvent order, and the diagrams:
   `shared` (some diagram draws A's output and B's participant as one glyph),
   `link` (a dashed link joins the two glyphs), `separate` (drawn on a common
   diagram as two glyphs), `apart` (never on a common diagram).
3. **Hand-offs the curators record that the network does not make**
   (`missing_handoffs.py`): A outputs E and B uses E, both in the pathway,
   exactly (no set inference), cofactors excluded.

## Findings

- **Every reaction is wired as curated.** No reaction is missing, and no
  input, output, catalyst or regulator is extra or missing (up to set
  membership).
- **Links the network makes are faithful.** Of 7,912 hand-offs, 6,946 are
  direct. Of those, 6,313 are drawn as one glyph or joined by a dashed link
  and supported by Neo4j. The rest pass a curated entity between reactions
  that sit on different diagrams (cross-sub-pathway), almost all in
  precedingEvent order. The 365 hand-offs that run against precedingEvent
  (reverse) are drawn as shared glyphs: curated loops, as in specs/069.
  Hand-offs through derived edges (966) are the generator's own inference by
  design: assembly into a complex no reaction forms (specs/030),
  specialization bridges (specs/067), set variants (specs/046).
- **Links the curators drew are missing.** Of 7,003 curated exact hand-offs,
  1,380 are not in the network. 430 of those are drawn as separate glyphs and
  490 never on a common diagram, so leaving them unjoined agrees with the
  drawing. **459 are drawn as one shared glyph.** 450 of these are drawn on a
  sub-pathway's own diagram, which the generator never reads: it reads only
  the covering (top-level) diagram, where such a sub-pathway is a box. Top
  pathways: DNA Repair 106, 446728 75, Signaling by Rho GTPases 50, DSB Repair
  38, WNT 30, TP53 28, 170834 23, PIP3/AKT 21. They are mostly reverse (229) or
  unordered (200): recycling returns, the loop shape specs/069 found intended.
  Traced: in PIP3/AKT, PTEN released by R-HSA-6807126 is a dead-end copy, while
  PTEN Regulation's diagram (R-HSA-6807070) draws it as the glyph R-HSA-6807106
  consumes.
- 8 more are drawn on the pathway's own diagram but missing. Of these, 5 are
  regulator glyphs, which the generator's diagram reader does not read
  (inputs and catalysts only); 3 not yet traced.

## Arm dsub: read sub-pathway diagrams (pre-registered 2026-10-18)

Generator flag `LNG_DIAGRAM_SUBPATHWAYS` (LNG branch `diagram-subpathways`,
`ddbd828`, default off): `diagram_shared_product_pairs` also reads the own
diagram of every pathway under the catalog pathway, still restricted to the
pathway's own reactions, and merges each drawn pair exactly as it merges the
covering diagram's. Variant build `<date>_ddbd828_dsub`, run with
`scripts/run_arm.sh dsub --catalog <id>`, both axes, against the canonical
scoring (`results/6fb4329`).

This is a faithfulness fix: the curators drew these hand-offs. Rule decisions
set RAF aside (specs/054 D1).

**Gate (registered):** curator held-out net ≥ 0, experimental net ≥ −5, no
pathway below −10, RAF set aside. A faithfulness fix is adopted unless it
costs accuracy; it does not need a significant gain. Otherwise not adopted,
and the losses are traced. The usual gain gate (held-out net > 0 at
McNemar p < 0.05) is reported alongside it.

Because the change adds hand-offs, root pins can change. Compare on all cases
valid in both runs as well as holdout_report's paired cases (specs/069).

**Predictions:**
1. The audit re-run on the variant reports at most 20 of the 459 shared-glyph
   hand-offs still missing.
2. Most added hand-offs close recycling loops (reverse or unordered). Where
   the 2019 answers assume one pass (as in DNA Damage Bypass, specs/069),
   cases will be lost. No sign is predicted for the gate.
3. The exposed pathways include the tuning pathways WNT, TP53 and PIP3/AKT,
   so experimental can move; no sign predicted.
