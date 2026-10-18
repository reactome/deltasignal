# specs/070 research: the curator-fidelity audit

Canonical build `20261017-1246_5e663d9`, Release 97 Neo4j and diagram layouts.
All files in `bench/analysis/070/`.

## 1. Reaction wiring: faithful

`participant_audit.py` → `participants_r97.tsv` (header only). All 4,975
reactions under the 92 catalog pathways are in the networks, and for every
role (input, output, catalyst, positive and negative regulator) the entities
joined to the reaction's copies are exactly the curated ones, a set standing
for its members. Negative control: a fake input added to R-HSA-5651992 is
reported as missing twice (DNA Damage Bypass and DNA Repair).

## 2. Hand-offs the network makes: faithful or derived by design

`handoff_audit.py` → `handoffs_r97.tsv`: 7,912 (A, B, role, route) hand-offs.

| route | n |
|---|---|
| direct (A's output node is B's participant node) | 6,946 |
| variant_pool / variant_split | 503 |
| assembly (component into a complex no reaction forms) | 279 |
| specialization (specs/067) | 131 |
| assembly + set_member, other mixed | 53 |

Direct hand-offs, order × diagram:

| | shared | link | separate | apart |
|---|---|---|---|---|
| forward | 4,666 | 440 | 69 | 407 |
| reverse | 365 | 3 | 2 | 6 |
| unordered | 825 | 72 | 16 | 75 |

Neo4j participants: 6,793 `same`, 661 `set`, 184 `component`, 274 `none`
(all `none` are derived routes: specialization, or assembly through a nested
complex the closure did not reach). The network does not join reactions the
curators did not join: every direct hand-off passes a curated entity from
its producer to its user, 91% of them on a drawn glyph or dashed link. The
`apart` ones cross sub-pathway diagrams, almost all in precedingEvent order.
The `separate` ones (87) are the only direct hand-offs where a common diagram
draws two glyphs; not yet traced individually.

Generator-derived links outside hand-offs, which no curator draws:
- `depletion`, 727 edges in the catalog (Hedgehog 177929 81, DNA Repair 63,
  RAF 58, ERBB2 55, ...): catalyst → input for dephosphorylations, and for
  ubiquitinations the ligase → every network copy of the free substrate.
  An inference from mass action, not a curated participant.
- `assembly` 174,212 and `dissociation` 65,911 boundary edges, from
  hasComponent (specs/030, specs/015).

## 3. Hand-offs the curators record that the network does not make

`missing_handoffs.py` → `curated_handoffs_r97.tsv`: of 7,003 exact curated
hand-offs (A outputs E, B uses E, same pathway, cofactors excluded), 5,623
are present and 1,380 missing.

| diagram | missing | of which reverse / unordered / forward |
|---|---|---|
| shared glyph | 459 | 229 / 200 / 30 |
| separate glyphs | 430 | 185 / 240 / 5 |
| never on a common diagram | 490 | 98 / 387 / 5 |
| dashed link | 1 | 1 / 0 / 0 |

Of the 459 shared ones, 450 are drawn only on a sub-pathway's own diagram
(`subpathways_r97.tsv` lists the pathways under each catalog pathway), 8 on
the pathway's own diagram, 1 elsewhere. The generator's
`diagram_shared_product_pairs` reads only the covering diagram; a sub-pathway
with its own diagram is a box there. Only precedingEvent pairs then join
those reactions, and recycling returns are rarely precedingEvent pairs.

Mechanism (`missing_mechanism.py` → `missing_mechanism_r97.tsv`):

| diagram | copy | seam | A output absent |
|---|---|---|---|
| shared | 349 | 74 | 36 |
| separate | 393 | 4 | 33 |
| apart | 455 | 13 | 22 |

- `copy`: A's output and B's input are two plain copies of E. In 814 of the
  1,198 such cases A's output copy feeds nothing at all.
- `seam`: B reads E through a specs/046 variant pool (one node per member,
  fed from root copies) while A outputs E as one node; A's copy is a dead end.
  Traced: RAF R-HSA-5672980 → 5672972 (5672717), Class I MHC
  R-HSA-983140 → 983156 (983131, 206 variant copies), 1912422
  R-HSA-9604834 → 9604829 (9604877). These are 3 of the 8 own-diagram misses.
  The other 5 are regulator glyphs, which the generator's reader skips
  (inputs and catalysts only).
- `A output absent`: A's output node for E carries a different entity id
  (a member); not traced.

Traced example of `copy`, PIP3/AKT (R-HSA-1257604): PTEN (R-HSA-199420) has
five network copies. R-HSA-6807126 outputs one that feeds nothing; R-HSA-6807106
reads the root copy. PTEN Regulation's own diagram (R-HSA-6807070) draws one
PTEN glyph as 6807126's output and 6807106's input.
