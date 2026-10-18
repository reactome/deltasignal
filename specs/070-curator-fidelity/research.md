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

## 4. Arm dsub: the drawn hand-offs, joined

Variant build `20261018-0815_ddbd828_dsub` (`LNG_DIAGRAM_SUBPATHWAYS=1`),
results `results/40adace/dsub`, against canonical `results/6fb4329`.

**Prediction 1 holds.** The audit re-run on the variant leaves 9 of the 459
shared-glyph hand-offs missing (the seams and regulator glyphs); 583
hand-offs are added. Network size falls (nodes 211,767 → 204,187, edges
500,965 → 490,333), as merges do.

**Gate fails.**

| reading | fixed | broken | net | p |
|---|---|---|---|---|
| curator held-out, paired (gate; RAF set aside, none in held-out) | 23 | 194 | −171 | < 1e-4 |
| curator tuning, paired | 169 | 401 | −232 | < 1e-4 |
| experimental, paired (RAF set aside) | 21 | 40 | −19 | 0.020 |
| curator, all valid in both | 213 | 840 | −627 | |
| experimental, all valid in both | 30 | 59 | −29 | |

Curator held-out in release 89.91% → 88.31%; experimental 606 → 577 of 849
(71.38% → 67.96%). 1,384 curator and 73 experimental cases change pin set.

Per pathway (all valid in both, `perpathway.py`):

| pathway | fixed / broken | net |
|---|---|---|
| TP53 | 146 / 511 | −365 |
| DSB Repair | 6 / 201 | −195 |
| RUNX1 | 2 / 95 | −93 |
| IL-2, WNT, TGF-beta, Hedgehog | | −8 to −6 each |
| RUNX2 | 12 / 1 | +11 |
| PIP3/AKT | 42 / 0 | +42 (experimental 13 / 0) |

The added hand-offs are recycling returns (TP53: 47 unordered, 30 reverse,
42 of 78 the TP53 tetramer R-HSA-3209194 returned to earlier reactions; DSB
24 of 40 reverse; RUNX1 12 of 18 R-HSA-8865295).

- **TP53: non-convergence.** Cases whose solve does not converge rise from
  51 to 663 of 1,836; 361 of the 511 broken are non-converged in the arm.
  The returned tetramer closes a large positive cycle and the iteration cap
  picks the class (specs/066's failure, at scale).
- **DSB Repair, RUNX1: knockouts drain the loop.** Converged; the broken cases
  are knockouts (PARP1, PARP2, FEN1, XRCC1, LIG3, POLQ, RAD52; PTPN11, SRC)
  whose readout now reads 0 where the answer is no change or up. A knockout
  on a closed positive cycle collapses it to the all-zero state (specs/014).
  Traced at this level only.
- **PIP3/AKT: the predicted fix.** PTEN released by R-HSA-6807126 now returns
  to R-HSA-6807106; 42 curator and 13 experimental cases fixed, none broken.

So the drawn networks are more faithful and score worse, because the solver
cannot handle the curated recycling loops they contain: it fails to converge
on them or collapses them to zero. This is the specs/069 finding at catalog
scale. The loss is the solver's loop handling, not a wrong hand-off.
