# specs/063: sequestration — groundwork before any rule (2026-10-15)

No rule is registered yet. This records what a rule would have to separate.

## R0: what the target cases are

The curation review (`~/deltasignal-curation-review_2026-10-14.tsv`, from
`bench/analysis/curation_gaps/`) labels 348 held-out no-path cases
"sequestration?": the perturbed gene's route ends at a complex nothing in the
pathway uses. Expected directions on the canonical results (`64d23cb`):

| pattern | cases |
|---|---|
| inhibitor-like (knockout UP, overexpression DOWN) | 274 |
| activator-like (knockout DOWN, overexpression UP) | 72 |
| other | 2 |

The 274 fit "A binds B into a dead end, so A lowers free B". Most of the 72
fit the double negative: A sequesters an inhibitor (TP53 → BH3-only:BCL2,
14 cases). PTPRK, FLRT1 and part of PTCH do not obviously fit either.

## R1: dead-end binding complexes are everywhere

`bench/analysis/063/deadend_census.py` on `20261012-1401_01206ca`
(`census_01206ca.tsv`): a binding (≥2 non-cofactor inputs, one non-cofactor
output) whose complex feeds only dissociation sinks.
- 6,295 such complexes in 77 pathways, 15,429 partner edges. Most are
  legitimate terminal products (a pathway's final signalling complex). A
  rule on all of them is close to the general consumption-by-binding rule
  refuted on 2026-07-16.
- Narrowed to bindings with one input used nowhere else (the would-be
  sequesterer) and exactly one input used elsewhere: 1,477, in 56
  pathways. It covers 8 of 24 target complexes (USP18, IL13RA2, DLK1,
  SOCS6, SH2B3, SOCS1, CER1:BMP, BH3:BCL2) and misses follistatin, FSTL3,
  LEFTY, CER1:NODAL, CDK6:RUNX1, PMEPA1, SMAD7 (their sequesterer has other
  uses in the pathway).

## R2: Reactome's event order does not separate them

`bench/analysis/063/followers.cypher`: of the 532 reactions producing a
dead-end complex, 49 have a reaction listing them as `precedingEvent`, and
483 do not. The route-gap cases of the curation list (KCTD1, SOCS1 in IL-4,
MET endocytosis, CBL in KIT) are mostly in the 483, next to the
sequestration cases. Only the PDGF ligand–receptor complex and the
NODAL ligand–receptor complex have a curated follower.

## What this means for a rule

No structural signature found so far separates a sequestering complex from
a terminal product. Candidates still to test, before registering:
- curated semantics on the complex or reaction (GO terms, a
  `NegativeRegulation` whose regulator is the complex elsewhere in Reactome,
  literature reference type);
- Reactome-wide use: a complex that is a dead end in this pathway but used
  in another pathway is a route that continues elsewhere (cross-pathway,
  specs/056), not sequestration.
