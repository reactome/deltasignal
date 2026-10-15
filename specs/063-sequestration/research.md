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

## R3: Reactome's annotations do not separate them either

`bench/analysis/063/annotations.cypher` (scratch query recorded here): for
the 24 target complexes and a random 300 of the other dead-end complexes
(seed 63), whether the complex is used anywhere in Reactome (input,
catalyst, regulator) and the GO process of the reaction producing it.
- Targets: 0 of 24 used anywhere; 3 have "negative regulation of BMP
  signaling pathway" (the I-SMAD complexes).
- Sample: 12 of 300 used anywhere (4%); 5 have any GO process, 3 a negative
  regulation one.
So the dead end is curated (not our generator), and "sequestering" is not
annotated. specs/063 therefore measures a narrow and a broad rule.

## R4 (2026-10-15): both arms refuted

Code `47614c3`, canonical build, against the control `64d23cb`, RAF set aside
where it applies (it changes nothing on curator held-out).

| arm | curator held-out | tuning | experimental | convergence (curator / exp) |
|---|---|---|---|---|
| control | — | — | — | 1,652 / 215 |
| `seqsole` (primary) | **−1,180** (121 / 1,301), 25 of 71 pathways move | −154 | −20 (1 / 21) | 1,488 / 187 |
| `seqall` (record) | **−1,686** (296 / 1,982), 46 of 71 | −239 | −36 (3 / 39) | 1,205 / 164 |

Both fail every part of the gate. The prediction for `seqall` (fails) held;
for `seqsole` I gave even odds, and it is a large loss.

**Where the targets went** (`seqsole`): USP18 +56 (every case, as predicted),
SOCS6 +8, SH2B3 +8, IL13RA2 +4, CER1 +4, SOCS1 −12. So the mechanism works
where the partner is a quiet root.

**Where the losses come from.** RUNX1 −607 and DSB Repair −414 carry 86% of
the held-out loss; 979 of the 1,301 broken cases are expected NOT and now
move. DSB loses 414 to two edges: UBXN1:BRCA1:BARD1 -| BRCA1:BARD1 (the
partner has 15 other uses inside DSB's giant loop) and p-PHF6:SMARCA5 -|
SMARCA5. The complex is built FROM the partner, so the edge is an own-product
depletion: any perturbation that lowers the partner lowers the complex, which
de-represses the partner by its own fall (the specs/060–061 motif). Most
perturbations in these pathways move the partner, not the sequesterer, so
the rule mostly injects feedback, and convergence drops (curator 1,652 →
1,488).

**Not adopted.** A version without the self-feedback is the obvious next
step: under multiplicative AND the complex's fold is the sequesterer's fold
times the partner's, so the part of the complex's change that is not the
partner's own is the sequesterer's fold. Depleting the partner from the
sequesterer directly (A -| B) carries exactly that. Registered as amendment 1
in spec.md before any code.
