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
- **DSB Repair, RUNX1:** converged; the broken cases are PARP1, PARP2, FEN1,
  XRCC1, LIG3, POLQ, RAD52 (DSB) and PTPN11, SRC (RUNX1). This bullet first
  said "knockouts drain the loop to 0"; that was wrong. The overexpressions
  break too, reading UP, and the mechanisms differ. §5 has the traces.
- **PIP3/AKT: the predicted fix.** PTEN released by R-HSA-6807126 now returns
  to R-HSA-6807106; 42 curator and 13 experimental cases fixed, none broken.

So the drawn networks are more faithful and score worse, because the solver
cannot handle the curated recycling loops they contain: it fails to converge
on them or collapses them to zero. This is the specs/069 finding at catalog
scale. The loss is the solver's loop handling, not a wrong hand-off.

## 5. Arm dsub: the losses traced

One case per pathway, solved under the code defaults on the canonical bundle
and on the `dsub` bundle (`bench/analysis/070/dump_case.jl`, pins from each
run's case log), every node compared. The escape point is the first node that
moves in `dsub` and not in the canonical solve.

**DSB Repair (−195): an open carrier loop's release read as supply.**
PARP1 knockout. Every alt-EJ (MMEJ) gene now reaches 12 HR and SSA readouts
in both directions (knockout DOWN, overexpression UP; curators: no change).
The signal leaves MMEJ at RBBP8 homotetramer (R-HSA-5684097). MMEJ's last
step, "LIG3 ligates remaining SSBs in MMEJ" (R-HSA-5687675), releases RBBP8,
and the MMEJ diagram draws that glyph as the input of MMEJ's first step
(R-HSA-5687465). The join makes the released copy the free RBBP8 node, which
also feeds HR end processing ("MRN complex binds CDK2 and RBBP8",
R-HSA-5684081; the two reactions are never on a common diagram). The free node
reads OR(SIRT6 deacetylation, MMEJ release), so stalling MMEJ halves it and
HR stops. The generator lists the loop in `carrier_loops.csv` but as OPEN
(SIRT6 also produces free RBBP8), and the default `closed_twin` acts only on
closed loops. The faithful build has 83 carrier loops in DSB Repair (24
closed, 59 open); the canonical has none. Over the catalog the join adds 173
closed and 191 open loops (DNA Repair 114, Rho GTPases 39, WNT 53, DSB 83).
Under conservation a carrier released unchanged is not a supply: if MMEJ never
starts, RBBP8 is never bound and its free level is unchanged.

**RUNX1 (−93): a phosphorylation cycle closed into a pool, and the enzyme
knocked out.** PTPN11 knockout. Free RUNX1 (R-HSA-8865295) is made by
translation and by "PTPN11 dephosphorylates RUNX1" (R-HSA-8937767); in the
canonical build the dephosphorylated RUNX1 was a dead-end copy. Joined, SRC
phosphorylation and PTPN11 dephosphorylation form a two-path cycle-balance
pool (specs/039) whose return path's carrier is PTPN11 itself. With PTPN11 at
0 the return rate is 0, all RUNX1 ends phosphorylated (π absorbed), and free
RUNX1 and every target read 0; SRC overexpression does the same through the
forward rate. That is the mass-action answer for a missing phosphatase. The
curators answer no change. Before any rule, the ground truth needs checking.
The 2019 answers may predate the SRC and PTPN11 reactions, or may not treat
p-7Y-RUNX1 as inactive.

**TP53 (−365): MDM2 and the tetramer pooled, the gene's supply excluded.**
The join returns the TP53 tetramer from USP7 deubiquitination (R-HSA-3215310)
and from the p14ARF release (R-HSA-6804996), and p-MDM2 from its own steps.
That makes five cycle-balance pools (canonical: none), MDM2 (R-HSA-69486)
among them. MDM2 knockout (100 cases broken, MDM2 both directions): the
pinned root (R-HSA-3700996) zeroes translation (R-HSA-3700992), but every pool
form stays at baseline. Traced in the pool update: the supply list of every
MDM2 form is empty. specs/039 finding 6 excludes from a pool's supply any
producer the pool's own forms reach. TP53 transcribes MDM2, so MDM2
translation is downstream of the MDM2 pool, and the gene's only supply is
dropped. With no supply the pool reads s = 1, so no perturbation of MDM2
can reach TP53. The other broken TP53 cases are non-converged (663 of
1,836, from 51; 361 of the 511 broken). The specs/066 LM fallback skips
components that hold pool nodes, and TP53's main component now holds five
pools. The traced MDM2 knockout itself converged.

Three mechanisms, none of them a wrong hand-off:
1. an open carrier loop's release read as supply (DSB);
2. a cycle pool absorbed when its reverse enzyme is knocked out (RUNX1), a
   ground-truth question;
3. a pool whose supply is in feedback with it reads no supply (TP53/MDM2),
   and non-convergence the fallback does not reach.
