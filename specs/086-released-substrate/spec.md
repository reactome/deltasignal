# specs/086: a substrate released from the complex that degrades it (analysis, design pending)

**Status:** analysis, 2026-10-21. No rule registered yet; the design needs Adam's approval.

Adam, 2026-10-21: "yes" to counting the pattern behind the APC inversion.

## The defect (traced)

**APC in Signaling by WNT (R-HSA-195721): both directions inverted.** APC knockout reads
β-catenin and the canonical WNT targets DOWN; overexpression reads them UP. The experiment
and the curators agree on the opposite (7 experimental cases, and 12 curator readouts
under each direction).

**What Release 97 curates.**
- One producer of free cytosolic β-catenin (CTNNB1, R-HSA-448839): "Beta-catenin is
  released from the destruction complex" (R-HSA-201685).
- One consumer: "Association of beta-catenin with the destruction complex" (R-HSA-195304).
- No synthesis in the pathway.

**Our network** has two cytosolic copies of β-catenin:
- a **starting copy** (no in-edge) that binds the destruction complex;
- a **released copy**, made by R-HSA-201685, which translocates to the nucleus and drives the
  targets.

The only route between them is 8 reactions through the APC-containing complex: binding,
five CK1/GSK3 phosphorylations, PP2A dephosphorylation (the complex catalyses its own
step), and release. So the nuclear signal is the traffic through the destruction complex,
and an APC knockout removes it.

**What would be right.** Treated as a conserved pool of free and complex-bound β-catenin,
the signs come out right:
- APC knockout: no binding, all β-catenin free, targets UP;
- APC overexpression: more is bound, targets DOWN.

That takes two things:
- the pool (as specs/081);
- the released copy reading the pool's free state, because both copies are free β-catenin.
  If it read the release traffic instead, the inversion would remain.

## How common (`count_released.py`, canonical `20261020-2038_23ef3e1`)

The pattern: a single-protein entity whose starting copy is an input (or catalyst) of a
reaction, and whose released copy of the same entity is made by a reaction and carries on
downstream, joined by a route of at most 10 reactions through forms that carry the protein
(by reference entity).

**9 entities in 7 pathways**, route length in reactions:

| pathway | entity | reactions |
|---|---|---|
| Signaling by WNT | CTNNB1 | 8 |
| Signaling by Rho GTPases | SRGAP2, MKL1 | 2, 2 |
| Chromatin modifying enzymes | Me2sR3-H3 (×2) | 3 |
| TP53 | USP7 | 2 |
| L1CAM interactions | L1CAM | 3 |
| Mitotic G2-G2/M | CDC25B | 2 |
| DNA Repair | LIG1 | 3 |

## Reading

The pattern is rare, and β-catenin is by far the longest and best-curated case. A general
rule would touch 9 entities, most in pathways the benchmark scores little. The
alternative is a rule written for this one cycle, which would be fitted to one case. Not
registered until Adam chooses.
