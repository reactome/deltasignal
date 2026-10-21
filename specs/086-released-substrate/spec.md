# specs/086: a substrate released from the complex that degrades it (analysis, design pending)

**Status:** analysis and rule 2026-10-21. Arm `rlp` failed its predictions (1 of 20 variant chains). Amendment 1 (all chains, post hoc), arm `rlp2`: **gate passes**, held-out 0, tuning +73 (84 / 11), experimental +8 (10 / 2); APC 5 of 7 fixed. Adoption is Adam's.

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

## The rule (pre-registered 2026-10-21, before the implementation and the arm)

Adam, 2026-10-21: "ok" to a general rule for the pattern.

**Generator `LNG_RELEASE_POOLS` (default off; LNG branch from main `43d7919`).** A
*release pool* is found when all of these hold:
- **S:** a starting copy (no in-edge) of a single protein R, which is an `input` or
  `catalyst` of a reaction.
- **P:** another copy of the same entity, made by a reaction, which has an out-edge (it
  carries on; this is what separates it from specs/081, where P leads nowhere).
- **Route:** S reaches P through at most 10 reactions, via forms that carry R by reference
  entity.
- **Modified on the way:** some form on the route carries R with a different modification
  signature from S. A partner released unchanged is a carrier (specs/077) and is excluded.

What it writes:
- **States:** M, the form the release reaction consumes, is the base (at rest the protein is
  mostly bound). S is the second state.
- **Intermediates:** the forms between S and M.
- **Paths:** forward S → … → M through the route's reaction copies; return M → S through
  the release reaction.
- **Carrier:** P is listed as the pool's carrier, with the release reaction copies as its
  release steps. Under `DS_CYCLE_CARRIERS=on` P then reads its other producers, not the
  release traffic.
- **Edges:**
  - a `released_twin` OR edge S → P, so P's other producer is S and P reads the pool's free
    state;
  - a structural `recycle` edge P → S (specs/081), so the cycle is one component.
- A candidate any of whose nodes an existing pool already claims is skipped and counted.

**DeltaSignal:** `CYCLE_MAX_STEPS` goes from 6 to 10 (the forward path for β-catenin is 7
steps). The generator caps ordinary pools at 6, so no existing pool table has a longer path
and canonical is unaffected; this is checked by the arm's baseline. No other solver change.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `rlp` | variant `--env LNG_RELEASE_POOLS=1` | main plus the `CYCLE_MAX_STEPS` change | that solver on canonical (`dflt086`) |

The variant must equal canonical except in the pathways the count lists; checked by
`compare_build.py` before the arm.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10 (RAF set aside).

## Predictions

1. APC in WNT: at least 5 of the 7 experimental cases fixed (knockout UP, overexpression
   DOWN at β-catenin and the canonical targets).
2. Signaling by WNT curator net at least +15.
3. No pathway other than the 7 counted moves.
4. No prediction for the held-out net: the other 8 entities include held-out pathways
   whose cases we have not traced.

## Amendment 1 (2026-10-21, after arm `rlp`, before its arm): every variant chain in one pool

Adam, 2026-10-21: "ok". **Post hoc:** written after `rlp` showed the pool captured 1 of the
20 destruction-complex variant chains (`research.md`).

- **Rule.** As registered, but for each (S, P) the pool holds every release copy of the
  release reaction that makes P:
  - each copy's consumed complex M_v is a base state, as specs/081 keeps p130's four resting copies;
  - each M_v gets its own shortest route from S (forward path) and its own release (return path);
  - intermediates are the union of the routes' forms;
  - P is the carrier, with all the release copies as its release steps.
- **Arm `rlp2`:** a variant rebuilt from the amended generator, against `dflt086`, with the
  same gate.
- **Predictions:**
  1. APC: at least 5 of the 7 experimental cases fixed.
  2. Signaling by WNT curator net at least +10.
  3. Nothing outside WNT and L1CAM moves.
