# specs/057: a ubiquitination depletes what it changes, not what it carries

**Status:** registered 2026-10-14, before any build or arm. **Gate fails** (research.md R1: held-out −39, p 0.0011; experimental 0). Not adopted; the flag stays off. The MET fix works as predicted (+30). Losses: CBL/EPS15/SH3GL1 (−38) were right only through the removed edges; SMAD2/3 (−24) come from a defect in the registered test, found after the arm (Reactome's "Ub-X" complex products, 24 of 301 edges); RNF168 (−16) is a non-converged coin.

## Why

Found by tracing the held-out gap (Adam, 2026-10-14: reach the 90s on the
curator axis by fixing what fails inside each pathway). In Signaling by MET,
the MET knockout reads every readout UP (20–100×) and the overexpression
reads every one DOWN, 20 of 20 wrong; the HGF knockout reads UP, 10 of 10
wrong. The response to the MET pin is smooth and converged (no switch, unlike
RAF in specs/054). Traced on the canonical build (`20261012-1401_01206ca`):

1. **The only MET root pinned is p-Y1003,4Y-MET** (R-HSA-8874689), the form
   that enters the CBL branch. Plain MET (R-HSA-419603) is not a root: its
   one in-edge is a depletion from "unknown ubiquitin ligase" (R-HSA-5250898).
2. **That depletion comes from "LRIG1 undergoes ubiquitination"**
   (R-HSA-8875431): input MET:LRIG1, output MET:Ub-LRIG1. LRIG1 is
   ubiquitinated; MET passes through unchanged.
3. **The CBL branch depletes GRB2.** "CBL monoubiquitinates activated MET"
   (R-HSA-8875183) takes HGF:p-4Y-MET dimer:GRB2:p-Y-CBL and returns it
   with MET mono-ubiquitinated; GRB2 is unchanged. The generator emits a
   depletion from the p-Y-CBL complex pool to free GRB2 (R-HSA-74686). So
   more CBL-branch MET lowers GRB2 (0.51 at MET 2×) and every readout through
   GAB1:GRB2; a knockout of it, or of HGF, de-represses GRB2 to the 10×
   ceiling, which drives the MET recycling loop up even with MET at 0.

The cause in the generator (`_emit_substrate_depletion_edges`, ubiquitin
pass): the substrate is "an input leaf whose gene also appears in the
output". When a complex is modified in place, every partner's gene appears in
the output, so every partner is a "substrate".

Trace, editing only the bundle (not an arm): removing the GRB2 edge makes the
HGF knockout read every readout ~0 (right, 10 cases); removing it and also
pinning plain MET makes both MET directions right (20 cases). Scripts:
`bench/analysis/057/`.

Verified in Neo4j (Release 97) for both reactions: the passed-through
partner's entity (same stId) is among the output's leaves.

## Census (before any build)

`bench/analysis/057/ub_census.py` on the canonical build, with
`ub_reactions.cypher`:
- 98 of the 126 ubiquitinations in catalog pathways count a passed-through
  partner as a substrate;
- **301 of the catalog's 727 depletion edges** target such a partner (and no
  changed substrate of another ubiquitination in the same pathway), in 32
  pathways: DNA Repair 50, TGF-beta 46, DSB Repair 36, EGFR 34, WNT 23,
  NOTCH1 13, Hedgehog 10, and 25 more;
- those pathways hold 5,886 of 18,481 valid held-out cases (21 pathways),
  10 tuning pathways, and 820 of 845 experimental cases.

The approximation: leaves are read to depth 3 (`hasComponent|hasMember|
hasCandidate*0..3`), as the generator reads them.

## The rule (generator, `LNG_UB_SUBSTRATE_CHANGED`, default off)

A ubiquitination's substrate is an input leaf whose gene appears in the
output **and whose own stId does not**: it leaves the reaction in a changed
form. A partner whose exact entity is among the output's leaves is carried,
not depleted. Nothing else changes: the phosphatase pass, the pooling of a
catalyst's variants, and every guard are as before. With the flag off the
generator is byte-identical to `01206ca` (asserted by a test).

## Arms

| arm | build | solver | against |
|---|---|---|---|
| control | canonical `20261012-1401_01206ca` (LNG `01206ca`) | defaults (`results/64d23cb`) | — |
| `ubs` | variant build `--env LNG_UB_SUBSTRATE_CHANGED=1` from the LNG branch | defaults, this branch | control |

The solver (`src/`) is identical to `64d23cb`; `bench/` differs from it only
by specs/056's `--pair-on-valid` and context code, which a build without
context bridges does not exercise. Both axes. Pin sets change by design (MET
gains a root), so pairs are taken with `holdout_report.py --pair-on-valid`
and the cases valid in one arm only are reported separately.

## Gate (the specs/049 amendment 5 faithfulness gate, RAF set aside per specs/054 D1)

Against the control:
- curator held-out net > 0 at McNemar p < 0.01;
- experimental net loss (RAF set aside) not significant (p ≥ 0.05);
- no pathway (RAF set aside) losing more than 10 cases on either axis;
- every experimental loss outside RAF traced to a named mechanism.

Both readings reported (with RAF and without). Held-out is SELECT
(specs/051). If the gate passes, the flag is recommended as the generator
default (Adam's call), and a new canonical build follows. If it fails,
nothing changes, and the losses are traced.

## Predictions (written before the build)

- MET: the 20 MET cases and the 10 HGF knockout cases go right (30).
- The 301 edges are spread over 32 pathways; I do not predict a sign for the
  rest. The risk is that some of these edges carry real biology by accident
  (a partner whose free pool genuinely falls when the complex is degraded);
  specs/048's history (topology ubiquitin, 2026-05-31) says ubiquitin
  depletion helped the experimental axis when it landed.

## Out of scope

- A pin rule for a form whose only in-edges are depletion (28 gene–pathway
  pairs in the canonical build). Under this rule MET's case dissolves (its
  depletion edge is one of the 301), so it is measured, if at all, after this.
- Degradation as such: what a ubiquitination does to the changed substrate
  is unchanged.

## Amendment 1 (2026-10-14, after R1, before any amended build): the "Ub-X" complex test

R1 found the registered test wrong for Reactome's "Ub-X" complex products
(e.g. Ub-SMAD2 [cytosol] = SMAD2 + Ub): the real substrate's stId reappears
in the output, so 24 of the 301 removed edges hit a real substrate.

**Rule (generator `LNG_UB_SUBSTRATE_UBX`, default off, read only with
`LNG_UB_SUBSTRATE_CHANGED=1`):** an input leaf whose gene is in the output is
the substrate if its stId is absent from the output's leaves, **or** it is a
direct component of an output complex (to depth 3) that also has Ub as a
direct component.

**Not registered: a degradation rule for the CBL class.** Proposed in R1, but
checked before registering: in "CBL monoubiquitinates activated MET" the CBL
complex is both input and catalyst, so "the catalyst depletes the input
complex holding the substrate" is a self-depletion the generator skips. A CBL
knockout could only raise signaling if CBL binding drained the active
receptor, the input-consumption rule refuted on 2026-07-16. Not pursued.

**Arm `ubx`:** variant build `--env LNG_UB_SUBSTRATE_CHANGED=1 --env
LNG_UB_SUBSTRATE_UBX=1`, same solver, against the canonical results; same
gate, same pairing. Also reported against `ubs`.

**Predictions:** the 24 edges return; SMAD2/3 in TGF-beta recover (about
+24 against `ubs`); MET stays +30; CBL/EPS15/SH3GL1 (−38) and RNF168 (a coin)
do not change. Against canonical: about −15 held-out, so **the gate most
likely still fails**. If it fails, nothing is adopted; the corrected test
stays on record as the faithful form for any future degradation work.
