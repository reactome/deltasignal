# specs/082: PTEN overexpression in PIP3 is decided by a knife-edge loop (analysis)

**Status:** analysis, 2026-10-20. No rule is proposed or measured.

From specs/079: the largest group of experimental errors on which the curators and the
experiment agree is PTEN overexpression in PIP3 activates AKT signaling (R-HSA-1257604).
That is 9 cases, all expected DOWN, and all predicted UP at 83.5×.

## What the solver does (`pten.jl`, `pten3.jl`; canonical build, current defaults)

The PTEN pins (gene, mRNA, protein ×2, K48polyUb-PTEN) set to a range of levels; readout
p-T,p-S-AKT [cytosol]:

| PTEN | 0 | 0.5 | 1.0 | 1.1 | 2 | 10 | 80 |
|---|---|---|---|---|---|---|---|
| readout | 100 | 0 | 1 | 100 | 100 | 100 | 100 |

A 1% rise (1.01×) already rails it. The readout is a switch, and its direction does not
follow PTEN: 0.5× reads 0 and 1.1× reads 100.

## The loop (Neo4j Release 97)

| step | edge | value at PTEN 1.01× |
|---|---|---|
| p-AKT1 catalyses "AKT phosphorylates MKRN1" (R-HSA-8948757) | catalyst | p-AKT1 50.5× |
| p-S109-MKRN1 (R-HSA-8948758), an E3 ligase, ubiquitinates PTEN | depletion (generator ubiquitin pass; PTEN is the real substrate) | p-MKRN1 83.5× |
| free PTEN (R-HSA-199420, the translated copy) | | 0.05× |
| PTEN dephosphorylates PIP3 (phosphatase depletion) | depletion | PIP3 at the 10× de-repression ceiling |
| PIP3 → PDPK1 → p-AKT | | p-AKT 100× |

The loop is AKT → MKRN1 ⊣ PTEN ⊣ PIP3 → AKT, a double negative and so positive feedback.
It is curated. Each fold-multiplying hop has gain 1, so the loop's gain is exactly 1 at
baseline. That is the specs/014 knife-edge: any perturbation rails it to the 10× and 100×
caps, and which way it goes is not decided by the biology.

Other links: FRK's stabilising phosphorylation of PTEN (R-HSA-8847968, R-HSA-8847977)
collapses with the free PTEN, and the PTEN:TRIM27 depleter falls to 0.22×.

## Reading

This is the RAF class (specs/054): a curated positive feedback loop with no saturation
in the model, so its readouts are a switch. RAF was set aside from rule decisions. PIP3
is a tuning pathway, and these are 9 of the 101 experimental errors on which every
reference agrees. Any fix is a loop-gain question (specs/014's sigmoid-elasticity arms), not a
network defect. The network is faithful here. Not pursued further in this spec.

## How common (`switch_census.py results/3e96e7c/dflt`)

These are (pathway, gene, readout) triples whose knockout and overexpression read the same
direction off baseline (|fold − 1| > 0.15 both), so the readout is a switch rather than a
response:

| axis | triples | cases | wrong | share of all errors |
|---|---|---|---|---|
| curator | 79 | 158 | 107 | 3.9% of 2,763 |
| experimental | 13 | 26 | 13 | 5.3% of 244 |

Wrong cases by gene:
- **curator:** CBFB in RUNX2 23, PTEN 14, ATM 12 and CHEK2 12 in TP53, SH3GL1 in EGFR 6,
  RAF1 5, SRC 4, RBL2 4;
- **experimental:** PTEN 8, RAF1 3, MDM2 in Cell Cycle Checkpoints 2.

Each is a candidate for the same trace: find the loop, check it is curated, and see whether
its gain is 1.

**Rails vs U-shapes.** The census mixes two things:

- **Rails: the same value both ways** (PTEN 83.5 / 83.5). 28 wrong curator cases, of
  which PTEN 14, RAF1 4, ACVR2A 3, SH3GL1 3, TP53 in Cell Cycle Checkpoints 2. And 12 wrong
  experimental cases: PTEN 8, RAF1 3, MDM2 1. These are knife-edge loops.
- **U-shapes: two different values in the same direction.**
  - CBFB in RUNX2: knockout frees RUNX2, so RUNX2:SOX9, RUNX2:SATB2 and RUNX2:RBM14
    read 5.6×; overexpression reads 1.9×. Curators: NORMAL. Partner competition, not a loop.
  - ATM and CHEK2 in TP53: knockout reads 0.56× and overexpression 0.23× at TP53 targets,
    both converged. Curators: NORMAL.
  - Not traced here.

## The other rails, traced (2026-10-21; `rails.jl`, `sh3.jl`, canonical `20261020-2038_23ef3e1`, current defaults)

Readout fold under the gene's pins at a range of levels:

| case | 0 | 0.5 | 0.9 | 0.99 | 1.0 | 1.01 | 1.1 | 2 | 80 | shape |
|---|---|---|---|---|---|---|---|---|---|---|
| EGFR, SH3GL1 at R-HSA-180523 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | **knife-edge** |
| RAF, RAF1 at R-HSA-1268261 | 100 | 100 | 100 | 100 | 1 | 0.003 | 0.003 | 0.002 | 100 | **knife-edge** (RAF, set aside) |
| Activin, ACVR2A at R-HSA-173511 | 50 | 2.14 | 1.03 | 1.001 | 1 | 0.999 | 1.007 | 1.87 | 50 | U-shape |
| Cell Cycle Checkpoints, TP53 at R-HSA-6804937 | 100 | 1.55 | 1.05 | 1.003 | 1 | 0.997 | 0.98 | 100 | 100 | U-shape |

**Two different things:**
- **Knife-edges: PTEN, SH3GL1, RAF1.** Any change rails the readout.
  - SH3GL1: a 1% rise drops the EGF:p-EGFR dimer to 0.05×. Through the curated
    PIP3 → GAB1:PI3K → PIP3 amplification loop (R-HSA-177927 / R-HSA-177939, catalyst on
    PIP3's own production), the route collapses to 0.
  - Each is a curated positive loop of gain 1 at baseline (specs/014). In PTEN's and EGFR's,
    depletion edges close the loop, which specs/014's elasticity (activator edges only)
    does not reach.
- **U-shapes: ACVR2A, TP53 in Cell Cycle Checkpoints, and CBFB, ATM/CHEK2 above.** These
  are graded near baseline, so not a knife-edge, but both directions move the readout the
  same way. A different mechanism (partner competition is one case), not traced here.

**How much each is worth** (from the census, wrong cases):
- knife-edges: about 28 curator and 12 experimental, PTEN the bulk; RAF is set aside;
- U-shapes: about 79 curator and 1 experimental.
