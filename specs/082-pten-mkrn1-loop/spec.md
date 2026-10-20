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
