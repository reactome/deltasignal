# specs/183 research (2026-11-27)

Registration `5f200052`. `bench/analysis/183/pakt_183.py` on the Release 98 predictions and both Lun halves (FIT).
Counts here are over every measured protein (635 at p-AKT), not only the pairs the model scores, so they are larger
than specs/181's 106. Output outside the repo (`analysis/183/pakt.json`).

## Verdict

**The p-AKT fall is a general consequence of over-expressing an active protein, kinase or phosphatase, and it is
outside the curated network.** No loop is named as a candidate.

## 1. Who lowers it

- 159 of 635 proteins lower p-AKT S473 (25%); 63 raise it.
- Kinases 128 of 540 (24%), **phosphatases 31 of 95 (33%)**. Tyrosine kinases 32 of 60.
- It tracks how active the protein is at the other readouts:

| other readouts the protein changes | proteins | lower p-AKT |
|---|---|---|
| 0 | 61 | 3% |
| 1–3 | 289 | 17% |
| 4–9 | 187 | 28% |
| 10 or more | 98 | 58% |

The proteins are unrelated by family: casein kinases, PKA subunits, DUSPs and PTPNs, receptor and SRC-family tyrosine
kinases, MAP4Ks, lipid kinases.

## 2. What moves with it (proteins that lower p-AKT against the rest)

| readout | lowers p-AKT: DOWN | rest: DOWN | lowers p-AKT: UP | rest: UP |
|---|---|---|---|---|
| p-GSK3β | 23% | 1% | 7% | 5% |
| p-4EBP1 | 20% | 3% | 0% | 0% |
| p-S6 | 23% | 16% | 6% | 3% |
| p-PDK1 | 6% | 0% | 13% | 6% |
| p-p70S6K | 12% | 3% | 34% | 20% |
| p-ERK1/2 | 13% | 7% | 29% | 17% |
| p-p38 | 3% | 0% | 41% | 21% |
| p-JNK | 1% | 0% | 37% | 16% |
| p-MEK1/2 | 2% | 0% | 28% | 10% |
| p-p90RSK | 3% | 3% | 28% | 8% |

- **The PI3K branch as a whole falling: not supported.** p-PDK1 does not fall and p-p70S6K rises more often than it
  falls. What falls with p-AKT is AKT's own substrate site (GSK3β S9) and 4EBP1: the fall is at AKT and passes
  downstream, so it is a real loss of AKT activity, not an antibody effect.
- **Feedback from a raised growth pathway: not distinguishable.** ERK and p70S6K rise more often in these proteins, but
  so do p38, JNK, MEK and RSK; every arm rises in active proteins. Most p-AKT falls (71%) come without an ERK rise.
- **A general cost of over-expression: supported.** Phosphatases lower p-AKT at least as often as kinases, the
  proteins share no family, and the rate follows activity from 3% to 58%.

## 3. It is real

74% of the falls repeat in the EGF well (117 of 159), and 108 of 159 are larger than 1.5 δ (54 above 3 δ).

## 4. What the model holds

Read from the prediction file (the tracer was not needed): of the 159 proteins, 106 are scored at p-AKT. 42 have a
curated route, and the model calls UP for 36 and DOWN for 5; 64 have no route. Where Reactome connects these proteins
to AKT it is as activators, which is the curated biology of a growth-factor receptor.

## Reading

- The stress and MAPK arms rise and AKT S473 falls together as a protein becomes active on over-expression. That is a
  cell-state response to the perturbation, not a pathway from the protein.
- A model of curated routes cannot produce it, and adding a feedback edge would be fitting one screen: the same edge
  would have to explain phosphatases.
- It is the same finding as specs/176 from the other side: a readout's own responsiveness and the protein's activity
  carry the screen.
