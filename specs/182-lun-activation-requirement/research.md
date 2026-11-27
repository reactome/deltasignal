# specs/182 research (2026-11-27)

Registration `d6ade545`. Classes from build `20261126-2322_d6d026b_top98` (`bench/analysis/182/activation_classes.tsv`,
structure only); predictions `results/906698aa/lun178_r98`; both Lun halves are FIT. Output outside the repo
(`analysis/182/act.json`).

## Verdict

**The registered criterion is met, and the result is not informative.** Class C holds 7 routed pairs. Read as a
whole, structure does not tell an over-expressed protein that acts from one that does not.

## Classes

Of the 428 scored proteins: A (the starting form is a catalyst or regulator) 287, B (binds or modifies itself first)
110, C (another enzyme first) 24, N 7. Pooling over networks puts a protein in A if it is A anywhere, which leaves C
almost empty.

## Registered measurements (pooled, pairs with a curated route)

| class | pairs | proteins | readout changed | model calls | calls right |
|---|---|---|---|---|---|
| A | 1,128 | 149 | 25% | 757 | 18% |
| B | 309 | 37 | 41% | 214 | 34% |
| C | 7 | 7 | 0% | 7 | 0% |
| N | 43 | 4 | 53% | 29 | 48% |

- **A − C: +0.18 [+0.13, +0.23]**, positive on each half (+0.21, +0.16), so the registered criterion passes. It
  rests on 7 calls, all wrong.
- Withholding class-C calls: 7 fixed, 0 broken, recall unchanged. It passes specs/161's rule and removes 7 of the 699
  false calls.
- **The large contrast goes the other way.** A, where Reactome says the starting form acts, is the worst class and
  holds three quarters of the calls; A − B is −0.16 [−0.28, −0.05].
- That contrast is the tyrosine kinases again (class B is where receptors that bind a ligand and phosphorylate
  themselves fall). Without tyrosine kinases: A 10%, B 19% pooled, and on half B alone A 9%, B 5%.
- Against measured engagement (specs/179): A 23 engaged / 42 not; B 13 / 6; C has no protein with a measured direct
  pair. A is engaged less often than B.

## Post hoc: the class taken in the network the route is in

Not registered. This fills class C (41 pairs, 23 proteins) and the signal disappears: calls right A 18%, B 25%,
C 25%; A − C −0.07 [−0.21, +0.11]. Without tyrosine kinases 11%, 13%, 14%.

## Reading

- Whether Reactome draws a protein's starting form as a catalyst says nothing about whether over-expressing it moves
  anything in HEK293T. Most kinases are drawn that way, and most of them do nothing here.
- With specs/181, every structural handle tried on this screen is now negative: more edges, route length, predicted
  magnitude, and activation class. What separates right from wrong is measured (specs/176 protein activity, specs/179
  first-hop engagement) or the tyrosine-kinase class (specs/177).
- The practical consequence for an over-expression query is a statement, not a rule: the answer assumes the protein is
  active, and a measured first hop is what tells you whether it is.
