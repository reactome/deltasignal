# specs/095: a self-inhibitor built from another copy of the input

**Status:** pre-registered 2026-10-23, before the implementation, the census and the arm.

From specs/094: PTEN overexpression in PIP3 is the largest class of experimental errors on
which every reference agrees (9 cases, all read 83.5× UP, expected DOWN). specs/094 refuted
the loop-elasticity explanation. This spec traces it to something else.

## The trace (`pten*.jl`, canonical `20261022-2011_233665c`, current defaults)

The bench pins five PTEN nodes at 80×: the gene, the root copy of PTEN mRNA (R-HSA-2318746,
`08cc6d72`), two root copies of PTEN protein, and K48polyUb-PTEN. Then:

1. The PTEN protein that dephosphorylates PIP3 is a third, **translated** copy
   (R-HSA-199420, `57a5dff3`). It reads **0**, so PIP3 reads 9× (the de-repression
   ceiling), AKT 100×, and the readout p-MDM2 83.5×.
2. It reads 0 because its translation ("PTEN mRNA translation is negatively regulated by
   microRNAs", R-HSA-8944497) reads 0. Its input is the pinned root mRNA copy (80×), but it has
   **12 negative regulators**, each a PTEN mRNA:miRNA(:RISC) complex, each at 64.5×.
3. Those complexes are built from a **second copy of PTEN mRNA** (`fbc1de1c`), the product of
   "PTEN gene transcription …" (R-HSA-8944104), whose input is the pinned PTEN gene.
4. specs/022 (`DS_SELF_INHIBITOR_WEIGHT`) exists for an inhibitor built from its own reaction's
   input, and specs/040 rule B for one sharing a leaf with it. Neither fires. The rule requires
   the input node itself to reach the inhibitor, or a common carrier of the leaf to reach both,
   and the pinned root copy reaches only its own translation. So the 12 complexes act at full
   strength and translation goes to 0.

With `DS_SELF_INHIBITOR_WEIGHT=off` the case reads 0.17× (DOWN, right), but only because
transcription then collapses under its own 5 TF:gene inhibitors, so the second mRNA copy is 0.
That is not a mechanism.

Biologically the 12 complexes are PTEN mRNA bound by miRNAs. More PTEN mRNA should mean more
translation, which is the specs/022 reading: the inhibitor's change that its shared input
explains is counted at weight w/n.

## The rule (`DS_SELF_INHIBITOR_COPIES`, default `off`; `on`)

In the specs/022 pair test, an inhibitor I of target t that contains the entity of one of t's
inputs A (containment table) is a self-inhibitor when A reaches I, **or, under `on`, when another
node with A's stable id (another copy of the same entity) reaches I**. The damping is unchanged:
specs/022's weight, shared over the reaction's n self-inhibitor slots, weaken-only. The solve
reports `self_inhibitor_copies` and how many slots the copy clause added. Rule B's leaf clause is
unchanged.

## Census (after the implementation, before the arm)

`census.jl` → `census.txt`, the solver's own pair test on canonical: the copy clause adds **115
inhibitor slots in 5 pathways**, all new (none widens an existing slot): R-HSA-446728 60, TP53
transcription (R-HSA-3700989) 24, R-HSA-2559583 18, **PIP3 (R-HSA-1257604) 12, exactly the 12
PTEN mRNA:miRNA complexes**, RAF 1.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `sic` | canonical `20261022-2011_233665c` | `--server DS_SELF_INHIBITOR_COPIES=on` | `results/8162f5f/dflt091` |

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

## Predictions

1. PTEN overexpression in PIP3: at least 6 of the 9 experimental cases are fixed.
2. TP53 in PIP3 (4 experimental cases, TP53 drives PTEN transcription): no prediction.
3. No prediction for held-out or tuning nets.
