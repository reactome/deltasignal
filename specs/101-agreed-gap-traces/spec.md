# specs/101: balanced variant choices (generator `LNG_BALANCE_CHOICES`)

**Status:** pre-registered 2026-10-24, before the implementation, the build and the arm.

Adam, 2026-10-24: "keep going" (closing the experimental gap through agreed errors).

## The defect (trace T2 in research.md; `trace_mdm4.jl`, canonical `20261024-0801_7a754fc`)

MDM4 knockout in Transcriptional Regulation by TP53 (R-HSA-3700989) reads TP53 only 1.06–1.10×
(curator 8 of 51 right; experimental 3 cases wrong, all UP expected). TP53 is depleted by the
MDM2 ligase pool, and the MDM2:MDM4 ligase reads 0.80× with MDM4 at zero.

Cause: "MDM2 forms homo- or heterodimers" (R-HSA-6804741) is curated as input set
{MDM2, MDM4} (R-HSA-6804750) + p-MDM2 ×2 → output set {MDM2 dimer, MDM2:MDM4} (R-HSA-6804745).
The output set's members have other reference signatures than the input set's, so
`bind_output_slots` does not pair them; the output set is a free slot and the generator emits
the cross product: four copies, two of them unbalanced:
- p-MDM2 ×3 → MDM2:MDM4 (MDM4 made from nothing; feeds the heterodimer 1.6× under MDM4 KO);
- MDM4 + p-MDM2 ×2 → MDM2 dimer (MDM4 lost).

Census (`balance_census.py`): 2,207 output edges in 85 reactions, 37 of 92 pathways, put a protein
into an output that the copy's inputs do not contain while another copy of the same reaction does
supply it (VEGF homodimerization, NOTCH cleavage, SMAD2/3:SMAD4, PAK activation, ...).

## The rule (generator, default off)

`LNG_BALANCE_CHOICES=1`: for each reaction's copies (choices), score each copy
d = |proteins its outputs contain that its inputs do not| + |proteins its inputs contain that its
outputs do not|, proteins compared by reference entity, inputs only (not catalysts or
regulators). Drop a copy when BOTH a copy with the same output variants AND a copy with the same
input-side choice score lower. A copy whose input choice has no better-balanced alternative is
kept (incomplete curation keeps its route). Reactions with no free output slot are untouched.

## Arms

Variant build with `LNG_BALANCE_CHOICES=1` from the generator branch; arm `bal` (solver defaults)
against `dflt100` (canonical, the control is the canonical build itself; a `--variant` control
build with the flag off must equal canonical).

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.
Both readings (RAF kept / set aside) reported.

## Predictions

1. Under MDM4 knockout the MDM2:MDM4 ligase node reads 0 (only the MDM4 copy feeds it).
2. MDM4 knockout in TP53: curator at least 10 of the 43 wrong cases fixed; experimental at least
   1 of the 3.
3. Only the 37 census pathways can move.
4. No prediction of sign for the other pathways (held-out mostly).
