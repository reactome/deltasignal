# specs/100: a root copy of an RNA reads its transcribed copies

**Status:** pre-registered 2026-10-23, before the implementation, the build and the arm.

Adam, 2026-10-23: "ok" to closing the experimental gap with MP-BioPath through agreed errors;
"ok" to tracing TP53 in PIP3.

## The defect (`tp53pip3.jl`, canonical `20261022-2011_233665c`)

TP53 in PIP3 activates AKT signaling (R-HSA-1257604) is wrong on **all 28 curator cases and 4
experimental cases**: TP53 knockout reads DOWN (expected UP), overexpression reads UP (83.5–100×,
expected DOWN). TP53 acts through PTEN transcription, and PTEN mRNA (R-HSA-2318746) has two node
copies:
- a **root** copy, the input of translation ("PTEN mRNA translation is negatively regulated by
  microRNAs", R-HSA-8944497), which makes the PTEN that dephosphorylates PIP3;
- a **transcribed** copy, the output of "PTEN gene transcription …" (R-HSA-8944104, which TP53
  stimulates), which feeds only the 12 PTEN mRNA:miRNA complexes that repress that translation.

| TP53 | root mRNA | transcribed mRNA | miRNA complexes | translation | PTEN on PIP3 | PIP3 | readouts |
|---|---|---|---|---|---|---|---|
| 1× | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 80× | 1 | 80 | 80 | 0 | 0 | 9 | 100, 50, 100 |
| 0 | 1 | 0 | 0 | 10 | 7.7 | 0.13 | 0, 0.5, 0 |

The root copy never follows TP53, so TP53 reaches AKT only through the repressors, with the sign
reversed. specs/095's damping cannot help: the input it compares against (the root copy) does not
move. One mRNA is one pool, whichever reaction made a node copy of it.

## Census (`census.py` → `census.txt`, Neo4j referenceEntity class)

Entities with a root copy and a produced copy of one stable id in one pathway: 444 protein pairs
in 46 pathways (the copy-joining direction specs/053 and 085 refuted), but **6 RNA pairs in 4
pathways**: PTEN mRNA in PIP3 and in TP53 transcription (R-HSA-3700989), R-HSA-3209107 in
R-HSA-2559583, three in R-HSA-446728.

## The rule

- **Generator:** a table `rnas.csv` (stable ids of entities whose reference entity is a
  `ReferenceRNASequence`), written for every pathway, as `drugs.csv` and `pathogens.csv` are.
- **Solver** (`DS_RNA_TWIN`, default `off`; `on`): a root copy (no in-edge) of an RNA entity in that
  table gains an OR in-edge (`rna_supply`) from each produced copy (one with an `output` in-edge)
  of the same stable id. One-way, as `recycled_twin` (specs/059) and `closed_twin` (specs/049).
  The edge exists only in the solver: the benchmark's root detection is unchanged, so a root copy
  the protocol pins stays pinned (its pin wins). Reports `rna_twin_rule` and `rna_twin_edges`.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `rnt` | variant build with `rnas.csv` (otherwise canonical) | `--server DS_RNA_TWIN=on` | the same build with main's defaults |

The variant build must equal canonical in every edge before the arm.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

## Predictions

1. TP53 in PIP3: at least 20 of the 28 curator cases and at least 3 of the 4 experimental cases
   are fixed.
2. Only the 4 pathways with an RNA pair can move.
3. PTEN overexpression in PIP3 (specs/095's 9 experimental cases) stays right: its root mRNA copy
   is pinned.
4. No prediction for held-out (TP53 transcription R-HSA-3700989, R-HSA-2559583 and R-HSA-446728 are
   held-out pathways).


## Build (before any arm)

The first variant build (`20261023-2301_998b4b9_rna`) listed, besides 135 RNA species, 2,393
variants of RNA-containing complexes (miR-26A RISC): the shared table writer counts a variant
that chose a listed member as listed, right for drugs, wrong here. Fixed in the generator (LNG
`135af3c`, test fails on the old code) so only nodes that ARE an RNA species are listed. Rebuilt:
`20261023-2341_135af3c_rna2` equals canonical in every edge (0 of 92 differ); `rnas.csv` lists 135
RNA species in 19 pathways, no complex.

Arms: `rnd` (main's defaults on this build, the control) and `rnt` (`--server DS_RNA_TWIN=on`).
