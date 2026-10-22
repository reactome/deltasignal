# specs/089: a ubiquitination does not deplete its own ligase (curated active unit)

**Status:** pre-registered 2026-10-22, before the implementation, the build and the arms.

Adam, 2026-10-22: "yes" to a spec for the CKS1B drain, after specs/088 was not adopted
("we could revisit later as it has positive attributes").

## The defect (specs/077, 078, 088)

The generator's ubiquitin pass (`_emit_substrate_depletion_edges`) takes as a
ubiquitination's substrate every input leaf whose gene is also in the output. A ligase
that is carried through unchanged therefore counts as being ubiquitinated, and the
catalyst drains its own free subunits. In Mitotic G1:
- R-HSA-187575 "Ubiquitination of phospho-p27/p21" drains free CKS1B;
- R-HSA-1363331 "Ubiquitination of p130 (RBL2) by SCF (Skp2)" drains free CKS1B, and
  releases CUL1:SKP1:SKP2:CKS1B unchanged.

specs/078 removed these drains by dropping every leaf that leaves the reaction unchanged
(`LNG_UB_SUBSTRATE_CHANGED`, with or without `_UBX`). That fixed all 4 CCNE1 experimental
cases, but also removed the drains on the substrate's PARTNERS (EGFR, GRB2 under CBL),
through which curators' CBL-knockout answers are right: held-out −44 (EGFR −29, MET −10).
specs/088 (`fed_single_ub`) fixed the MDM2 knockout but lost 6 Mitotic G1 cases to this
drain.

## The rule (generator `LNG_UB_LIGASE_UNIT`, default off)

Only the **ligase** is taken out of the substrate set, not every carried partner. A
substrate leaf is dropped when its own stId is also an output leaf (it is not itself
modified, so autoubiquitination is untouched) and either:

- **A:** it is a leaf of the catalyst activity's curated `activeUnit` (Reactome marks the
  ligase inside the catalyst complex: CUL1:SKP1:SKP2:CKS1B in R-HSA-187575); or
- **B:** the catalyst has no `activeUnit`, and the leaf belongs to an output entity that
  is also a direct component of an input complex: the ligase released unchanged
  (CUL1:SKP1:SKP2:CKS1B in R-HSA-1363331).

Everything else in the ubiquitin pass is unchanged.

## Structural count (`census.py` → `census.tsv`, Neo4j Release 97, before the implementation)

102 ubiquitinations in the catalog pathways with an identified substrate. The rule drops
ligase leaves in **42** (A: 36, B: 6), in 28 pathways. For example:
- Mitotic G1: CKS1B (R-HSA-187575, A); CKS1B, CUL1, SKP1, SKP2 (R-HSA-1363331, B).
- TP53 and Cell Cycle Checkpoints: p-MDM2 and MDM4 in "MDM2 ubiquitinates TP53"
  (R-HSA-6804879, B).
- CBL: only p-Y-CBL in "CBL monoubiquitinates activated MET" (R-HSA-8875183, A). In the EGFR
  CBL reactions CBL is itself modified, so nothing is dropped there, and the drains on
  EGFR and GRB2 are kept in every CBL reaction.
- RUNX2: STUB1, SMURF1, WWP1, SCF(SKP2) (A), the class specs/078 saw as a threshold shift.
- TGF-beta (8 reactions), WNT (5), Hedgehog (4), DSB/Fanconi (4).

The number of network depletion edges removed is counted by diffing the variant build
against canonical, and committed, before any arm.

## Arms

| arm | build | solver | gated |
|---|---|---|---|
| `ulu` | variant `--env LNG_UB_LIGASE_UNIT=1` | main defaults | **yes** |
| `ulu88` | the same | `--server DS_CARRIER_LOOP_OPEN=fed_single_ub` (specs/088) | yes, separately |

Both against `results/e4e7192/dflt086b` on canonical `20261021-2218_f99c4ab`, read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`, case by case.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below
−10 (RAF set aside). If `ulu` passes, the flag is recommended as the generator default; if
`ulu88` also passes and beats `ulu`, specs/088 is recommended with it. Adam decides.

## Predictions

1. `ulu`: at least 2 of the 4 CCNE1 experimental cases specs/078 fixed are fixed.
2. `ulu`: EGFR and MET together move by no more than 5 curator cases (net), because the CBL
   drains on the receptor complex are kept.
3. `ulu88`: Mitotic G1 does not lose the 6 curator / 2 experimental cases specs/088 lost,
   and both MDM2-knockout experimental cases stay fixed.
4. No prediction for the held-out net (RUNX2, TGF-beta and WNT move in ways not traced).

## Network count (after the build, before any arm)

Variant build `20261022-1541_d31f95e_ulu` (LNG `d31f95e`, branch `ub-ligase-unit-089`)
against canonical `20261021-2218_f99c4ab`, by entity
(`specs/078-ub-ligase-subunits/count_depletion.py`):
- depletion edges 727 → 654: **73 removed**, none added, in 18 pathways (TGF-beta 27, WNT 9,
  NOTCH 6, RUNX2 6, Hedgehog 5, ERBB2 2, PIP3 2, TP53 2, others);
- every free CKS1B drain is removed: Mitotic G1 (from R-HSA-1363326 and R-HSA-187565) and
  S phase (from R-HSA-187565);
- the only other change is 24 `variant_pool` edges removed (ERBB2 6, Mitotic G1 8, S phase 8,
  Cell Cycle Checkpoints 2): depleter pools no longer needed. Nothing else moves.
