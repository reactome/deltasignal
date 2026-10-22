# specs/088: ubiquitin is not a displacer (amends specs/077's fed_single)

**Status:** pre-registered 2026-10-22, before the implementation and the arm.

Adam, 2026-10-22: "ok" to the rule.

## The defect (traced, `mdm2.jl`, canonical `20261021-2218_f99c4ab`)

MDM2 knockout in Cell Cycle Checkpoints (R-HSA-69620) reads UP at 100×, the same as
overexpression, at PolyUb-TP53 Tetramer (R-HSA-3209186) and at
PolyUb,p-MDM2:PolyUb,p-MDM4 (R-HSA-6804937). The latter is a complex containing MDM2. The
experiment and the curators say DOWN. Sweeping the MDM2 pins (p-S166,S188-MDM2,
R-HSA-6793669) gives:

| MDM2 fold | 0 | 0.5 | 0.9 | 0.99 | 1 | 1.01 | 1.1 | 2 | 80 |
|---|---|---|---|---|---|---|---|---|---|
| readout | 100 | 1.0 | 0.83 | 0.98 | 1 | 1.02 | 1.2 | 100 | 100 |

Under knockout:
- "MDM2 forms homo- or heterodimers" (R-HSA-6804741) makes p-MDM2:MDM4 (R-HSA-6804932)
  at 0, as it should.
- p-MDM2:MDM4 is also the carrier released by "MDM2 ubiquitinates TP53"
  (R-HSA-6804879: MDM2:MDM4:TP53 + Ub → PolyUb-TP53 + MDM2:MDM4).
- The ubiquitination's inhibitor, p14ARF:p-MDM2 (R-HSA-6804995), contains MDM2 and falls to 0.
  The step is de-repressed to the cap, and the release feeds p-MDM2:MDM4 at 33× with no MDM2
  present.

specs/077 drops exactly this kind of release, a carrier whose free form has another
producer. But `fed_single` (amendment 1) keeps a release whose reaction has a second
non-cofactor input, as a displacer, and here that input is **ubiquitin**. Ubiquitin is a
donor co-substrate, not a displacer. The generator already treats it that way when it
finds pool steps (`r_steps`, `_UBIQUITIN_STIDS`).

## Structural count (`count_ub_displacers.py`)

Of the 23 releases `fed_single` keeps as displacements on canonical:
- **10 are kept only because of ubiquitin:**
  - EGFR 7: CBL complexes released after ubiquitinating the receptor (R-HSA-182928);
  - Cell Cycle Checkpoints 2: p-MDM2:MDM4 and the MDM2 dimer;
  - Mitotic G1 1: SCF(Skp2) released after ubiquitinating p130.
- **13 have another displacer:** R-HSA-168928 12, intrinsic apoptosis 1 (BAD).

## The rule

`DS_CARRIER_LOOP_OPEN=fed_single_ub`: as `fed_single`, except that the ubiquitin entities
are not counted as a release reaction's second input. These are R-HSA-68524 (nucleoplasm),
R-HSA-113595 (cytosol) and R-HSA-9660007 (lysosomal lumen), the generator's list. The
solve reports the rule in `carrier_loop_open_rule`.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `fsu` | canonical `20261021-2218_f99c4ab` | `--server DS_CARRIER_LOOP_OPEN=fed_single_ub` | `results/e4e7192/dflt086b` |

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10 (RAF set aside).

## Predictions

1. MDM2 knockout in Cell Cycle Checkpoints: the experimental cases at PolyUb-TP53 and
   PolyUb-MDM2:MDM4 no longer read 100× (at least 1 of the 2 fixed).
2. Only EGFR, Cell Cycle Checkpoints and Mitotic G1 move.
3. No prediction for EGFR or held-out. The CBL releases are the class specs/057 found right
   for the wrong reasons.
