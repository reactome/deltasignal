# specs/093: a release that moves the carrier across compartments is a return, not a let-go

**Status:** pre-registered 2026-10-23, before the implementation, the build and the arm.

Adam, 2026-10-23: "ok" to following the CBL/MET lead from specs/092.

## The defect (specs/092, traced)

specs/092 found 10 curator CBL cases in Signaling by MET (R-HSA-6806834) whose only route the
solver cuts. They were all right before specs/077 (20/20) and are 10/20 after. The cut is
specs/077's open-loop rule (`DS_CARRIER_LOOP_OPEN=fed_single`): the free form
HGF:p-Y1234,Y1235,Y1349,Y1356-MET dimer [plasma membrane] (R-HSA-6806977) has a producer
outside its loop (HGF activation), so it stops reading its own release. The dropped release is
**"RAB4:GTP promotes recycling of activated MET receptor to the plasma membrane"**
(R-HSA-8875659): the receptor comes back from the recycling endosome. That is not a partner
letting go of a carrier. It is the receptor returning to where it signals, a real supply that
CBL's ubiquitination competes with (degradation against recycling).

## Census (`census.py` → `census.txt`, Neo4j compartments, canonical `20261022-2011_233665c`)

Of 94 distinct (pathway, free form, release) open-loop pairs, **3** move the free form across
compartments (its compartment is not among the release's inputs'):

| pathway | free form | release |
|---|---|---|
| Signaling by MET (R-HSA-6806834) | HGF:p-MET dimer [plasma membrane] | RAB4:GTP recycles activated MET (from the recycling endosome membrane) |
| L1CAM interactions (R-HSA-373760) | L1CAM [plasma membrane] | Reinsertion of L1 into the plasma membrane (from the cytosol) |
| R-HSA-453274 | p-S435-GTSE1 [nucleoplasm] | GTSE1 translocates to the nucleus (from the cytosol) |

## The rule

- **Generator** (`carrier_loops.csv`, additive): a column `moves` is True for a loop row whose
  release reaction carries the free form across compartments, by the census's test (Neo4j
  compartments). Always written; it changes no edge.
- **Solver** (`DS_CARRIER_LOOP_MOVES`, default `drop`; `keep`): under `fed_single` (and
  `fed`, `fed_single_ub`), an open release with `moves` = True is kept, as a displacement is.
  A table without the column reads `moves` = False, so old builds solve as before. The solve
  reports `carrier_loop_moves` and how many releases it kept.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `mvk` | variant build with the `moves` column (generator branch, otherwise canonical) | `--server DS_CARRIER_LOOP_MOVES=keep` | the same build with main's defaults (`mvd`) |

The variant build must equal canonical in every edge (`compare_build.py`) before the arm; `mvd`
must reproduce canonical's defaults.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

## Predictions

1. Only Signaling by MET, L1CAM interactions and R-HSA-453274 move.
2. At least 6 of the 10 CBL cases in MET specs/092 lists are fixed.
3. No prediction for L1CAM or GTSE1.
