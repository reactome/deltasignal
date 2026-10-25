# specs/104: a pool's drive does not read another copy of its own entity (pre-registration)

Adam, 2026-10-25: "I want you to figure out the last two". Found while tracing the agreed WNT errors
after specs/103 failed its gate.

## Why (traced, canonical `20261024-2211_7f4b28a`, APC overexpression in Signaling by WNT)

APC and AMER1 overexpression read the ubiquitinated p-β-catenin:SCF(β-TrCP1) complex (R-HSA-2130284)
at 0.014. Curators and experiments both say UP. The readout's route leaves the β-catenin release
pool (specs/086) at its last intermediate, p-S33,S37,T41,S45 CTNNB1:p-AXIN:…:phospho-ub-APC
(R-HSA-195297), which reads 0.014 while every other intermediate reads 1.11.

- An intermediate reads its source fold times the drives of every step except its exit (specs/039
  amendment 2). R-HSA-195297 is the only one whose product leaves out step 7, PP2A
  dephosphorylation (R-HSA-3601585), so 0.014 / 1.11 = 1/80 is step 7's drive: 80×.
- Step 7 has 40 reaction copies. 20 take the pool's own intermediate. The other 20 take root copies
  of the same entity, R-HSA-195297 copies with no in-edge (variant combinations nothing produces).
  They contain APC, so they read 80× under APC overexpression.
- `path_drives` holds the pool's own nodes at baseline, since a drive is a rate per unit of source.
  It does not hold another copy of the same entity, so the overexpressed copy enters the drive as a
  rate. The pool then reads the exit as 80× faster and empties the intermediate the readout leaves
  from.

Footprint: 7,662 step-reaction inputs or catalysts in 9 pathways are non-pool copies of a pool entity
(EGFR 5,166 edges, IFN 936, DNA repair 546, RAF 480, RHO GTPases 392, Mitotic 112, WNT 20, S Phase 6,
Mitotic G1 4).

## The rule

Solver flag `DS_POOL_DRIVE_TWINS` = `off` (default, today's solver) | `hold`. Under `hold`, a node
that is an activator (input or catalyst) of a pool's step reaction, is not itself a pool node, and is
a copy of an entity the pool holds as a state or intermediate (same stable id before `::`) is held
at baseline while the pool's drives are evaluated, as the pool's own nodes already are. Nothing
else changes: such nodes keep their solved values everywhere else.

## Gate (as specs/077)

Curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10; read with RAF set aside too.

## Predictions

1. APC and AMER1 overexpression: R-HSA-195297 and the readout R-HSA-2130284 read between 0.85 and 1.5
   (no longer 0.014). Since the β-catenin pool is conserved, the ubiquitination flux stays near the
   release flux, so these two agreed cases go from DOWN to NORMAL: still wrong (expected UP).
2. Experimental net between −2 and +2. The rule is a faithfulness fix; it is not expected to close the gap.
3. Only the 9 listed pathways move.
