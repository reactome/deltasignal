# specs/110: a complex that modifies itself, counted once

**Status:** pre-registered 2026-10-26, before any arm.

Adam, 2026-10-26: after specs/109, work on the large loops, starting with complexes that modify themselves.
Rules are not chosen from experimental or curator traces (specs/106); they are judged by the loop-rule gate
(specs/109).

## What the traces show

Reactome often curates one complex as both the input and the catalyst of a reaction: the complex modifies or
processes itself. Checked in Neo4j: the TGFB1:p-TGFBR:ZFYVE9:SMAD2/3 complex phosphorylating its SMAD2/3
(R-HSA-170868), p-2S-SMAD2/3:MTMR4 dephosphorylating it (R-HSA-2187401), DNA-PK autophosphorylation and the
processing steps of NHEJ (R-HSA-5693575, 5687183, 5687360, 5693604). The generator writes both an `input` and a
`catalyst` edge for such a pair (7,898 of 22,397 catalyst pairs on the canonical build), and the AND reads the
complex twice, squaring its fold (specs/012: a 5× complex drives its product 25×).

Inside a cycle this makes the loop gain about 2 at that step, not 1: a positive loop through it does not sit on
the knife-edge (specs/014), it runs away. Of the material-closed positive remnants (specs/106 amendment 2), 11 of
43 on the benchmark and 33 of 56 on the sample contain such a pair. Different loops among them need different
further treatment: a modification pool where the generator's regeneration test drops the path (TGF-beta: the
receptor and ZFYVE9 are consumed as two entities and released as one), a carrier released at the end (Ku in NHEJ).
Each of those first needs the step to read the complex once.

## Rule

`DS_DEDUP_ACTIVATORS=1` (solver; exists since specs/012, default off, no code change): parallel activator edges
from the same node into one reaction are one slot, so a complex that is input and catalyst of its own
reaction contributes its fold once. The rationale is the same as specs/083's (adopted): a species read in two
roles of one reaction is one species. A self-modifying complex converts itself: the step is first order in it.

specs/012 measured this on 2026-09-18, on the curator accuracy alone: held-out −15 (3 fixed / 18 broken,
p 0.0015), tuning −4; it was recorded as refuted. The solver and networks have changed since (about 50
adopted rules and five canonical builds), and the decision rule has changed (specs/109's gate). This is a
re-measurement, and the earlier result is stated here so it cannot be forgotten in the reading.

## Evaluation (pre-registered)

- **Behaviour checks** (specs/107 harness, B2–B5) on the benchmark (`20261024-2211_7f4b28a`) and the sample
  (`20261025-1318_7f4b28a_sample101`). Selection, made once by structure before any solve: per pathway, the
  nodes of cyclic components (the network as the solver sees it, every default) that contain a node which is
  both an `input` and a `catalyst` of a reaction node in the same component; perturbations are roots reaching
  them (at most 8, by stable id), readouts those nodes and up to 10 terminals they reach.
- **Benchmark arm** `dedup110`: `run_arm.sh dedup110 --server DS_DEDUP_ACTIVATORS=1` against `dflt102`, RAF set
  aside for the gate, both figures reported.
- **Gate** G1–G5 as registered in specs/109 (behaviour checks pass on both catalogs; ≥ 90% of un-railed scored
  cases keep direction; held-out losses read by specs/108 against `mild108`; at most 4 unjustified losses and
  McNemar p ≥ 0.05; experimental not significantly worse; raw scores reported). The rule touches reactions
  outside loops too, so every held-out loss is read, not only those in loops.

## Predictions

- D1. G1 holds, and rails among the selected pairs fall on both catalogs (any amount).
- D2. Curator held-out net between −20 and +5; the gate fails on G4 unless most losses are rail-dependent
  (specs/012's 18 broken cases are the prior).
- D3. Experimental moves by fewer than 10 cases either way.
