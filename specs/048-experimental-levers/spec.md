# specs/048: two experimental levers from specs/047 (pre-registration)

**Goal.** Close part of the experimental gap to MP-BioPath (68.79% against
75.74%, canonical build `20261009-2307_132d4b8`, solver `0f0a11e`).
specs/047 showed that the gap comes from how the networks are built, not from
the solver. This spec takes the two narrowest traced classes, re-traced on the
new canonical build, and pre-registers one arm for each. **Experimental is the
primary axis.**

Both arms are scored against the canonical build at solver main, code
defaults, with the `root_cycle` pin protocol.

## Lever A: pathogen-derived nodes held at baseline

**Re-traced on the canonical build.** The specs/047 trace 2 still holds.
- PDPK1 knockout reads 100× at all 7 DOWN readouts in PIP3 (199844, 198335,
  9614997, 200163, 202072, 202074, 198373); truth is DOWN.
- "PDPK1 phosphorylates AKT at T308" (R-HSA-198270) has a negative regulator,
  `N:M:PDPK1` (R-HSA-9755778). That is the SARS-CoV-2 N dimer and M protein
  bound to PDPK1.
- A PDPK1 knockout zeroes the complex and de-represses the reaction, and the
  AKT loop rails.
- In a cell without the virus, the complex does not exist.

**Census** (`~/deltasignal-catalogs/analysis/048/pathogen_census.py`, Neo4j
release 97). The census runs over every entity in the canonical catalog.
- **Rule:** a leaf is pathogen-derived when it has a species and none is Homo
  sapiens; a complex is, when any component is; a set is, when every member
  is. This is the specs/032 drug rule, with species in place of the `Drug`
  class.
- **Count:** 65 of 16,359 entities, 20 of them leaves, in 9 pathways.
- **Pathways:**
  - PIP3 (R-HSA-1257604), VEGF (194138) and CD28 (388841): SARS-CoV-2 N and M.
  - DDX58/IFIH1 (168928): 45 entities, viral RNAs and proteins.
  - IFN α/β (909733), Class I MHC (983169), NK-cell receptors (198933),
    death receptors (140534), and the E-protein complexes in 446728.
- **Experimental reach:** of the ten experimental pathways, only PIP3 is
  touched.

**Design.**
- **Generator:** writes `pathogens.csv` beside `drugs.csv`, by the rule above,
  with variant ids resolved as `drugs.csv` resolves them. The file is always
  written, even when it lists nothing, and the networks are otherwise
  unchanged.
- **Solver:** `DS_PATHOGEN_MODE=propagate` (default) or `inert`. Under `inert`,
  listed nodes are pinned at baseline, as drugs are under `DS_DRUG_MODE=inert`:
  fold 1, never 0, and an explicit observation wins. The solve reports
  `pathogen_rule` and `pathogens_held`. The bench refuses an `inert` arm whose
  build has no `pathogens.csv`.

**Predictions** (arm `path`, against the same build under `propagate`):
- **A-P1.** At least 5 of the 7 PDPK1-knockout DOWN cases become correct.
- **A-P2.** The 3 PDPK1-overexpression UP cases stay correct.
- **A-P3.** Experimental net is +3 to +7, all in PIP3. The other nine
  experimental pathways move by 0.
- **A-P4.** Curator held-out net is within ±15. No direction is predicted: in
  DDX58/IFIH1 the viral RNAs are the pathway's own triggers, so holding them
  at baseline is faithful to an uninfected cell, but that may not be what
  the 2019 curator ground truth assumes.

## Lever B: modification cycles that close through a pooled set

**Re-traced on the canonical build, which replaces specs/047 trace 3.**
- The scaffold loop of trace 3 is gone. Under variant nodes, the scaffold
  released by R-HSA-5672980 is a dead end, and the binding step reads
  `R-HSA-5672717::pool` instead.
- RAF is now 28/49, against MP-BioPath's 45 and 39 on the previous canonical
  build.
- 12 of the 21 failures are one class. ARAF, BRAF, RAF1 or NRAS knockout reads
  p-T,Y MAPK readouts (1268261, 5674340, 5674341) at 36–70×; truth is DOWN.

**Trace: NRAS KO → R-HSA-1268261 = 53.8×**
(`~/deltasignal-catalogs/analysis/048/trace_up.jl`, output in `nras_up.txt`).
The rising path is a curated modification cycle of MAPK:
1. **R-HSA-5672972:** MAPKs (R-HSA-169291) bind the activated RAF complex.
2. **R-HSA-5672978:** RAF phosphorylates MAP2K.
3. **R-HSA-5672973:** MAP2K phosphorylates MAPK.
4. **R-HSA-5672980:** dissociation releases p-T,Y MAPKs (R-HSA-169289).
5. **R-HSA-5675376:** DUSPs dephosphorylate them, giving MAPK monomers and
   dimers (R-HSA-5675361).
6. The released MAPK3 and MAPK1 nodes feed `R-HSA-169291::pool` (over the
   variant cap), which re-enters step 1. The pool node reads 53.8×.

This is MAPK ⇄ p-MAPK, the shape specs/039's balance pools solve. The RAF
build reports 4 pools, all RAS GDP/GTP, and **0 dropped candidates**: MAPK was
never a candidate. There are two reasons, both in `find_pools`:
- **B1. Membership hop.** The R-graph is built from reaction steps only. The
  return from the MAPK3 node into `169291::pool` is a `variant_pool`
  membership edge, not a reaction, so for pool detection the cycle is open,
  while for the solver it is closed. This is a specs/046 seam.
- **B2. Regeneration by entity id.** Amendment 5 keeps a multi-step path only
  if every non-R, non-small input it consumes is output again, compared by
  stable id. Step 1 consumes R-HSA-5672718 (RAS:GTP:activated RAF dimers) and
  R-HSA-5672716 (MAP2K dimers). Step 4 releases R-HSA-5672712 (RAS:GTP:
  activated RAF dimers, a different stable id with the same proteins) and
  R-HSA-5672721 (p-2S MAP2K dimers, the same proteins, modified). Nothing is
  lost, but the ids differ, so the path would be dropped.

**Design** (generator only; the solver is unchanged):
- **B1:** `LNG_POOL_JOIN_MEMBERS=1`. When a step copy's R-containing input is
  a pool node (`variant_pool` or `set_member` target), the step is taken from
  each member node of that pool that contains R. The transition row's
  `source_uuid` is the member. The solver already re-evaluates pool-fed
  inputs of a pool's steps with the pool at baseline (specs/039 amendment
  5), so the pool node is read correctly in the drive.
- **B2:** `LNG_POOL_REGEN=protein`. A consumed partner counts as regenerated
  when every protein it carries is output again by a step of the same path,
  in any form. `entity` (default) keeps amendment 5. Small molecules and
  ubiquitin stay exempt.
- **One arm with both flags.** B1 alone recovers nothing in RAF, because B2
  then drops the path. Other pathways gain pools too; they are counted
  before scoring.

**Structural checks** (must hold before any score is read):
1. The build completes for all 92 pathways.
2. RAF has a pool whose states include a MAPK3 node and a p-MAPK3 node, and
   likewise for MAPK1.
3. Coverage (valid cases, `Pinned:` count) is within 1% of the control on
   both axes.
4. The new pools are counted per pathway (pools, states, multi-step paths)
   and recorded before the score is read.

**Predictions** (arm `mpool`):
- **B-P1.** At least 6 of the 12 RAF knockout→UP cases become correct (DOWN).
- **B-P2.** RAF experimental rises to ≥ 34/49.
- **B-P3.** Experimental net is ≥ +5. Changes outside RAF are not predicted
  in direction: ERBB2, PIP3 and WNT may gain pools.
- **B-P4.** No curator held-out direction is predicted. Pools change wherever
  MAP kinase cycles or set-closed cycles exist (MAPK cascades appear in many
  signalling pathways).

## Decision rule (both arms, separately)

- **Adopt** if:
  - experimental net > 0;
  - curator held-out net ≥ −15;
  - no single pathway loses more than 10 cases on either axis;
  - the structural checks hold.
- Any loss beyond a floor is traced to a case before anything else.
- **Report:** both axes, held-out and tuning, McNemar p, the pathways and
  distinct readouts moved, convergence, and each prediction marked met or
  missed.
- **The arms are independent.** Each is decided against the control on its
  own. If both are adopted, the combination is measured once before it
  becomes canonical.

## Not in scope

These specs/047 classes are left for later specs:
- the PTEN mRNA silo (trace 1);
- the Mitotic G1 E2F release silo (trace 4);
- TP53 DAXX over-coupling;
- THEM4, 7 PIP3 failures, not yet traced.
