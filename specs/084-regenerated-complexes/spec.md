# specs/084: a complex the pathway only regenerates follows its components' starting copies

**Status:** pre-registered 2026-10-21, before the implementation and the arm. Arm `cct`: **gate fails** (held-out −8, all PDGFRB knockout through specialization-bridged copies); ATM fixed as predicted (+6 curator, +4 experimental). Not adopted. Candidate amendment in `research.md`.

Adam, 2026-10-21: "ok" to the regenerated-complex rule as specs/084.

## The defect (specs/080 correction, traced)

In Cell Cycle Checkpoints (R-HSA-69620), the DSB signalling complex that carries ATM
(R-HSA-5683605, "DNA DSBs:p-MRN:p-S1981,Ac-K3016-ATM:KAT5:…:BRCA1-A complex"):
- recruits CHEK2 (R-HSA-5683735);
- phosphorylates it ("Phosphorylation and activation of CHEK2 by ATM", R-HSA-69891);
- is released unchanged.

Nothing in the pathway makes the complex, which comes from DNA Double-Strand Break
Response. It is a closed carrier loop (specs/049), and `closed_twin` feeds a closed free form
only from a starting copy of the same entity, of which this pathway has none. So the complex
reads baseline, and ATM's starting copy (p-S1981,Ac-K3016-ATM, R-HSA-5693527), which is one
of its components, never reaches it. ATM knockout and overexpression read NORMAL at
p-CHEK2 and p-CDC25C; the curators say DOWN and UP.

## Structural count, before the implementation (`count_component_twins.py`, canonical `20261020-2038_23ef3e1`)

- Closed carrier free forms with no starting copy of their own entity: **193**.
- **49** of them (16 distinct entities, 7 pathways) contain entities that do have a starting
  copy in the pathway. Keeping only the maximal ones (an entity inside another contained
  entity is skipped) gives **202** component-root edges, exactly one starting copy per
  contained entity.
- The pathways: DSB/Cell Cycle Checkpoints (ATM complex), Base Excision Repair (R-HSA-73942,
  8 entities), Signaling by NOTCH (R-HSA-453274, 3), MHC class I antigen presentation, and
  R-HSA-186797, R-HSA-373760 and R-HSA-5693567.
- The other 144 have no starting copy of any kind and are untouched.

## The rule (`DS_CARRIER_COMPONENT_TWIN`, default `off`; `on`)

Under `DS_CARRIER_LOOP_MODE=closed_twin`, when a closed carrier free form has no starting
copy of its own entity:
- it gains an activating AND in-edge (type `carrier_component`) from the starting copy of
  each maximal contained entity that has one (containment table; cofactors excluded; never
  from a listed free form);
- so it follows the product of its components' starting copies, as the boundary expansion
  builds a starting complex from its components;
- one way only: nothing feeds back into the components.

The solve reports `carrier_component_rule` and `carrier_component_edges`. The benchmark's
root detection is unchanged: the edges exist only in the solver, as `closed_twin`'s do.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `cct` | canonical `20261020-2038_23ef3e1` | `--server DS_CARRIER_COMPONENT_TWIN=on` | `results/8190d83/dflt` |

`dflt` is scored-prediction identical to today's main: specs/083 changed no class. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade`, on both axes, case by case.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10 (RAF set aside).

## Predictions

1. Cell Cycle Checkpoints: ATM knockout and overexpression at p-CHEK2 (R-HSA-5683784) and
   p-CDC25C (R-HSA-156496) read DOWN and UP, at least 6 of these 8 cases fixed.
2. Fewer than 2 pathways lose more than 5 curator cases.
3. No prediction for the held-out net. 7 pathways move at most, and Cell Cycle Checkpoints
   is a tuning pathway.
