# Re-measure: variant-node sharing on the current catalog

`LNG_SHARE_VARIANT_NODES` (introduced in specs/016, LNG #90, off by default)
keys Phase-1 UUIDs on (entity, **Reactome** reaction, role) instead of
(entity, **virtual** reaction, role), so a participant that does not vary
across a reaction's set-member variants stops being photocopied once per
variant.

## Why re-measure

The only measurement (specs/016) ran on `cat_fix`, the **v1** loop-fix catalog
that specs/018 superseded, and the adversarial review then showed its one
visible gain (+18 tuning) was 12 cases of cross-catalog relabelling churn plus
6. The decision to flip the default needs a clean measurement on the current
generator with a same-catalog control.

## What the change does to the graph (measured on HDR, R-HSA-5693567)

BCDX2 (R-HSA-5685316) appears in six (reaction, role) contexts. Four of them
are outputs of reactions that expand into 8 variants each, so BCDX2 was
photocopied 8 times per reaction: **33 nodes → 5**. The other participants of
those reactions are photocopied identically (each of R-HSA-5686410's four
feeders also appears ×8) — the whole reaction is duplicated when only one
participant differs between the copies.

Wiring, before and after:

| | before | after |
|---|---|---|
| BCDX2 as output of R-HSA-5686410 | 8 nodes, each with 1 producer, each feeding 4 dissociation sinks | 1 node with **8 producers, OR-clustered**, feeding 4 sinks |
| BCDX2 as input to R-HSA-5685341 | 1 node | 1 node (unchanged) |

### Semantics, case by case

1. **Variants uniform** (the common case): before, 8 nodes each at v; after, one node at mean(v…v) = v. **Identical.**
2. **Perturbation upstream of all variants**: all 8 at the same value → mean is that value. **Identical.**
3. **Perturbation hits one variant's distinguishing set member**: before, one copy at 0 and seven at 1, and readout resolution takes the **max** over copies, so the effect is **masked**; after, one node at mean = 7/8, so a small effect survives. **This is the one genuine semantic change** — max-over-copies becomes mean-over-producers.
4. **AND multiplication of copies** (the dangerous case): does not arise here — the consumer side was already a single node; the copies only fed dissociation sinks. Confirmed by inspection, not assumed.

## Pre-registration (committed before the regeneration)

Regenerate `cat_share2` from the current generator with
`LNG_SHARE_VARIANT_NODES=1`; control is the current production catalog
`cat_fix2` from the same generator without it. Cross-catalog by necessity
(uuids differ), so **both a same-catalog reference and the churn bound apply**:
specs/019 C3 measured cross-catalog churn at up to 96 cases, all in TP53, zero
held-out.

- **P1 (size)**: nodes fall ≥ 25%; oracle root→terminal reachability
  **unchanged in every pathway** (it was identical in 92/92 on the old
  catalog); cyclic component census unchanged.
- **P2 (accuracy)**: held-out within ±20 and not significant (p > 0.05) — the
  claim is neutrality, and a significant move either way falsifies it.
- **P3 (the case-3 risk)**: false change does not rise by more than 20; any
  rise is concentrated in pathways with high variant multiplicity.
- **P4 (experimental)**: conditioned net ≥ −5.
- **P5 (churn discipline)**: TP53 tuning figures are **not** quoted as an
  effect; only held-out and non-TP53 tuning are.
- **Decision**: flip the default ON if P1 holds and P2–P4 hold — i.e. adopt on
  size and faithfulness grounds given accuracy neutrality. If held-out moves
  significantly in either direction, the change is no longer a pure
  correctness fix and gets argued on its merits.
