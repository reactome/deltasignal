# specs/059: a recycled copy reads its root twin

**Status:** registered 2026-10-14, before any code or arm.

## Why

specs/057 R1 and specs/058 R1 lost the same 24 TGF-beta cases (SMAD2/3).
Traced (Neo4j Release 97 and the canonical bundle R-HSA-170834):

- Reactome has one SMAD2 [cytosol] (R-HSA-3006353). SARA recruits it to the
  receptor (R-HSA-170835, the signaling route, through the set SMAD2,3
  [cytosol]); SMURF2 binds it (R-HSA-2176445); MTMR4 regenerates it by
  dephosphorylating p-SMAD2 (R-HSA-2187401).
- The network splits it into copies. The copy SARA takes (7632c373) is fed
  only by MTMR4's output, which comes from it through the receptor: a
  recycled copy. The root copies (in-edges: only the SMAD2:SMURF2 depletion)
  feed SMURF2 binding and MAPK phosphorylation. The only route from a root
  copy to the recycled copy is SMAD2:SMURF2 -| SMAD2: more SMAD2 at the root
  means less signaling SMAD2.
- Canonical scores the cases because no SMAD2 copy is a true root, so the
  specs/024 fallback pins the recycled copy. Any change that makes the root
  copies true roots (specs/057, 058) loses them.

This is the identity-twin split specs/049 met for carrier loops, where
`closed_twin` (adopted, held-out +122) gives a closed free form an in-edge
from its gene's root copy. specs/053 merged such copies in the generator and
lost 493 (a merge welds cycles both ways); this rule adds a one-way edge.

## Census (canonical build)

32 entity copies in 17 pathways are recycled (every producer that is not a
depletion or dissociation edge is downstream of the copy) while a root copy
of the same entity (in-edges, if any, all depletion) does not reach them
along positive (non-depletion) edges. TGF-beta 4 (SMAD2, SMAD3, SMURF2,
R-HSA-173478), R-HSA-3247509 5, R-HSA-1236394 3, R-HSA-453274 3, 13 more
pathways 1–2 each.

## The rule (solver, `DS_RECYCLED_TWIN=on`, default `off`)

Applied after the carrier-loop rule, on the network the solve runs. For each
node with an entity id:
- a **root copy**: every in-edge (if any) is depletion, and it has an
  out-edge;
- a **recycled copy**: not a root copy, not a listed carrier free form, has
  at least one producer that is not depletion or dissociation, every such
  producer is downstream of the copy itself, and no root copy of the same
  entity reaches it along non-depletion edges.

Each recycled copy gains an activating OR in-edge (`recycled_supply`) from
every root copy of its entity: the gene's own level is its supply, read
alongside its regeneration. An observation still wins. The solve reports
`recycled_twin_rule` and `recycled_twin_edges`. Influence scores use the same
network.

## Arms

| arm | build | solver | bench | against |
|---|---|---|---|---|
| control | canonical | defaults (`results/64d23cb`) | default | — |
| `rtwin` | canonical | `--server DS_RECYCLED_TWIN=on` | default | control |
| `rtwin58` | canonical | `--server DS_RECYCLED_TWIN=on` | `--bench DS_PIN_ROOT_IGNORE=depletion` | control, and `deproot` |

Both axes, `--pair-on-valid`, RAF set aside.

## Gate (each arm against the control)

The specs/049 amendment 5 faithfulness gate with RAF set aside: curator
held-out net > 0 at p < 0.01; experimental net loss (RAF set aside) not
significant; no pathway losing more than 10 on either axis; every
experimental loss outside RAF traced. Both readings reported. If an arm
passes, its settings are recommended as defaults (Adam's call).

## Predictions

- `rtwin`: SMAD2/3 TGF-beta unchanged (the fallback already pins the
  recycled copy, which now also reads a root copy that is not pinned).
  Small elsewhere; no sign predicted.
- `rtwin58`: the specs/058 gains (MET OE +10, CBL, PTCH, BRCA1/BARD1, MDM2)
  stay, and SMAD2/3 recover most of their −24, so held-out about +30 against
  the control. Whether p < 0.01 is reached is not predicted.

## Amendment 1 (2026-10-14, after R1, before the arm): add `suppress_only` — POST HOC

R1: `rtwin58` fails only on TGF-beta −12, the SMAD2/3 knockouts, traced to a
depleter built from its own target (SMAD2:SMURF2 -| SMAD2) de-repressing
when the knockout zeroes it. Chosen after seeing that loss, so post hoc.

**Arm `rtwin58s`:** `rtwin58` plus `--server
DS_DEPLETION_OWN_PRODUCT=suppress_only` (existing switch, specs/016; it
applies to every own-product depleter in the catalog, not only SMAD2).
Same gate, against the canonical results and against `rtwin58`.

**Predictions:** the 12 SMAD2/3 KO cases recover (TGF-beta 0 or better).
Elsewhere unknown: specs/016 measured `suppress_only` on an older catalog and
did not adopt it; it loses de-repression where the complex fell because of
its other partner (TP53 KO → less MDM2:TP53 → more free MDM2). If this arm
passes, the decision is Adam's and must weigh that it was chosen after R1.
