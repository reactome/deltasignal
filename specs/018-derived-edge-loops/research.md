# Research: where our loops come from — a diff against MP-BioPath's hand-curated networks

**Date**: 2026-09-19 · **Branch**: `feat/018-derived-edge-loops` · **Status**: findings only; no rule built yet

Adam's steer: understand how MP-BioPath's curators resolved loops — not to break
loops ourselves — and derive solver-side rules that stay extensible toward
representing the biology correctly. Everything below is read-only analysis of
two catalogs we already hold: `cat_os` (ours, Reactome v97, 81 scored
pathways) and `mpb_catalog` (MP-BioPath's published networks converted by
`bench/analysis/mpbiopath_network_adapter.py`, 75 pathways). Scripts:
scratchpad `mpb_loop_diff2.py` (edge-for-edge diff) and the census below.

## 1. The two networks join at the reaction level

MP-BioPath's node ids are Reactome dbIds, and Reactome stable-id numbers equal
dbIds (checked on five nodes in Neo4j), so the networks join on the numeric
part. Their network is **bipartite like ours**: 10,918 of their 13,464 node
ids are our entity stable ids and 3,138 are our reaction ids. (A first join
that found zero overlap was a sort-locale bug, not a fact.)

## 2. Edge-for-edge, by role and by cycle membership (75 shared pathways)

Our 36,035 base-id edges; 20,453 lie inside a strongly connected component.
For each class: present in MP-BioPath / both nodes present but edge dropped /
a node absent.

| our edges | kind | present in MPB | nodes present, edge dropped | node absent |
|---|---|---|---|---|
| **cycle** · input | ent→rxn | **48%** | 36% | 16% |
| **cycle** · catalyst | ent→rxn | **64%** | 23% | 14% |
| **cycle** · output | rxn→ent | **49%** | 35% | 16% |
| acyclic · input / catalyst / output | | 35% / 42% / 37% | | |
| **cycle** · assembly | ent→ent | **21%** | 53% | 27% |
| **cycle** · dissociation | ent→ent | **0%** | 81% | 19% |
| **cycle** · depletion | ent→ent | **0%** | 83% | 16% |
| cycle · regulator (neg) | ent→rxn | 6% | 51% | 43% |

The reaction backbone inside our cycles is present in their networks *more*
often than our acyclic backbone is — the curators kept the reactions. What
they essentially never have is our **derived** edge classes: dissociation and
depletion (0%), assembly (21%). Their networks have **26** cyclic components
across the 75 pathways (largest 141, in DNA Repair); ours have 258 (largest
1,127).

## 3. Remove our derived classes and the giants dissolve

Census on `cat_os`, cyclic components over the raw edge lists:

| graph | cyclic comps | nodes in cycles | largest | TP53 | DSB | WNT |
|---|---|---|---|---|---|---|
| all edges | 258 | 12,749 | 1,127 | 1 comp, **836** | 1, **1,127** | 7, 211 |
| − dissociation | 258 | 12,749 | 1,127 | 836 | 1,127 | 211 |
| − dissociation − assembly | 208 | 8,339 | 675 | 1, **56** | 3, 290 | 6, 168 |
| − dissoc − assembly − depletion | 190 | 6,667 | 675 | 2, **34** | 2, **126** | 1, 168 |
| reaction backbone only (input/catalyst/regulator/output) | 190 | 6,667 | 675 | 2, 34 | 2, 126 | 1, 168 |

Dissociation edges point at sink copies and close nothing. **Assembly and
depletion edges are the glue**: only 2,077 assembly edges (5% of all
assembly edges) and 847 depletion edges (25%) sit inside a cyclic component,
yet removing them takes TP53 from 836 nodes to 34 and DSB from 1,127 to 126.
A handful of member → complex and complex ⊣ subunit edges weld separate
reaction-level cycles into one giant component. This is also why every rule
applied uniformly to an SCC failed today (specs/014, 016, 017): the SCC is
not one loop, it is several small reaction cycles welded together by our own
edges.

What remains on the backbone alone is real Reactome cycling — catalytic
recycling and genuine feedback — and it is still large in places: Class I
MHC 675 nodes, Cell Cycle Checkpoints 483, HOX 425, HDR 237. Catalyst edges
are 39% inside cycles (10,238 of 26,465); those are the recycling closures
the loop taxonomy named in June, and MP-BioPath's curators kept 64% of them.

## 4. What this suggests, as a solver rule (not built, not measured)

The derived edge classes are ours: assembly (member → complex, from
hasComponent), depletion (complex ⊣ free subunit, our inference). Their
*feed-forward* meaning is right and load-bearing — a complex cannot exceed its
scarcest member; an abundant complex drains its free subunit. Their
*cycle-closing* instances are the double count: a member that a complex's own
cycle released rebuilding the complex, a complex draining the subunit it is
made from (the own-product case of specs/016, seen again).

Candidate rule, solver-side and edge-class-typed so it does not touch the
network: **a derived edge whose source and target lie in the same SCC is a
recycling closure and does not carry loop signal** — it is read at its
component-entry value (the `supply` mechanism `DS_SCC_BREAK_CATALYST` already
has for catalysts), and SCC detection is re-run without those edges so the
welded component falls apart into its reaction-level cycles, which the
existing damped fixed point then solves separately. Generalising the existing
knob: `DS_SCC_BREAK_ROLES=catalyst,assembly,depletion`.

Pre-registration for that rule belongs in its own spec once Adam has read
this; predictions to state first: TP53's solve becomes converged and
label-independent (836 → 34-node components); the AKT-KO cases stay correct
(the route AKT → MDM2 → TP53 is backbone, not derived); held-out ≥ 0 with
false change *down* in the loop-heavy pathways; and the catalyst-break
history (TP53 +104 / held-out −65 when catalysts alone were frozen) says the
catalyst role must be measured separately from assembly + depletion.

## Results

### Pre-registration (committed before the arms ran)

Arms on `cat_os`, production defaults, from a worktree pinned to the committed
SHA (logged per arm); baseline `ab_fixed_point_ctrl.tsv` — the current tree,
`--max-edges 40000`, 24,100 cases, bit-identical to the September-18 dump on
the 23,908 shared keys. Role lists: **A** `assembly,depletion` (our derived
classes), **C** `catalyst` (Reactome's recycling closures, now with component
recomputation), **ACD** all three. Curator held-out / tuning with
concentration, false change, per-pathway net, McNemar p; experimental axis for
each; the 102 TP53 AKT1/AKT2-KO cases (100 correct in the control); relabel
churn on `cat_perm` for any list that wins.

- **P1 (structure, A).** TP53's largest component after recomputation ≤ 60
  nodes (836), DSB's ≤ 150 (1,127); catalog-wide cyclic components ≤ 210
  (258). The TP53 AKT1-KO solve converges. Closures: ~2,000 assembly, ~850
  depletion.
- **P2 (accuracy, A).** Held-out vs control **≥ +20**, p < 0.05, ≥ 5 pathways
  and ≥ 10 genes, false change **down** in the loop-heavy pathways (the
  welds were the rails); TP53 loses ≤ 10 of the 100 AKT-KO cases (the
  AKT → MDM2 → TP53 route is backbone, not derived); experimental conditioned
  net ≥ −5.
- **P3 (C).** Catalyst closures alone reproduce the June pattern in sign —
  TP53 tuning up, held-out **negative** — because signal-carrying catalysts
  are read at entry; if instead held-out is ≥ 0, recomputation was what the
  June knob lacked.
- **P4 (ACD).** Not better than A on held-out (the catalyst half costs what it
  costs in P3).
- **P5 (decision).** Adopt A only if P2 holds in full; a gain carried by one
  pathway (> 50% of the net) is reported as such and not adopted. If A loses
  held-out, the closure reading (entry value) is the wrong semantics for
  derived edges and the next candidate is reading them feed-forward from the
  recomputed order without `supply` — recorded, not built.
