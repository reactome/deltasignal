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

### Scoring, first three arms (2026-09-20 00:20; worktree `560a85c`, control = current-tree `fixed_point`, 24,100 cases)

| role list | headline | held-out (fixed/broke, p) | tuning | false change (loop-heavy / all) | experimental (conditioned) | AKT-KO TP53 (of 100) |
|---|---|---|---|---|---|---|
| **A** `assembly,depletion` | **84.59% / mF1 0.8107** (control 83.47 / 0.8009) | **+226** (251/25, p < 1e-4), 6 pathways, 42 readouts, 45 genes, dominant perturbation 5% | +44 | 637 → **278** / 1,437 → **1,006** | +1 (17/16) | **4** |
| **C** `catalyst` | 83.16 / 0.7967 | −65 (12/77) | −11 | 637 → 571 / −77 | −1 | 12 |
| **ACD** all three | 84.00 / 0.8034 | +100 (241/141) | +26 | 637 → 276 / −435 | 0 | 4 |

TP53 probe under A: 1 component (836) → 2, largest **23**, from 2 assembly + 5
depletion closures.

**P1 — holds.** TP53 ≤ 60 (23); DSB's giant dissolves (the +190 below is that
component); components recomputed; closures counted.

**P2 — holds on every clause but one.** Held-out +226 with p < 1e-4, 6
pathways, 45 genes, false change down 359 in the loop-heavy pathways (the
welds were the rails, as predicted), experimental +1. **The AKT clause fails
badly**: 92 of the 100 correct AKT1/AKT2-KO cases go UP → NORM. The route
AKT → MDM2 → TP53 is backbone, but the *step that raises TP53* is the depletion
edge MDM2:TP53 ⊣ TP53 — and that edge closes a cycle (TP53 → MDM2 transcription
→ MDM2 → MDM2:TP53), so under `depletion` it is read at its entry value and the
de-repression never happens. Depletion closures carry genuine regulation, not
only our double count. TP53 nets −4 overall because DAXX / ELL / BRCA1 cases
(181) are fixed while the AKT cases (92) and others break.

**P3 — holds**: `catalyst` alone is held-out −65, DSB −48 — the June pattern
with recomputation added; recomputation was not what the June knob lacked.

**P4 — holds**: ACD (+100) is worse than A (+226); the catalyst half costs
what P3 says.

**P5 — the concentration clause fires.** DSB Repair is **+190 of the +226**
(84%); the other five held-out movers are DNA Damage Bypass +15, EGFR +8, HOX
+7, TGF-β +6, Hedgehog 0. By the rule written before the arm, a gain carried
by one pathway is reported as such and not adopted. The mechanism explains
the concentration — welded giant components exist in two pathways, DSB
(held-out) and TP53 (tuning), so that is where a weld-breaking rule can act —
but the rule is the rule. What DSB's +190 is: **198 cases DOWN → NORM with
truth NORM** — the welded 1,127-node component was collapsing to zero under
knockouts of PARP1, PARP2, FEN1, POLQ, RTEL1, XRCC5/6, MUS81 (13–15 cases
each); with the weld broken those knockouts stay local. That is the
over-coupling error class the project identified as its largest, removed in
the pathway where it was largest.

### Pre-registered P6: `assembly` alone (committed before the arm)

The AKT failure says depletion closures are load-bearing regulation; the
census says assembly alone removes most of the weld (TP53 836 → 56; DSB
1,127 → ~290). Arm **As** = `DS_SCC_BREAK_ROLES=assembly`, same control, plus
the experimental axis, plus relabel churn on `cat_perm` for **A** and **As**.
Predictions: As keeps ≥ 90 of the 100 AKT-KO cases (the depletion edge feeds
back again); DSB keeps ≥ 120 of its +190 (its weld is mostly assembly); held-out
≥ +100 with the same distribution caveat; false change still down ≥ 250 in the
loop-heavy pathways; experimental ≥ −5; relabel churn 0 for both A and As on
the recomputed components. Decision: if As holds all of those, it is the
candidate — still DSB-concentrated by mechanism, to be stated as such; the
depletion half is then a separate question (which depletion closures are our
double count and which are regulation).

### The generator-side cause, and pre-registered P7 (committed before the regeneration)

Adam: "if these edges are the breaking apart of complexes in terminal outputs
and root inputs, they should not create loops as they would get different
uuids." They should, and they do not: boundary expansion's `_leaf_uuid`
reuses a leaf's *existing* node whenever the protein already has one anywhere
in the network, including a copy that a reaction **produces** downstream of
the root complex. Census on `cat_os`: of 2,077 cycle-carrying assembly edges,
**1,994 have a produced source and a root-complex target**; 40,208 acyclic
assembly edges are proper unproduced leaves into root complexes. The depletion
closures are a different thing — MDM2:TP53 ⊣ TP53 is genuine regulation inside
a genuine feedback loop, which is why silencing them broke the AKT cases.

Fix in the generator (LNG branch `fix/boundary-leaf-no-produced-reuse`, 3e6d00c):
a boundary leaf reuses only an *unproduced* node; a produced copy gets a fresh
leaf uuid. `LNG_BOUNDARY_LEAF_REUSE=any` restores the old behaviour.

**P7 — regenerate the full catalog with the fix (`cat_fix`), solve with the
production solver (no `DS_SCC_BREAK_ROLES`), compare to the current-tree
control on `cat_os`.** Cross-catalog (uuids differ), so churn is part of what
is measured; conditioning drops are expected only where a gene's node set
changed.
- Structure: TP53's largest component ≤ 60 nodes, DSB's ≤ 300; catalog cyclic
  components ≤ 215; cycle-carrying assembly edges < 100 catalog-wide.
- Accuracy: held-out ≥ +150 (the assembly half of A's +226), DSB ≥ +120,
  false change down ≥ 250 in the loop-heavy pathways, the 100 AKT-KO TP53
  cases stay ≥ 90 (depletion untouched), experimental conditioned ≥ −5.
- Interpretation: the same DSB concentration as A, for the same reason; the
  fix is adopted on *correctness* grounds (the network no longer contains
  cycles Neo4j does not have) if the accuracy clauses hold, and reported as
  a DSB-concentrated accuracy gain.

### P6 scored — `assembly` alone (2026-09-20 00:55; worktree `560a85c`)

| | headline | held-out (fixed/broke, p) | tuning | false change (loop-heavy / all) | experimental | AKT-KO TP53 (of 100) | relabel churn |
|---|---|---|---|---|---|---|---|
| **As** `assembly` | **84.74% / mF1 0.8136** (control 83.47 / 0.8009) | **+214** (241/27, p < 1e-4), 5 pathways, 40 readouts | **+92** (110/18) | 637 → **374** / 1,437 → **1,127** | 612 → 612 (0) | **100** | 15 of 24,100 (fixed point 14) |
| A `assembly,depletion` (for comparison) | 84.59 / 0.8107 | +226 | +44 | 278 / 1,006 | +1 | 4 | 17 |

Per pathway (As): DSB **+190**, TP53 **+46**, HDR +20, DNA Damage Bypass +15,
PIP3 +14, ERBB2 +12, HOX +7, TGF-β +6, EGFR −4, Checkpoints −2.

- AKT clause **holds** (100 of 100 kept; the depletion edge MDM2:TP53 ⊣ TP53
  feeds back again). DSB ≥ 120 **holds** (+190). Held-out ≥ +100 **holds**
  (+214). False change ≥ −250 in loop-heavy **holds** (−263). Experimental
  ≥ −5 **holds** (0). **Relabel churn 0 fails** (15 vs 14): the reaction-level
  components that remain after unwelding are still iterated by Gauss-Seidel,
  whose result depends on visit order (specs/013); unwelding does not touch
  that, `DS_SCC_SWEEP=jacobi` does.
- Concentration, stated: DSB is 190 of the held-out +214 (89%). The tuning
  +92 is spread over six pathways with TP53 +46. The rule acts where a derived
  edge welded a giant, and there are two such pathways.

`assembly` alone is strictly better than `assembly,depletion` on every axis
except the loop-heavy false-change count, and it does not touch the AKT
route. It is the candidate solver rule — pending P7, which asks whether the
same result is obtained by not constructing the welding edges in the first
place, in which case no solver rule is needed for assembly at all.

### P7 scored — the generator fix, regenerated catalog, production solver, no rule (2026-09-20 01:30)

Regeneration: LNG `3e6d00c` clean, 93 pathways, 737 s. Census of `cat_fix`
vs `cat_os`: cyclic components **207** (258), nodes in cycles 8,449 (12,749),
cycle-carrying assembly edges **93** (2,077); largest components TP53 **56**
(836), DSB **290** (1,127), WNT 168 (211), MHC 675 (675 — a reaction-backbone
component, untouched as it should be). Every structural clause holds.

| | control (`cat_os`, current tree) | **`cat_fix`, production solver** |
|---|---|---|
| all pathways | 83.47% / mF1 0.8009 | **84.81% / 0.8141** |
| held-out, conditioned (18,358) | 86.15% | **87.36%** — **+222** (255/33, p < 1e-4) |
| held-out, unconditioned (19,000) | | +260 (306/46), **15 of 71 pathways moved**, 56 readouts |
| tuning, conditioned | 77.60% | 79.53% — **+81** (120/39) |
| false change (23,908 shared keys) | 1,437 | **1,043** (loop-heavy 637 → 322) |
| experimental, conditioned (627) | | **0** (4/4); unconditioned −8 (4/12, p 0.08; WNT −4, Mitotic G1 −4) |
| AKT-KO TP53 (of 100) | 100 | **98** |
| conditioning drops | | 1,550 (fresh boundary leaves carry the gene, so gene node sets grew) |

Per pathway: DSB **+210**, TP53 +49, DNA Damage Bypass +22, HDR +20, PIP3
+14, WNT +9, pluripotent stem cells +9, HOX +7, TGF-β +6, FGFR3 +3;
**Mitotic G1 −28**, Intrinsic Apoptosis −6, RUNX1 −3.

- Held-out ≥ +150 **holds** (+222 conditioned). DSB ≥ +120 **holds** (+210).
  False change ≥ −250 in loop-heavy **holds** (−315). AKT ≥ 90 **holds** (98).
  Experimental ≥ −5 conditioned **holds** (0).
- Concentration, stated: DSB is 210 of the held-out +260 (81%); but fifteen
  held-out pathways moved, three times as many as under the solver rule, and
  twelve of them moved up.
- The one new loss: **Mitotic G1 −28** — RBL1-OE (17) and RBL2-OE (15) cases
  that read DOWN correctly through the welded component now read NORMAL. A
  correct answer that depended on a wrong edge; to be traced before it is
  called a regression rather than a revealed gap.

**Decision (per the pre-registered rule): the generator fix is adopted on
correctness grounds** — the network no longer contains 1,900 cycles that
Reactome does not have — and its accuracy effect is reported as a
DSB-concentrated +222 held-out with broad small gains. It makes the solver's
`assembly` closure rule redundant (93 closures remain catalog-wide);
`DS_SCC_BREAK_ROLES` stays available, off by default, as the measured record.
Everything measured this month on the welded catalog is now re-evaluated
against `cat_fix` (next section).
