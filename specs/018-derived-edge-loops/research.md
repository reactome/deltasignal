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

## Re-evaluation of the welded-era A/Bs on the fixed catalog (pre-registered 2026-09-20 01:50)

Adam: "now that we did that and it is a real fix shouldn't we re-evaluate a
lot of the other things we were looking into?" Yes: every A/B this month whose
gain or loss lived in a loop-heavy pathway was partly measuring which side of
a welded knife-edge a coin landed. Control for all arms below: `ab_fix.tsv`
(`cat_fix`, production solver). Solver knobs are single-catalog (zero churn);
sharing and composition need regenerations (`cat_fix_share`, `cat_fix_comp`)
and are cross-catalog. Code: worktree `560a85c`; curator and experimental
axes for every arm.

| arm | knob | welded-era result | prediction on `cat_fix` |
|---|---|---|---|
| el | elasticity 0.5/0.2/0.95 | held-out +29 / tuning −61 | gain shrinks (the rails it damped were welds): held-out within ±15, tuning cost mostly gone |
| jacobi | `DS_SCC_SWEEP=jacobi` | −0.26pp, churn 14 → 0 | cost shrinks below 0.1pp; still 0 churn |
| dedup | `DS_DEDUP_ACTIVATORS=1` | held-out −15 (specs/012) | the squaring was inside welds: now within ±10, i.e. no longer load-bearing |
| floor0 | `DS_DEPLETION_H_MIN=0` | reverses the +28 floor win | still negative (the EGFR mechanism was acyclic): ≤ −15 |
| hill_log | `DS_AND_MODE=hill_log` | −16 vs hill_sat (specs/010) | still negative; hill_sat stands |
| catbreak | `DS_SCC_BREAK_ROLES=catalyst` | −65 | still negative |
| ownprod | `DS_DEPLETION_OWN_PRODUCT=suppress_only` | −15 / inert | still within ±10 |
| asm | `DS_SCC_BREAK_ROLES=assembly` | +214 on the welded catalog | ≈ 0 now (93 closures remain) — confirms the fix absorbed the rule |
| share | regen `LNG_SHARE_VARIANT_NODES=1` on the fix | −1 held-out / −64 tuning (TP53 basin flip) | tuning cost gone (no weld to flip); held-out within ±10; if so, adopt as correctness |
| comp | regen `LNG_COMPOSITION_EDGES=1` + `DS_COMPOSITION_MODE=limit_novel` | +12 held-out; DSB −62 | DSB loss gone (it was the weld); held-out > +50 if IFN α/β's +100 survives |

Decision rule: each arm keeps or reverses its earlier verdict on held-out
sign with p < 0.05, both axes reported; a change of verdict is recorded as
the fix's effect, not as a new discovery.

### Re-evaluation scored (2026-09-20 04:30; control `ab_fix`, worktree `560a85c`, both axes)

| arm | held-out (fixed/broke, p) | tuning | false change | experimental (cond.) | welded-era verdict | **verdict on `cat_fix`** |
|---|---|---|---|---|---|---|
| el (elasticity 0.5/0.2/0.95) | **+29** (30/1, <1e-4; PTK6 +13, MET +5, EGFR +4) | **−85** (TP53 −94) | −44 | +3 | +29 / −61 | **unchanged** — the gain was never a weld artifact and neither was the TP53 cost: elasticity damps the real MDM2 loop. Not adoptable. Prediction ("gain shrinks, cost gone") wrong on both counts. |
| jacobi | −10 (41/51, n.s.) | −72 (TP53) | +18 | −14 | −0.26pp | **unchanged, worse**: not adoptable; label independence still costs TP53. |
| dedup activators | **−15** (3/18, 0.0015) | +39 (TP53 +51, PIP3 −14) | −14 | 0 | −15 | **unchanged**: the squaring is load-bearing with or without the welds. |
| depletion floor OFF | **−25** (8/33) | **−126** (TP53 −125) | −5 | +3 | floor was +28 | **stands, larger**: the floor is not a loop artifact. |
| hill_log | +2 (n.s.) | −72 (TP53) | +8 | 0 | −16 | **unchanged**: hill_sat stands. |
| catalyst closures | **−104** (5/109; DSB −99) | −33 | −38 | −4 | −65 | **unchanged, worse**. |
| own-product depletion | −5 | 0 | −1 | −2 | inert | **unchanged**: inert. |
| assembly closures | −5 / +8 ALL | +13 | −12 | 0 | **+214** | **absorbed by the fix**, as predicted: 93 closures left do nothing. |
| **sharing** (regen on the fix; 4,094 conditioning drops) | −1 (11/12) | **+18** (TP53 +18) | +11 | +3 | −1 / **−64** | **REVERSED**: the −64 was the TP53 basin flip and it is gone. Sharing is now a free correctness change (nodes −35%, reachability identical) with no accuracy cost on either axis → **adopt** (generator default), pending Adam. |
| **composition** `limit_novel` (regen on the fix, vs same-catalog ctrl) | **−185** (137/322) | −174 | **+426** | −5 | +12 (DSB −62) | **REVERSED the other way**: DSB **−241**, TP53 −151, IFN α/β +100 unchanged. With the welds gone, composition edges are the derived class that closes cycles (cycle-carrying by construction: complex → containing complex → … → complex); they re-weld what the fix unwelded. Refuted as an edge on this network; off, and the LNG flag should stay off. `assembly`-semantics variant: held-out −337. |

Regeneration churn control: `cat_fix_comp` with composition skipped vs `cat_fix`
differs by +12 (one tuning pathway, 6 readouts) — the size of relabelling noise
on the remaining reaction-level cycles.

**Summary.** Of ten welded-era verdicts, eight stand (six knobs remain
not-adoptable or inert, the depletion floor and hill_sat remain right, the
assembly rule is now redundant) and two reverse: **sharing becomes adoptable**
and **composition edges become clearly harmful**. Elasticity is the surprise:
its +29 / −85 is the same trade after the fix as before, which locates it on
the genuine MDM2–TP53 feedback loop rather than on the welds, and says the
remaining loop work is about that loop and the reaction-level recycling
cycles (MHC 675 nodes, Checkpoints 483, HOX 425), not about our edges any
more.

## Review of P7 and the re-evaluation (2026-09-20, two independent reviewers) — verified errata

**E1. The first fix over-reached, and Mitotic G1 −28 is a genuine regression.**
Refusing every *produced* node as a boundary leaf also severed **752 acyclic**
feed-forward assembly links (of 2,749 produced-source assembly edges, 1,997
were cycle-carrying, 752 were not). Traced: RBL1-OE → RBL1:E2F4/5:TFDP ⊣
R-HSA-8964513 → R-HSA-68639 → (assembly, produced source, acyclic) → root
complex R-HSA-68653 → … → readouts; in `cat_fix` 0 of those readouts are
reachable from RBL1. Mitotic G1's largest component is 58 in **both** catalogs
— no weld was involved; the write-up's guess ("read DOWN through the weld")
was wrong. `no_path` grew 1,482 → 1,656, 65 of them former passes (Mitotic G1
25, DSB 11, WNT 8, DNA Damage Bypass 6, EGFR 4, TP53 4). **Corrected fix (LNG
`79feca7`)**: reuse is refused only for a node the root complex can *reach*
(reachability over every edge emitted so far, bridges and depletion edges
included — the output-only definition let 28 welds survive); every acyclic
link is kept. Re-measured as P8 below.

**E2. "Adopted on correctness grounds" was carrying weight the concentration
rule forbids.** DSB is 186 of the +222 conditioned (84%). The pre-registration
wrote the escape hatch in, and the correctness claim was weakened by E1. With
the targeted fix the correctness claim is exact (only cycles Neo4j lacks are
removed); the accuracy claim remains DSB-concentrated and is stated as such.
Conditioned per-pathway nets (the P7 table mixed unconditioned values in):
DSB +186, TP53 +52 (tuning), HDR +20 (t), PIP3 +14 (t), DNA Damage Bypass +11,
WNT +8 (t), HOX +7, TGF-β +6, FGFR1–4 +3 each, Chromatin +3, IL-3/5 +1, RUNX1
−3, EGFR −1, pluripotent 0, Intrinsic Apoptosis 0; Mitotic G1 −12
conditioned.

**E3. The gain is almost entirely "stop saying DOWN".** Of 255 conditioned
held-out fixes, **239 are DOWN → NORMAL with truth NORMAL**; 30 of the 33
broke are correct DOWNs lost to NORMAL. Held-out predicted-DOWN falls 3,088 →
2,726 against 3,357 true DOWNs: DOWN recall fell. Macro-F1 still rose (0.8009
→ 0.8141), so this is not accuracy inflation, but it must be said.

**E4. Conditioning excluded the treated cases.** The 1,550 curator drops (a
gene's node set grew by a fresh leaf) are the fix's *direct* targets: dropped
subset +20 overall (60/40, p 0.057), **tuning dropped −18** (9/27, p 0.004),
held-out dropped +38. Experimental: 222 dropped, **0 fixed / 8 broke
(p 0.008)**, Mitotic G1 −4, WNT −4 — the only significant experimental signal,
negative, and outside the quoted "conditioned 0". Both subsets are now
reported.

**E5. Sharing's "+18 tuning" is 12 relabelling churn + 6.** The churn control
(`cat_fix_comp` with composition skipped vs `cat_fix`: +12, one pathway, 6
readouts, MDM2/MDM4-KO NORMAL→UP) flips the *same* 12 cases sharing flips
(overlap 12 of 12); conditioned, sharing's tuning is **+6** (6/0, p 0.03,
ATM-KO). "The basin flip is gone" is not shown — the MDM2/MDM4 basin still
flips under relabelling, favourably this time. What survives: sharing changed
41 of 24,100 predictions, held-out −1, no accuracy cost on either axis; its
adoption rests on structure (nodes −35%, reachability identical), not on the
+18.

**E6. Pre-registration integrity.** P7's structural clauses (TP53 ≤ 60, DSB ≤
300, components ≤ 215, closures < 100) were written *after* the CSV census
that produced them (LNG commit 3e6d00c at 00:30:10 already states "836 → 56,
1,127 → ~290"; the pre-registration commit is 00:30:37) — they were checks on
the regeneration, not predictions, and are relabelled as such. The accuracy
clauses were set after P6's near-equivalent solver rule had scored +214 /
DSB +190, so they were low-risk. Re-evaluation timing was clean; four of its
ten predictions failed (el both counts, dedup, jacobi, comp) and were scored
as failed.

**E7. Small errors.** MHC's largest component is 680 in `cat_os` (not 675),
675 in `cat_fix`. There are 102 AKT1/AKT2-KO cases: control 100/102, fix
98/102.

**What survived review**: every number reproduces (83.47/0.8009 → 84.81/0.8141;
+222 conditioned, p 1e-43; +260 unconditioned; tuning +81; census exact); at
the stable-id level `cat_fix` = `cat_os` + 475 fresh leaf nodes with **identical
edge multisets in all 92 pathways** (no confound); `cat_fix_comp` minus
composition edges is census-identical to `cat_fix`, so the composition
comparison is zero-churn; "composition re-welds" is confirmed and understated
(with composition edges: 13,823 nodes in cycles vs `cat_os`'s 12,749, DSB
largest 1,525 vs 1,127); no coverage inflation (valid / not-in-network counts
identical).

### Pre-registered P8 — the targeted fix (LNG `79feca7`), committed before the regeneration

Regenerate (`cat_fix2`), production solver, vs the current-tree control on
`cat_os`; both axes; conditioned **and** dropped subsets reported.
- Structure: TP53 ≤ 60, DSB ≤ 300, cyclic components ≤ 215, cycle-carrying
  assembly edges < 100 (same as `cat_fix` — the acyclic links do not affect
  the cycle census); **stable-id edge multiset identical to `cat_os`**; fresh
  leaf nodes fewer than 475.
- Connectivity: `no_path` former-pass losses ≤ 10 (65 in `cat_fix`).
- Accuracy: **Mitotic G1 within ±5** (−28 in `cat_fix`); DSB ≥ +150
  conditioned; held-out ≥ +200 conditioned; tuning ≥ +60; AKT-KO ≥ 96 of 102;
  experimental conditioned ≥ −5 **and** dropped subset not worse than −4.
- Decision: adopt `79feca7` (correctness: only cycles Neo4j lacks removed) if
  the connectivity and Mitotic G1 clauses hold and held-out is not below
  `cat_fix`'s; the accuracy gain is reported DSB-concentrated regardless.
