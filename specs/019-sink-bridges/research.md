# Research: acyclic sink bridges

## Findings before the rule was built (fixed catalog `cat_fix2`, 2026-09-20)

Error anatomy: 24,100 cases, 3,686 wrong (15.3%). By class: **missed change /
`no_path` 1,526**; false change 1,150; propagator missed 364 + wrong
direction 309; gene or readout not in network 337. Worst pathways: Mitotic G1
520/884 wrong (430 `no_path`), Interferon α/β 354/448 (256 `no_path`).

Oracle on `cat_fix2`: Mitotic G1 957 stable-id-connected root→terminal pairs,
530 uuid-connected, **430 severed** (= its `no_path` count); break kinds
dissociation_sink 264, simple_entity 145, reaction 17. Interferon α/β 642 →
293, 392 severed; break kinds reaction 266, dissociation_sink 78 — the
"reaction" breaks sit downstream of the sink break at R-HSA-909683 →
R-HSA-9710965.

Simulation (bridges added to copies of the bundles, oracle re-run):

| | sinks | acyclic bridges | cycle-closing skipped | fan-out median / p90 / max | uuid-connected | cyclic comps |
|---|---|---|---|---|---|---|
| Interferon α/β | 30 | 13 | 23 | 3 / – / 3 | **293 → 631 of 685** | 2 → 2 |
| Mitotic G1 | 287 | 234 | 152 | 1 / – / 19 | **530 → 705 of 960** | 4 → 4 |
| **catalog (92)** | 35,230 | **10,547** | **56,259** | 1 / 8 / 160 | severed **47.5% → 37.9%** (18,303 → 14,615 of 38,550) | **210 → 210**, largest unchanged |

Five of six candidate bridges would close a cycle; the guard is the rule. The
remaining severance is `simple_entity` and reaction-level breaks, not sinks.

## Pre-registration (committed before the regeneration)

Regenerate `cat_fix2_sb` (LNG `2387c88`, `LNG_SINK_BRIDGES=1`) and
`cat_fix2_sb8` (fan-out cap 8). Arms from a pinned DS worktree (`main`
`72015ce`), production solver, same-catalog control via
`DS_SKIP_EDGE_TYPES=sink_bridge` (zero churn), curator and experimental.

- **P1 (structure)**: cyclic components 210 in every catalog; bridges within
  10% of the simulation (10,547; cap-8 fewer); oracle severance ≈ 37.9%.
- **P2 (the two pathways)**: Interferon α/β ≥ +150 (it has 256 `no_path`
  cases and its connected routes go 293 → 631); Mitotic G1 ≥ +100 (430
  `no_path`, routes 530 → 705).
- **P3 (held-out)**: ≥ +100 vs control, p < 0.05; Interferon α/β will
  dominate it and is stated as such; Mitotic G1 is tuning.
- **P4 (the risk)**: new routes carry perturbations to readouts the curators
  call NORMAL — false change rises, but by less than half the gross gain; the
  cap-8 arm has a smaller false-change rise than the uncapped one and a
  held-out net within 30 of it.
- **P5 (experimental)**: conditioned net ≥ −10 (Mitotic G1 is experimental-
  heavy: this is where the earlier bridges lost).
- **Decision**: adopt the better of the two arms as the generator default
  only if P3, P4 and P5 hold; a P3 gain carried entirely by Interferon α/β
  with Mitotic G1 not moving is a one-pathway fix and is reported as such.
  If held-out < 0, the sink design stands as measured (specs/015) and the
  acyclic guard was not the missing piece.

## Side question (Adam, 2026-09-20): are the depletion edges necessary?

Pre-registered before the arm: `cat_fix2`, production solver, client-side
`DS_SKIP_EDGE_TYPES=depletion` vs the same catalog with them (zero churn), both
axes. Prediction: removing all 4,930 depletion edges is **negative** — held-out
≤ −40 and tuning ≤ −80, concentrated in TP53 (the AKT-KO cases run through
MDM2:TP53 ⊣ TP53; expect ≥ 60 of the 100 lost) and EGFR/GRB2 (specs/011's
case); false change *falls* (fewer routes) while missed change rises more.
If instead held-out ≥ 0, the edges are not load-bearing and the class is a
candidate for removal on fidelity grounds (Reactome never asserted them).

### Depletion ablation scored (2026-09-21; `cat_fix2`, same catalog, zero churn)

| | with depletion (control) | **without (all 3,444 edges skipped)** |
|---|---|---|
| all pathways | 84.71% / mF1 0.8137 | **83.52% / 0.7960** |
| held-out | | **−90** (38 fixed / 128 broke, p < 1e-4), 13 pathways, 43 readouts, 33 genes |
| tuning | | **−196** (TP53 −241) |
| experimental, conditioned | | **−24** (21/45, p 0.0043), 5 of 10 pathways |
| false change | 1,150 | **974 (−176)** |
| AKT1/AKT2 **KO**, TP53 pathway (of 102) | 100 | **2** |
| the same 102 cases, **OE** direction | 84 | **2** |
| predictions changed | | 815; **704** collapse to NORMAL (269 true DOWN, 254 true UP) — of which 489 were correct before |

Per pathway: TP53 −241, MET −60, IFN-γ −34, EGFR −25; only ERBB2 +26 and WNT +9 positive.

**Every clause of the pre-registration held.** Depletion edges are load-bearing
on both axes: removing them does exactly what was predicted — false change
falls (fewer routes) and missed change rises more than twice as far. They are
the only edge class that carries "consumption": an abundant complex drawing
down its free subunit, a phosphatase consuming its substrate. Without them a
knockout upstream of a complex cannot raise the free partner, which is how
98 of the 100 AKT1/AKT2-**knockdown** cases in the TP53 pathway are answered
(AKT ⊣ MDM2 phosphorylation → less MDM2:TP53 → more free TP53). The
over-expression direction of the same 102 cases falls **84 → 2** on the same
ablation. (Label note: these are AKT perturbations scored across all readouts
of `Transcriptional_Regulation_by_TP53`, not a TP53 readout.)

They remain **our inference, not curation** (3,444 edges, **1.21%**; Reactome asserts
no such relation), and the two bounded-ness defects found this month were both
in this class (the floor at zero, specs/011, +28 held-out; the own-product
double count, specs/016, inert). The honest statement for a paper: *a
consumption term is required for knockout propagation through complexes;
Reactome does not record it; we derive it from reaction stoichiometry and
bound it symmetrically.*

## Results (2026-09-21; pinned DS worktree `72015ce`, same-catalog control, zero churn)

**Catalog identity, checked (and a correction).** Every regeneration mints
fresh uuid4s: `cat_fix2`, `cat_fix2_sb` and `cat_fix2_sb8` share **zero**
uuids. So only the **uncapped** contrast is same-catalog and zero-churn
(`sb_live` vs `sb_ctrl`, both on `cat_fix2_sb`, bridges skipped client-side);
the **capped** arm ran on its own regeneration and its comparison against
`sb_ctrl` is **cross-catalog** — the table below originally implied otherwise.
Measured relabelling churn between two regenerations that both lack bridges
(`sb_ctrl` vs `cat_fix2`): **2 cases in one pathway**, far below the ±14–22
seen on `cat_perm`, so the capped arm's figures are not materially confounded
— but its own control (`sb8_ctrl`, same catalog, bridges skipped) is reported
below where available.

| | control | **uncapped** (36,125 bridges) | **cap 8** (4,003 bridges) |
|---|---|---|---|
| all pathways | 84.70% / mF1 0.8136 | **84.85% / 0.8177** | 84.51% / 0.8115 |
| **held-out** | | **−6** (158/164, p 0.78) | **+23** (161/138, p 0.20) |
| tuning | | +42 (65/23) | −68 (51/119) |
| false change | 1,150 | **1,321 (+171)** | 1,300 (+150) |
| experimental, conditioned | | **0** (9/9) | **−15** (8/23, p 0.011) |
| cyclic components / nodes in cycles / largest | 210 / 8,407 / 675 | **210 / 8,407 / 675** | 210 / 8,407 / 675 |

Per pathway (uncapped): Interferon α/β **+100** (100 fixed, 0 broke — all four
gene sets, NORM → DOWN with truth DOWN), Mitotic G1 **+35**, against ROBO
receptors **−43** (44 of its 45 changes are new false changes: SLIT1/ROBO2
perturbations now reach readouts the curators call NORMAL), PDGF −13, and a
diffuse tail. 472 predictions change; 176 are **new** false changes across 17
pathways.

### Scored against the pre-registration

- **P1 (structure) — holds in the part that matters**: the acyclic guard works
  exactly as designed, cyclic components and nodes-in-cycles **identical** to
  the unbridged catalog in both arms. The bridge count is 3.4x the simulation
  (36,125 vs 10,547) because the emitter runs over the in-memory registry
  before final edge dedup, where the simulation used the written bundle.
- **P2 (the two pathways) — fails narrowly**: Interferon α/β +100 (bar: +150),
  Mitotic G1 +35 (bar: +100). The routes are restored — but many restored
  routes deliver a change where the curator recorded none.
- **P3 (held-out ≥ +100, p < 0.05) — fails decisively**: −6 uncapped, +23
  capped, neither significant.
- **P4 (false change < half the gross gain) — fails**: +171 against 223 gross
  fixes uncapped; +150 against 212 capped. The sub-clause held (the two arms'
  held-out nets are within 30).
- **P5 (experimental ≥ −10) — holds uncapped (0), fails capped (−15)**.

**Verdict: not adopted; `LNG_SINK_BRIDGES` stays off.** By the pre-registered
decision rule ("if held-out < 0, the sink design stands as measured and the
acyclic guard was not the missing piece"), this is the **third** measurement
of sink bridging (after 2026-05's −15pp and the silo bridge's −73/−77) and the
first with the cycle objection removed — the objection was not what was wrong
with it. The connectivity gain is real and large (Interferon α/β's severed
routes 293 → 631 in simulation, +100 cases in fact) and it is paid for
one-for-one in false change elsewhere.

The mechanism is the project's oldest finding restated: **on these networks,
added connectivity buys as much over-coupling as it buys reach.** A released
subunit is *one molecule among many copies*; bridging it to every consumer
asserts that the copy the complex released is the copy every downstream
reaction uses. Where that happens to be true (Interferon α/β's ISGF3 branch)
it is worth 100 cases; where it is not (ROBO/SLIT, PDGF, FGFR4, RET) it
manufactures change. Distinguishing the two needs the *identity* of the copy,
which is the node-identity problem (specs/005), not a connectivity heuristic.

**Superseded 2026-09-21 by the flag-expiry policy (specs/020): the emitter is
not merged and PR #92 is closed — the record below is the deterrent, not the
code. The Interferon alpha/beta +100 (0 broken) remains a live unclaimed
finding; it needs copy identity (specs/005), not this emitter.** Originally:
kept as an off-by-default generator flag with the measurement recorded, so the
question does not need re-opening a fourth time.

## Corrections after adversarial review (2026-09-21)

Every headline in this document reproduced independently. Two classes of
problem were found: one that invalidates the sink-bridge *experiment*, and
several errors of fact in the write-up.

### C1 [critical] The sink-bridge arms measured an intervention the design did not describe

`_emit_sink_bridge_edges` decided consumer eligibility with `out_deg[c] > 0`,
reading the same dict it increments. Once a sink received a bridge it became an
eligible "consuming copy" for every later sink of the same entity. On the
shipped catalog **24,967 of 36,125 edges (69%) pointed at another sink** and
consumed nothing (cap-8: 2,053 of 4,003). They concentrate exactly where the
arm lost:

| pathway | bridges | sink→sink | curator net |
|---|---|---|---|
| Interferon α/β | 14 | **1 (7%)** | **+100** (100 fixed, 0 broke) |
| ROBO receptors | 92 | **71 (77%)** | −43 (1 fixed, 44 broke) |
| Mitotic G1 | 3,472 | **3,238 (93%)** | +35 |

The one pathway whose bridges are nearly all genuine gained 100 cases and broke
nothing; the one that lost hardest is 77% artifact. So **the conclusion drawn
below — "added connectivity buys as much over-coupling as reach … it needs the
identity of the copy (specs/005)" — is withdrawn**: it generalises from edges
the design never intended to emit. The cheaper mechanism the review proposes is
also more likely: a sink *is* a readout node, and readout resolution collects
every uuid for a stable id, so a sink→sink edge writes straight into another
readout's value. The genuine-bridge count (11,158) is within 6% of the
simulation's 10,547 — and this also falsifies my explanation of the 3.4x
inflation as "the registry before dedup"; it was this bug.

Fixed (LNG `a81887e`): eligibility against a pre-emitter snapshot. Re-running.

### C2 [high] The emitted bridge set was a random draw per regeneration

Sinks were visited in `sorted()` order over uuid4 labels; the guard is greedy,
so order changes *which* edges survive, not just their order. Constructed case:
swapping two sink labels emits a different edge. At scale, the
cycle-closing-skipped count differs by ~5,400 between the two arms of the same
rule. Fixed: ordering by (stable id, first appearance in the edge list), both
functions of the Reactome data. Tests added; the three guarantees (no sink is
a bridge target; a bridge can make a later candidate cycle-closing; visit
order is data-determined) each now fail their mutation — all three passed
before.

### C3 [high, methodology] Cross-catalog churn is ~96 cases, all in TP53

Two regenerations that differ only in uuids and in bridges that were then
skipped client-side (`sb_ctrl` vs `sb8_ctrl`) differ in **96 predictions, all
in Transcriptional_Regulation_by_TP53, zero held-out**. An earlier pair of the
same kind gave 2 — so the earlier "2-case churn" was a lucky draw, not a
bound. This is the MDM2 loop's basin flipping under relabelling (specs/013,
018). Consequences to apply going forward:
- **Held-out figures across catalogs are safe** (0 churn in both pairs).
- **Tuning figures across catalogs are not**: the +109/+81 tuning numbers in
  specs/018's P7/P8 are within this noise band and should not be quoted as
  effects. The held-out +173/+222 there stand.
- Cross-catalog arms must either report a same-catalog control (as the capped
  arm now does) or omit TP53-dominated tuning claims.

### C4 [medium] Errors of fact, corrected above

- Depletion edges: **3,444**, not 4,930; **1.21%** of edges, not 1.5%.
- "815 changed; 489 collapse to NORMAL" conflated two quantities: **704**
  collapse to NORMAL, of which 489 had been correct.
- "AKT-KO → TP53" is AKT1/AKT2 perturbations scored over all readouts of the
  TP53 pathway, not a TP53 readout; the over-expression direction (84 → 2) was
  omitted and is now reported.
- The experimental axis has **no held-out split at all** — all 849 cases lie in
  the ten tuning pathways. "P5 (experimental) holds" is a tuning-set result.
- `sb8_ctrl` did exist by the time it was cited: capped arm vs its own control
  is held-out **+23** (161/138, p 0.20), tuning +24, false change +150.
- "Third measurement of sink bridging" overstates continuity. 2026-05's −15pp
  was an ablation of the whole dissociation class under a shared-node model —
  it is what *created* sinks; the silo bridge was solver-side, to in-degree-0
  copies, guarded by a reach cap rather than a cycle test. Fair phrasing:
  *third intervention in the family; the first targeting dissociation sinks,
  bridging to consuming copies, with a cycle guard, pre-registered.*
- Logging: `fan` recorded only emitted sinks (cap-8 hid 1,290 from its own
  statistic) and `skipped` was counted before the cap check — both fixed; a
  negative `LNG_SINK_BRIDGE_MAX_FANOUT` is now an error rather than a silent
  total cap.

### Verdict status

**The sink-bridge result is withdrawn pending a re-run** with the corrected
emitter (~11.2k genuine bridges). The depletion ablation is unaffected: the
review confirmed the client-side skip removes exactly the 3,444 depletion
edges with zero collateral, conditioning drops zero cases, and every figure
reproduces — depletion edges remain load-bearing.

## Re-run with the corrected emitter (2026-09-21; LNG `a81887e`, catalogs `cat_fix3_sb` / `cat_fix3_sb8`)

Structure: **11,069** bridges uncapped (simulation predicted 10,547 — within
5%), 3,416 at fan-out ≤ 8, **0 sink→sink**, cycle census **210 / 8,407 / 675**
identical to the unbridged catalog. Each arm has its own same-catalog control
(per C3), so no cross-catalog churn.

| | control | uncapped | cap 8 |
|---|---|---|---|
| all pathways | 84.71% / 0.8137 | 84.80% / 0.8173 | 84.87% / 0.8179 |
| **held-out** | | **−11** (169/180, p 0.59) | **−4** (155/159, p 0.87) |
| tuning | | +35 (68/33) | +38 (67/29) |
| false change | 1,150 | **+195** | +173 |
| experimental, conditioned | | **−2** (9/11) | **−2** (7/9) |

Per pathway: Interferon α/β **+100** (100 fixed, 0 broke), Mitotic G1 +31,
ROBO **−44** (1 fixed, 45 broke), RET −18, PDGF −17, IL-3/5 −13; 200 new
false changes in 18 pathways.

### Scored against the original pre-registration

P1 structure **holds** (bridge count matches the simulation, cycle census
unchanged). P2 **fails** (IFN +100 vs bar +150; Mitotic G1 +31 vs +100).
P3 **fails** (held-out −11 / −4, neither significant, bar +100). P4 **fails**
(false change +195 against 237 gross fixes; the sub-clause holds — the arms'
held-out nets differ by 7). P5 **holds** (experimental −2, bar −10).

**Verdict: not adopted, `LNG_SINK_BRIDGES` stays off** — the same verdict the
buggy run reached, now on a valid measurement.

### What the bug did and did not change

| | buggy (69% sink→sink) | corrected (0) |
|---|---|---|
| held-out uncapped | −6 | −11 |
| false change | +171 | +195 |
| Interferon α/β | +100 | +100 |
| ROBO | −43 | **−44** |
| Mitotic G1 | +35 | +31 |
| PDGF / RET | −13 / −10 | −17 / −18 |

**The review's inference was wrong even though the defect was real.** It
reasoned that ROBO lost because 77% of its bridges were artifacts and
Interferon α/β gained because only 7% were — but with *zero* artifacts ROBO
still loses 44 and Interferon α/β still gains exactly 100. The artifact edges
were nearly inert for a structural reason neither of us stated at the time: a
sink→sink edge points at a node with no outgoing edges, so it can only change
that sink's own value, and only where the sink is itself a readout. The defect
had to be fixed — 69% of the emitted edges were not the intervention — but it
was not what produced the result.

**So the conclusion withdrawn in C1 is reinstated, now on evidence that
supports it**: reconnecting a released subunit to the copies that consume it
buys reach and over-coupling in roughly equal measure. Interferon α/β is the
case where the released copy really is the copy the downstream reactions use
(+100, nothing broken); ROBO, RET, PDGF and IL-3/5 are cases where it is not,
and the bridge manufactures change. Which of the two a bridge will be is a
question about the identity of the copy — specs/005 — and cannot be decided by
a connectivity rule, a cycle guard, or a fan-out cap (the cap changed the
held-out net by 7 cases).
