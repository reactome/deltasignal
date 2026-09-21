# A curator-faithful oracle: read the pathway from Reactome and check what the solver can reach

**Created**: 2026-09-19
**Status**: METHOD BUILT. Two pathways: right at the entity level, ~half-severed at the node level. A simulated 2-hop composition bridge (fan-out 1–2) restores 392→43 severed routes in Interferon and 462→165 in Mitotic G1 — the non-broadcasting repair. IMPLEMENTED and validated on Interferon (392→52); catalog-wide it recovers 11% of severed curator routes, the protein→complex variant 8% more with hub fan-out to 840 (not built). ~33% of curator routes stay cut: instance multiplicity + released-subunit sinks. A/B FINAL (81 pathways): comp_live −85 held-out (false change −219, Interferon +74), comp_group −322. No aggregation over node copies is right; the fix is upstream (LNG_SHARE_VARIANT_NODES, implemented, validating on HDR). Regulator-dropping hypothesis REFUTED (all 61 set-typed regulators carried).
**Tool**: `scratchpad/reactome_direct.py` (to be promoted into `bench/analysis/`)


> **Reproduction note (2026-09-21).** The `LNG_SHARE_VARIANT_NODES` and
> `LNG_BOUNDARY_LEAF_REUSE` settings in the commands below **no longer exist**:
> variant sharing and downstream-free leaf reuse are unconditional, and setting
> either name is now a hard error (LNG #93, #94). That is deliberate — a stale
> value must not let a run quietly measure the default. It does mean these
> commands will abort rather than reproduce; the recorded numbers are the
> evidence. To re-measure the *old* behaviour you have to check out the
> generator commit named in the arm, not set a flag.

## 1. Why this exists

Adam, on the claim that the benchmark's resolution was the limit: *"I don't
think that is true. What the curators did was look at what they thought would
happen with the pathway... If you traversed all the tests from root input to
terminal output in Neo4j you could see what you think should happen and see if
our modelling is doing that."*

He is right. Comparing the model to 2019 curator calls conflates our defects
with theirs and with version skew. Comparing the model to a **direct reading of
the same Reactome release** isolates ours: wherever Reactome says a route
exists and the graph the solver runs on says it does not, that is a
construction defect regardless of what the curator wrote. It is the method that
found the dropped BlackBoxEvents, made systematic.

## 2. The oracle

For one pathway, straight from Neo4j (no LNG code involved):

1. Every `ReactionLikeEvent` under the pathway, with `input`, `output`,
   `catalystActivity→physicalEntity`, and `regulatedBy→regulator` (signed:
   `NegativeRegulation` = −1).
2. **Composition, the way a curator reads it:** `hasComponent` in the
   assembly direction only (component → containing complex), and set
   membership both ways (a set *is* its members). The dissociation direction
   (complex → component) is **excluded**: it is the broadcast route the
   −15pp A/B condemned, and including it inflates reachable pairs ~7×.
3. Roots = entities with no incoming reaction edge, restricted to protein-
   bearing (`R-HSA-`) and not in `cofactors.csv`; terminals = no outgoing.
4. Signed reachability root → terminal.

Then, against LNG's network for the same pathway, at two resolutions:

- **stable-id level** (sets resolved to their leaves): does a route exist?
- **uuid level** — the graph the solver actually runs on: does a route exist
  between any uuid of the root and any uuid of the terminal?

## 3. Two things the oracle ruled out on the way

- **Regulations on `CatalystActivity`.** If they existed and LNG read only
  `ReactionLikeEvent.regulatedBy`, every one would be silently dropped.
  Catalog-wide: 16,344 human reactions, 3,003 regulations on reactions,
  **0 on catalyst activities.** Not a bug.
- **"Missing" `DefinedSet`s.** A reaction-only first pass reported 9 missing
  roots in Interferon; all were sets, which LNG expands to members by design.
  Representation, not omission. The oracle now resolves sets before judging.

## 4. Results, strict oracle (recomputed 2026-09-19 with small molecules and cofactors excluded as carriers; Interferon unchanged, Mitotic G1 revised slightly)

| | Interferon α/β (R-HSA-909733) | Mitotic G1 (R-HSA-453279) |
|---|---|---|
| reactions / entities | 25 / 83 | 101 / 181 |
| curator-reachable root→terminal pairs | 719 (583 +, 136 −) | 962 (435 +, 527 −) |
| root entity absent from LNG | 34 pairs, **1 negative regulator** | 2 pairs, 2 catalysts |
| both present, unreachable at stable-id level | 41 | 3 |
| **LNG agrees, stable-id level** | **642 (89%)** | **957 (99%)** |
| **LNG agrees, uuid level (solver's graph)** | **293 (43%)** | **537 (56%)** |
| **severed by shard splits** | **392** | **423** |

**The network is right at the entity level and roughly half-severed at the
node level.** The stable-id agreement is partly an artifact — a dissociation
sink shares its stable id with the live node it is disconnected from, so a
route "exists" at that resolution while the solver cannot use it. The uuid
row is the true one. This is the mechanism `specs/015` traced case by case,
now measured as a fraction of everything a curator can reason to.

It reframes the silo history: the earlier bridges were not wrong that routes
were missing — they were wrong about *which node* to reconnect. Reconnecting
through leaf subunits (fan-out median 112) broadcasts; a route through the
containing complex (a complex sits inside few complexes) need not.

## 5. Simulation of the non-broadcasting repair

At uuid level, add an edge from every non-sink node of `X` to every node of
`Y` for each direct `Y hasComponent X`, then recount. No regeneration, no
solver change — a pure reachability test of whether the hierarchy is the gap.

| | Interferon α/β | Mitotic G1 |
|---|---|---|
| hierarchy edges added (uuid) | 99 over 23 relations | 167 over 112 relations |
| containing-complex fan-out | median 1, max 1 | median 1, max 2 |
| connected before → after | 293 → 300 | **539 → 836** |
| **severed before → after** | 392 → 385 | **462 → 165** |

**Mitotic G1: the hierarchy IS the gap** — 259 of 423 severed routes (61%)
come back through edges with fan-out 1–2, against the leaf bridge's 112. This
is the repair the silo record asked for and never got: it reconnects what a
curator reconnects and cannot hub-flood. It is a simulation; a real arm needs
LNG to emit `composition` edges along direct `hasComponent` (specs/015 §6) and
a held-out A/B on both axes with the concentration columns.

### Two hops

Allowing the containing complex to be up to **two** `hasComponent` steps away
(ISGF3 → ISGF3:KPNA1 → ISGF3:KPNA1:KPNB1, where the middle complex has no
node):

| | Interferon α/β | Mitotic G1 |
|---|---|---|
| hierarchy edges added (uuid) | 121 over 41 relations | 302 over 169 |
| containing-complex fan-out | median 1, **max 1** | median 1, max 2 |
| connected before → after | 293 → **642** | 537 → **796** |
| **severed before → after** | **392 → 43** | **423 → 164** |
| still severed, by break kind | 40 no stable-id path, 1 reaction, 1 entity | **131 released-subunit sinks**, 26 entity, 4 reaction, 3 no path |

**Interferon: essentially every severed curator route is restored** by 121
edges of fan-out 1. Mitotic G1's remainder is the released-subunit class
(section 4 of specs/015) — the broadcast route — and is correctly *not*
touched by this repair.

This is the repair with no prior against it: every bridge that lost was a
leaf-subunit broadcast at fan-out ~112; this is complex → containing complex
at fan-out 1–2.

### Implemented and validated on the real network

`LNG_COMPOSITION_EDGES=1` (default off): `_emit_composition_edges` in
`logic_network_generator.py`, `get_containing_complexes` (upward, ≤2 hops) in
`neo4j_connector.py`; 8 unit tests with Neo4j stubbed
(`tests/test_composition_edges.py`). Solver side: `composition` shares
`assembly`'s limiting-reactant semantics (`reaction_model.jl`, one token;
pinned in `test_propagator_invariants.jl`, 74 assertions).

Interferon α/β regenerated with the flag on:

| | before | after (real network) | simulated |
|---|---|---|---|
| composition edges | 0 | **15** (0.09% of 16,091 edges), fan-out median 2, max 4 | 121 (uuid-level) |
| uuid-level connected | 293 / 685 | **633 / 685** | 642 |
| severed | 392 | **52** | 43 |

The remaining 52 are 40 with no stable-id route at all plus a handful of
one-off breaks. **The connectivity ladder's end-to-end tier flipped**:
`test_cytosolic_isgf3_reaches_the_readout` went xfail → XPASS(strict) — the
exact signal that test was written to raise. Tiers 2–3 (containment export of
intermediates; assembly/dissociation sharing a uuid) still xfail: they encode
the leaf-bridge repair, which is the one measured to broadcast and is
deliberately not made.

### Mitotic G1 on the real network: the simulation was mis-specified

| | simulation | real regenerated bundle |
|---|---|---|
| composition edges | 302 | 179 |
| severed | 423 → **164** | 423 → **420** |

Three routes recovered where the simulation promised 259. Cause: the
simulation's `fetch_direct_containers` took *every* reaction participant as a
source, so it silently included **protein → containing-complex** edges (live
CDK4 node → every complex containing CDK4); the implementation, per the design
stated in §5, emits **complex → complex only**. Interferon's gap is
complex → complex (ISGF3 → ISGF3:KPNA1:KPNB1), so there the real bundle matched
the simulation (392 → 52 vs 43). Mitotic G1's gap is protein → complex, so it
did not. The Mitotic G1 simulation number was reported above as validating the
design; it validated a different, broader edge class. **The 11% catalog-scale
recovery is the honest figure for what was built.**

Confirmed to the route by re-simulating with the source class split:

| Mitotic G1 simulation sources | edges | fan-out | severed |
|---|---|---|---|
| complex → complex (what LNG emits) | 86 | 1 / 1 | **420** — identical to the real bundle |
| protein → containing complex | 209 | 1 / 2 | **167** |
| both | 302 | 1 / 2 | 164 |

The broader class — a live protein node → the complexes that contain it, ≤2
hops — is the *assembly* direction, which the May A/B measured **positive**
(+4.2pp) when emitted at root complexes only. Extending it to all complexes is
a plausible second fix and its hub-protein fan-out is the risk to size first.
The oracle now simulates the two source classes separately.

### Sizing the protein → complex variant catalog-wide — narrow, and a hub risk

Simulated on the 24 regenerated bundles *on top of* the real complex → complex
edges:

| | routes severed | recovered | edges | worst fan-out |
|---|---|---|---|---|
| real complex → complex | 4,924 of 13,880 (35.5%) | (11% of the original severed) | — | — |
| + protein → containing complex | 4,543 (32.7%) | **381 (8%)** | **7,821** | **840** (Cellular Senescence) |

Recovery is concentrated — Mitotic G1 420 → 164, Base Excision Repair 90 → 46,
WNT 323 → 282, R-HSA-453274 150 → 121 — and zero in most bundles. Fan-out 840
is the hub-flooding shape that cost −15pp; WNT and HDR emit 4,051 and 1,160
edges through shard multiplicity. **Not built as a blanket class.** If it is
ever tried, it is with a fan-out cap (Mitotic G1's recovery lives at fan-out
1–2), and the silo record's warning applies: do not sweep the cap on the
evaluation set.

**Where the two composition classes land:** ~19% of severed curator routes
recovered between them; **~33% of all curator routes remain cut** at the uuid
level. What those die on (classified sweep-wide below) is the next generator
defect — the four bundles already classified say released-subunit sinks (the
broadcast trade-off, left alone) and `simple_entity` / `simple_complex` shard
splits, i.e. instance multiplicity: one product node per variant reaction
(HDR's 33 BCDX2 copies; Base Excision Repair's 77 of 90). That is the 737-case
bucket from the gap anatomy and an LNG matcher change, not an edge to add.

### What the remaining cut routes die on — sweep-wide, 24 bundles

| frontier at the first uuid break | routes | share |
|---|---|---|
| **released-subunit `dissociation_sink`** | **3,643** | **74.0%** |
| `simple_entity` shard split (instance multiplicity of a protein) | 649 | 13.2% |
| `simple_complex` shard split (set-member shattering / complex copies) | 330 | 6.7% |
| no stable-id path in the bundle | 189 | 3.8% |
| `reaction` variant split | 101 | 2.1% |

Small molecules are already excluded as carriers, so the 3,643 are *protein*
sinks: the class whose reconnection measured −15pp (shared node), −73 and −77
(capped silo bridge), −204pp (full merge). A bounded bridge was considered and
rejected on paper: a bound tight enough to stop false change (≤ the ±15%
classification band) also stops the bridge from ever moving a call, so for a
3-class score it carries nothing. **Left alone, on purpose.**

The ~22% in the three split rows is one defect: **the matcher mints one
product node per variant reaction** (BCDX2 × 33 from 5 reactions; Base
Excision Repair 77 of 90; the ISGF3:KPNA1 member never materialised). That is
the gap anatomy's 737-case instance-multiplicity bucket, now located at the
node level and sized at ~1,080 curator routes across these 24 bundles. It is
an LNG matcher/decomposition change — share the product node across variant
reactions producing the same base entity — not an edge to add, and it is the
identified next construction fix.

### Where the construction thread lands (pending the A/B)

| | severed curator routes recovered | built? |
|---|---|---|
| complex → containing complex (`LNG_COMPOSITION_EDGES`) | 11% | yes, default off, validated on Interferon, A/B running |
| protein → containing complex | +8%, fan-out to 840 | no — hub risk |
| reconnect released subunits | (74% of the remainder) | no — measured harmful ×4 |
| one product node per entity, not per variant reaction | ~22% of the remainder | **next** |

**The connectivity ladder now documents both states.** Its Tier 4
(`test_cytosolic_isgf3_reaches_the_readout`) was `xfail(strict)`; it is now
conditional on the bundle: with no composition edges it asserts the readout is
*not* reached and xfails with "regenerate with LNG_COMPOSITION_EDGES=1"; with
them it must pass. Measured: old catalog 12 passed / 3 xfailed; regenerated
Interferon **13 passed / 2 xfailed**. LNG's full suite: 994 passed (the one
deselection is the pre-existing Class I MHC edge-ratio test). The oracle is
promoted to `bench/analysis/curator_oracle.py` with its graph functions
unit-tested without Neo4j.

### A risk surfaced by the full regeneration: shard multiplicity

Most regenerated bundles emit 14–68 composition edges at one edge per relation.
HDR (R-HSA-5693567) emitted **747 edges over 78 relations, max 33 per
relation**: the BCDX2 complex (R-HSA-5685316) exists as **33 node copies
produced by only 5 reactions**, each a full copy (out=11) — the Hungarian
matcher mints one product node per variant reaction. That is the
instance-multiplicity defect from the gap anatomy (737 cases), given a face.
All 33 feed the same container as `and`, so under limiting-reactant aggregation
the container is capped by the *least* of 33 copies of one entity — a
perturbation reaching one variant reaction collapses it while 32 copies sit at
baseline.

The right fix is one BCDX2 node (an LNG defect, separate). Until then the
least-wrong aggregation for same-entity shards is an empirical question, so it
is measured as its own arm rather than argued: **`DS_COMPOSITION_GROUP`**
(default off) groups composition inputs by the source's base stable id — max
within a group (node copies of one entity are alternatives), then the
limiting-reactant min across groups (distinct components are co-required).
Threaded through `Reaction.activator_group` → dedup compaction →
`IndexedReaction` → the assembly branch of the propagator. Pinned in
`test_propagator_invariants.jl` (82 assertions): default byte-identical; one
copy knocked out collapses the container under default and leaves it at
baseline under grouping; a distinct component knocked out still collapses it.

**The A/B is therefore three arms on one regenerated catalog**, production
solver defaults, held-out split, both ground-truth axes, concentration
columns and per-pathway net:

| arm | client `DS_SKIP_EDGE_TYPES` | server `DS_COMPOSITION_GROUP` | measures |
|---|---|---|---|
| `comp_ctrl` | `composition` | off | today's behaviour on the new catalog (zero label churn vs the other arms) |
| `comp_live` | — | off | composition edges, min over copies |
| `comp_group` | — | **on** | composition edges, max within entity / min across |

### Catalog-scale reachability — the two traced pathways were not representative

The promoted oracle run over the first 16 regenerated bundles against their
originals:

| | curator pairs | severed before | severed after | recovered |
|---|---|---|---|---|
| 16 bundles, small molecules allowed as carriers (superseded) | 11,608 | 4,981 (42.9%) | 4,532 (39.0%) | 449 (9%) |
| **19 bundles, corrected oracle** | **11,194** | **4,794 (42.8%)** | **4,261 (38.1%)** | **533 (11%)** |

The correction barely moved the totals, so the GTP/ATP sinks were not most of
the severing. Notable recoveries: Signaling by WNT 487 → 323 (1,590 edges,
fan-out max 32 — shard multiplicity again), Rho GTPases 3,696 → 3,378,
Cellular Senescence 158 → 127. Zero recovery despite edges: RAF, DDX58/IFIH1,
Base Excision Repair, HDR.

Interferon recovered 87% of its severed routes; catalog-wide the same edges
recover **11%**. Several bundles recovered nothing despite 56–94 edges (RAF
58→58, R-HSA-168928 207→207, R-HSA-73884 90→90); most of the recovery is one
pathway (R-HSA-194315, −408 of 449); HDR's 747 edges touch zero routes. The
two pathways traced by hand were chosen because they were the worst no-path
pathways — and they are where the composition-hierarchy gap *dominates*, not a
sample of it. Recorded so the A/B is read with the right expectation: a modest
effect, concentrated where the hierarchy gap is the mechanism.

**Caveat on the table above, found by classifying the still-severed routes:
the 42.9% is inflated by an oracle artifact and is being recomputed.** In
Rho GTPases **3,384 of 3,984** still-severed routes, and in DDX58/IFIH1 188 of
207, die on a dissociation sink whose entity is **R-ALL-29438 or R-ALL-113592 —
GTP and ATP**. The oracle excluded small molecules as *roots* but let them
*carry* signal mid-path. A curator does not route a perturbation through GTP,
and the solver pins cofactors inert, so those are not curator routes.
`drop_carriers` now removes non-protein nodes and declared cofactors from the
oracle graph entirely (unit-tested); the sweep is being rerun and the numbers
above will be replaced.

What survives the correction is real and of two kinds. **Base Excision Repair:
77 of 90** severed routes die on a `simple_entity` shard split — a protein
existing as several positional copies with the arriving copy disconnected from
the departing one: instance multiplicity, the same defect as HDR's 33 BCDX2
copies. **RAF: 14** die on `simple_complex` splits, the Interferon-shaped
set-member shattering. Those two are the next generator defects; released
protein sinks (the broadcast class) are the known trade-off and are left.

**Reachability is necessary, not sufficient.** A full-catalog regeneration
with the flag on is running; the A/B is then on **one** catalog —
`DS_SKIP_EDGE_TYPES=composition` versus live — so label churn is zero, at
production solver defaults so it measures the edges alone, held-out split,
both ground-truth axes, concentration columns.

**Interferon α/β at one hop: it is not** — 392 → 385. Classifying the still-severed
pairs by the first uuid break: **260 of 385 break at one reaction**,
`R-HSA-9710959` *p-STAT1 dimer binds KPNA1* → p-Y690-STAT2 (frontier kind
`reaction`); 78 at the ISGF3 dissociation sink; 41 have no stable-id path
either.

The reaction's input is a **set whose members are complexes** — "p-STAT1
dimer, ISGF3" = {ISGF3, p-STAT1-1 dimer} — plus KPNA1; its output is the
matching set {ISGF3:KPNA1, p-STAT1:KPNA1}. `_matching_leaves` decomposes a
complex *if it contains an EntitySet*; ISGF3 contains the STAT1 isoform set,
so the ISGF3 member is shattered to IRF9/STAT2/STAT1 leaves at the matching
layer, and the product **ISGF3:KPNA1 (R-HSA-9710958) is never materialised —
zero nodes, including variants.** ISGF3 (live, out=4) and ISGF3:KPNA1:KPNB1
(live, in=6) both exist; the chain has a missing middle. That is also why the
hierarchy simulation could not help here: there was no node to connect to.
Traced to the node: the live ISGF3 node is fed by *ISGF3 formation*
(R-HSA-909725) and its **only outgoing edges are four dissociation sinks**
(IRF9, STAT2, STAT1-1, STAT1-2). In LNG's network ISGF3 is a *terminal*; in
Reactome it continues, via the set "p-STAT1 dimer, ISGF3", into *binds KPNA1*.
Because no variant reaction consumes the ISGF3 complex node (the set member was
matched at leaf granularity), boundary expansion sees a terminal complex and
decomposes it to sinks — the two defects compound. **The entire nuclear branch
of interferon signalling is severed at one node.**

A two-hop hierarchy edge (ISGF3 → its grand-container) would bridge it at
fan-out 1 (sized below); the real fix is in LNG: when a set member is a
complex, the variant reaction should consume the complex *node* and produce
the product complex *node*, not their leaves.

### The regulator hypothesis — REFUTED, and a defect in the oracle fixed

The oracle first reported `PTPN6,PTPN11 [cytosol]` (a `DefinedSet` negative
regulator of *Phosphorylation of STAT1*) as absent from LNG's network. That was
the oracle's error: it matched node identity by exact stable id and so missed
`::variant::` nodes. A corrected catalog-wide check — every `Regulation` under
the 92 pathways, presence judged by a regulator-class edge from the entity or
any member into the reaction's node — finds **all carried**: 61 set-typed
regulators (29 negative, 32 positive), 742 complex-typed, 137 protein-typed,
24 small-molecule. **No set-typed regulator is dropped.** Adam's hypothesis
was worth one query; it does not hold. (An earlier version of this section
claimed one confirmed instance; that claim is withdrawn.)

### Predictions for the composition A/B, stated before results

Baseline for every comparison is `comp_ctrl` (same catalog, edges skipped
client-side), so label churn between arms is zero by construction.

12. **`comp_live` held-out net is positive but modest** — on the order of
    +20 to +60, not the hundreds the two traced pathways would suggest — and
    concentrated in the pathways whose severed routes the oracle showed
    recovering: Interferon α/β, Mitotic G1, Rho GTPases. Pathways where the
    sweep recovered nothing (RAF, DDX58/IFIH1, Base Excision Repair) move
    little either way.
13. **No broadcast signature:** catalog-wide false change (truth NORMAL,
    predicted a change) does not rise by more than ~1% in `comp_live`, and no
    held-out pathway loses more than ~10 net. Fan-out 1–2 edges should not
    reproduce the −15pp leaf-bridge failure; if false change jumps, the
    shard-multiplicity `and`-min (HDR's 33 copies) is the first suspect.
14. **`comp_group` ≥ `comp_live` where shard multiplicity is high** (HDR,
    R-HSA-194315) and ≈ elsewhere; if `comp_group` is worse overall, min over
    copies was the better-than-expected semantics and the grouping is dropped.
15. Both arms are reported with the concentration columns and per-pathway net;
    a gain that is one pathway or one perturbation is reported as such.

**A/B coverage note, found when the control arm scored 22,324 cases instead of
23,908:** DNA_Double-Strand_Break_Repair (1,584 held-out cases) is absent from
all three arms. Its regenerated bundle is complete, but 8,234 composition
edges took it from ~15k to **23,388 edges**, over the benchmark client's
default `--max-edges 20000`, which skips a pathway silently. The skip is on the
raw file, so it applies to control and live arms alike: the 80-pathway paired
comparison is valid and DSB Repair is unmeasured, not broken. Supplementary
DSB-only arms with the cap raised are chained behind the main three; the
merged result is what gets reported. (8,234 edges in one pathway is shard
multiplicity again — one more reason the instance-multiplicity fix is next.)

### Results — `comp_live` vs `comp_ctrl`, 80 pathways (DSB Repair merged in below)

| | comp_ctrl | comp_live |
|---|---|---|
| accuracy / macro-F1 (22,324 cases) | 83.74% / 0.8015 | **83.09% / 0.7874** |
| ALL net (fixed/broke, p) | | **−146** (316/462, <0.0001) |
| held-out net; pathways / readouts / genes | | **−83** (168/251, 0.0001); 40 / 169 / 150 |
| tuning net | | −63 (148/211) |
| false change (truth NORMAL → predicted change) | 1,161 | **966 (−195)** |
| largest per-pathway nets | | **Interferon α/β +74**, WNT +20 · TP53 −47, IFN-γ −33, DAP12 −30, BER −18, MHC −17 |

**P12 — half right, and the half that matters fails.** The gain is exactly
where the oracle said it would be — Interferon α/β +74 is the largest single
move, WNT +20 next — but the held-out net is **−83**, not positive: the edges
break more than they fix, broadly (40 pathways, 150 genes), by **missed
change** (316 fixed against 462 broke while false change *fell*).
**P13 — holds.** No broadcast signature: false change dropped 195 catalog-wide.
So the loss is not the leaf-bridge failure repeating; it is a container being
**capped by its least copy** — 31,159 composition edges, most of them from
node copies of one entity, all `and`-min'd — suppressing genuine signal. That
is precisely what `comp_group` (max within an entity's copies) was built to
test; its result decides P14.

**As built, `LNG_COMPOSITION_EDGES` is a net loss and stays default-off.**
Whether the *idea* survives depends on `comp_group`.

### DSB Repair (1,584 cases, run separately with the edge cap raised)

| | ctrl | live | **group** |
|---|---|---|---|
| accuracy | 79.23% | 79.10% (−2) | **71.84% (−117)** |
| macro-F1 | 0.7969 | 0.7843 | **0.7180** |
| false change | 269 | 245 (−24) | **330 (+61)** |

**P14 fails on DSB, decisively, and the way it fails settles the semantics.**
The 296 breaks under grouping are **226 NORM → UP on over-expressions** (one
elevated copy lifts the container) and **70 DOWN → NORM on knockouts** (an
untouched copy holds it at baseline). Under min-over-copies (`live`) the
failure is the mirror image — missed change. Neither aggregation over copies
can be right, because the copies are neither alternatives nor co-required: they
are **the same entity, minted once per variant reaction**. The correct
treatment is one node, and that is a generation-time fix with a precise scope —
variants of *one* reaction share one product node per base entity — not a
solver aggregation rule. `DS_COMPOSITION_GROUP` is retained as the measurement
that refuted its own hypothesis.

### Final A/B — all 81 pathways, 23,908 cases, one catalog, zero label churn

| | comp_ctrl | comp_live | comp_group |
|---|---|---|---|
| ALL net (fixed/broke) | — | **−148** (480/628) | **−507** (639/1,146) |
| held-out net (fixed/broke, p) | — | **−85** (332/417, 0.0021) | **−322** (375/697, <0.0001) |
| held-out pathways / readouts / genes moved | — | 41 / 186 / 185 | 41 / 244 / 199 |
| tuning net | — | −63 | −185 |
| false change (truth NORMAL → change) | 1,430 | **1,211 (−219)** | 1,325 (−105) |
| DSB Repair alone (1,584) | 79.23% | 79.10% | **71.84%** |
| largest per-pathway nets (live) | | Interferon α/β **+74**, WNT +20 · TP53 −47, IFN-γ −33, DAP12 −30 | |

**Scored against the pre-registered predictions:**
- **P12 — fails.** The gain landed exactly where the oracle said (Interferon α/β
  +74, the largest single move; WNT +20), but held-out is −85, distributed over
  41 pathways and 185 genes: a genuinely distributed *loss*.
- **P13 — holds.** False change fell 219. The leaf-bridge failure did not repeat.
- **P14 — fails decisively.** Grouping (max within an entity's copies) is worse
  everywhere: −322 held-out, −117 on DSB alone. The DSB breakdown shows why:
  226 false UPs from over-expressions (one elevated copy lifts the container)
  and 70 missed DOWNs from knockouts (an untouched copy holds it). Min over
  copies fails in the mirror direction. **No aggregation over copies is right.**
- **P15 — satisfied**: concentration and per-pathway net reported; the result
  is not a single-pathway or single-perturbation artifact in either direction.

**Verdict.** `LNG_COMPOSITION_EDGES` as built is a net loss and stays default
off. The *mechanism* it targets is real (Interferon +74 on 15 edges), but the
31,159 edges it emits are mostly duplicates of one entity — one node per
variant reaction — and a container fed by N copies of one molecule cannot be
aggregated correctly by any rule. The defect is upstream, in how many nodes a
molecule gets, not in the edge. `DS_COMPOSITION_GROUP` is retained as the
measurement that refuted its own hypothesis and is not recommended.

What the A/B therefore says to do: collapse the copies at generation time
(`LNG_SHARE_VARIANT_NODES`, implemented above), *then* re-measure composition
edges on the deduplicated catalog. Validation of the collapse on HDR follows.

### The instance-multiplicity fix — implemented, default off

`LNG_SHARE_VARIANT_NODES=1`: in Phase 1 of UUID assignment (`_register_phase1`,
factored out of `create_pathway_logic_network` so it can be tested), the same
entity in the same role across the **variants of one Reactome reaction** shares
one UUID. Variants that differ in *which* set member they use have different
entity ids and are untouched; copies across *different* reactions stay
positional, as the silo record requires; boundary entities keep their own
per-stId caches. 7 unit tests (`tests/test_variant_node_sharing.py`); the
existing uuid-position and reaction-connection suites are unchanged (133 pass
together). Validation on HDR — BCDX2's 33 copies, and the 747 composition
edges they generated — is running, alone and combined with composition edges.

### Sharing validated on HDR, and predictions for the sharing A/B (stated before results)

Regenerating HDR with `LNG_SHARE_VARIANT_NODES=1` (1,441 (entity, reaction,
role) registrations reused a sibling's uuid; 197 distinct):

| | cat_os | share | share + composition |
|---|---|---|---|
| BCDX2 (R-HSA-5685316) nodes | 33 | **5** | 5 |
| nodes | 1,707 | **787** | 787 |
| edges | 3,648 | 3,123 | 3,273 |
| composition edges | — | 0 | **150** (was 747) |
| oracle uuid-connected | 507/507 | 507/507 | 507/507 |

The collapse does what it was built to do: 54% fewer nodes, the composition
catalog on HDR shrinks 5x, and no root→terminal route is lost. Full-catalog
regeneration of both arms launched 2026-09-19.

**Pre-registered predictions.** Three arms on the 81-pathway benchmark, one
catalog per arm (sharing changes uuids, so this is necessarily a cross-catalog
comparison against `cat_os`; the composition contrast within the shared catalog
is single-catalog via `DS_SKIP_EDGE_TYPES`).

- **P16 (structure)** — catalog-wide node count falls by more than a third and
  the composition-edge count falls from 31,159 to under 10,000. Oracle
  uuid-severed fraction on the 19-bundle sweep does not rise.
- **P17 (share vs ctrl)** — held-out net ≥ 0 and false change does not rise:
  sharing adds no routes, it only merges copies. The honest risk, stated now:
  a shared product node now has N producing variant reactions combined by
  OR = mean, so a knockout that empties one variant is diluted by N−1
  untouched siblings. If P17 fails by *missed change* (not false change), the
  diagnosis is OR-mean dilution over variant reactions, the fix is the OR rule
  for same-Reactome-reaction variants, and sharing itself is not the defect.
- **P18 (share+composition vs share)** — held-out net > 0, Interferon α/β
  keeps most of its +74, and DSB Repair moves by less than ±10: with one node
  per entity the composition edge no longer feeds a container from duplicated
  copies, which is the mechanism that cost −85 in P12.
- **P19 (concentration)** — pathways moved / distinct readouts / distinct
  genes / dominant-perturbation share and per-pathway net reported for every
  contrast; a result carried by one pathway or one gene is not adopted.

Decision rule: adopt `LNG_SHARE_VARIANT_NODES=1` as the generator default only
if P17 holds; adopt composition edges only if P18 holds on top of it.
Experimental axis reported alongside as always.

### HDR alone, three arms (360 scored cases) — early signal, tuning pathway

| arm | accuracy | macro-F1 | vs ctrl |
|---|---|---|---|
| ctrl (`cat_os`) | 262/360 = 72.78% | 0.6485 | — |
| share | 262/360 = 72.78% | 0.6485 | **0 discordant predictions** of 360 |
| share + composition | 257/360 = 71.39% | 0.6523 | −5 (27 fixed / 32 broke), 10 readouts, 6 genes |

Collapsing BCDX2 from 33 nodes to 5 (and the pathway from 1,707 nodes to 787)
changed **no prediction** on HDR, although 120 cases resolved a smaller
perturbation set. Composition edges on the deduplicated HDR are mildly
negative (−5, p 0.60; earlier DSB Repair with composition on the duplicated
catalog was 79.23 → 79.10, and grouped 71.84). HDR is one of the paper's ten,
so this is a tuning-set signal only; the catalog-wide arms decide P16–P19.

### P16 scored: structure and reachability on the shared catalogs (all 92 bundles)

| | cat_os | share | share + composition |
|---|---|---|---|
| nodes | 107,893 | **70,402 (−35%)** | 70,402 |
| edges | 285,497 | 224,189 (−21%) | 229,814 |
| composition edges (CSV) | — | 0 | **5,625** (was 31,159) |
| shared (entity, reaction, role) registrations | — | 58,052 → 11,278 uuids | same |
| oracle root→terminal pairs | 41,520 | 41,520 | 41,520 |
| uuid-severed | 16,753 (40.3%) | **16,753 (40.3%)** | **12,233 (29.5%)** |

**P16 holds on every clause.** Sharing removes a third of the nodes and
five-sixths of the composition edges, and changes reachability in *no*
pathway — the severed count is identical in all 92 bundles, which is what
merging copies should do (a copy's routes are its siblings' routes). On the
deduplicated catalog the composition edges recover **27% of severed routes**
(4,520 of 16,753) against 11% on the duplicated one, and the recoveries are
where the oracle first pointed: Interferon α/β 392 → 52, DNA Repair (R-HSA-73894)
4,791 → 2,007, FOXO (R-HSA-8878171) 570 → 165, R-HSA-5617472 214 → 0, DSB Repair
50 → 33. Reachability is necessary, not sufficient; the A/B follows.

### Sharing A/B — results against P17–P19 (82 pathways, 24,100 cases; DSB Repair and CD28 now in scope with the edge cap at 40,000)

Curator axis, unconditioned on all 23,908 shared keys (`el_analysis.py`) and
conditioned held-out (`holdout_report.py`):

| contrast | ALL net | held-out | tuning | false change | pathways moved |
|---|---|---|---|---|---|
| **share vs ctrl** (cross-catalog) | **−86** (20/106) | **−8** (0/8; STAT3 in MET is 75%) | −78 | 1,430 → 1,437 (+7) | **5 of 82** |
| **share+comp vs share** (same catalog) | −2 (528/530) | **−79** (340/419, p 0.0046, 42/71 pw, 187 genes) | +77 | 1,437 → 1,210 (**−227**) | 52 |

Experimental axis (849 cases, all in the tuning ten): share vs os **+12**
(17/5, p 0.017; PIP3 +9, TP53 +3); share+comp vs share −8 (15/23, n.s.).

Per-pathway, share vs ctrl: **TP53 −92**, PIP3 +14, MET −6, EGFR −1, RUNX2 −1;
the other 77 pathways are prediction-identical. HDR, where BCDX2 collapsed
33 → 5, moved nothing.

**P17 — fails narrowly, and not by the predicted mechanism.** Held-out is −8
on one gene (STAT3 KO in MET: six NORM→UP false changes), false change is flat
(+7), and the experimental axis is +12. The pre-registered failure mode —
OR-mean dilution producing *missed* change — did not appear: the one large
move is a **sign inversion** in TP53: AKT1 KO and AKT2 KO (49 cases each) →
TP53 targets, truth UP (AKT phosphorylates MDM2; less AKT, more p53), control
predicts UP at 2.6x, sharing predicts DOWN at 0.49x. Node-level trace of the
AKT1 route is in progress; it is a mechanism, not noise, and it is one
pathway. PIP3 +14 is the mirror: PTEN OE → truth DOWN, control said UP,
sharing says DOWN.

**P18 — fails, and the failure identifies the real defect.** Dedup did not
change the composition edge's effect at all: −79 held-out against −85 on the
duplicated catalog, and the per-pathway list is the same to the case
(Interferon α/β +74, IFN-γ −33, DAP12 −30, BER −18, MHC −17, VEGF −12). So
the loss was never about copies. The transitions say what it is:

| pathway | discordant | dominant transition |
|---|---|---|
| IFN-γ | 35 | **DOWN→NORM, truth DOWN, 34** (PTPN11 OE, SOCS1 OE) |
| DAP12 | 30 | **DOWN→NORM, truth DOWN, 30** (six KOs) |
| IFN-α/β | 76 | NORM→DOWN, truth DOWN, 75 (JAK1 KO, IFNAR2 KO, PTPN11 OE) |
| TP53 | 176 | UP→NORM truth NORM 52 fixed; UP→DOWN truth DOWN 35 fixed; DOWN→NORM truth DOWN 15 broke |

Composition edges *mask* an existing DOWN. The cause is in
`compute_reaction_output_vec`: when a node has both an AND cluster and an OR
cluster, the default combination is **`max(and_result, or_result)`**. A
container complex's producing reaction is its OR cluster; the new composition
edge is an AND input from the component complex. When a perturbation drives
the producing reaction DOWN but does not touch the component, the component
sits at baseline and `max` returns baseline: the DOWN is discarded. Where the
container had *no* producing route (the severed Interferon α/β branch) the
composition edge is the only input and it works. The `DS_OR_COMBINE=gate`
alternative multiplies the AND result by the OR cluster's fold-change — which
is the semantics a composition edge needs: *container = component × producing
fold*. `specs/007` measured `gate` as changing **0 of 21,450** predictions on
the network without composition edges, so the attribution is clean.

**P19 — satisfied**: concentration reported for every contrast above.

**Pre-registered P20 (`DS_OR_COMBINE=gate` on the shared catalog, three arms):**
`share_gate` vs `share` changes 0 predictions (replicating specs/007 on this
catalog); `shareplus_gate` vs `share` is held-out > 0 with IFN-γ and DAP12's
DOWN→NORM cases restored and Interferon α/β keeping its +74; if instead the
gate broadcasts (false change rises catalog-wide), composition edges are
refuted as an edge and the hierarchy gap needs a different carrier.

### The TP53 −92 traced to a node: it is the loop knife-edge, not sharing

One case, AKT1 KO → TIGAR (truth UP), solved on both networks with node
activities read back (`trace_akt1.py`, `explain_node.py`; folds vs baseline):

| node | control (`cat_os`) | shared |
|---|---|---|
| p-AKT1 (pinned) | 0 | 0 |
| AKT phosphorylates MDM2 → p-MDM2 [cytosol] | 0.5 | 0.5 |
| p-MDM2 [nucleoplasm], the translocation product | **0.05** | **5.0** |
| MDM2 dimers | 0.19 – **100** (railed) | 1.9 – 5.0 |
| MDM2:TP53 (depletes TP53) | 0.27 | 2.65 |
| TP53 → TIGAR | **2.24 UP ✓** | **0.77 DOWN ✗** |
| solver | **converged = false**, 500 iterations | converged, 417 |

The nuclear p-MDM2 node has one producer (translocation, 0.5) and three
**depletion** inhibitors — its consumers, the MDM2 dimer nodes and MDM2:TP53.
Two of those dimer nodes are fed by "USP7 deubiquitinates MDM2", which sits on
a cycle: MDM2 dimer →(catalyst *and* substrate)→ MDM2 autoubiquitination →
PolyUb-MDM2 → "DAXX binds Ub-MDM2 and USP7" → DAXX:Ub-MDM2:USP7 → USP7
deubiquitinates → MDM2 dimer. A catalytic recycling loop (Type II in
`specs/008`), positive throughout, AND-multiplicative, so **gain exactly 1 at
baseline** — the knife-edge of `specs/014`. In the control network the loop
rails to 100x (the solve never converges); the railed dimers *suppress* nuclear
MDM2 through depletion (floored at 0.1), MDM2:TP53 falls, TP53 rises — the
right readout from a runaway. In the shared network the loop falls to the
all-zero root (0 × anything = 0 is absorbing); the zero-valued dimers
*de-repress* nuclear MDM2 10x, MDM2:TP53 rises, TP53 falls — the wrong readout
from a converged solve. Same loop, two basins, decided by topology details
the merge changed. Neither intermediate state is biology.

So the −92 is not a sharing defect and not something to fix in the generator:
sharing is structurally right (P16), prediction-identical in 77 of 82
pathways, +12 on the experimental axis, and its one large loss is the
already-root-caused loop coin-flip landing on the other side in one tuning
pathway. The designed remedy is loop elasticity (`specs/014`, ε on
cycle-closing activator edges so the loop gain is < 1 and the loop relaxes to
its external drive). With the loop relaxed, the expected chain is nuclear MDM2
≈ 0.5 → dimers ≈ 0.5 → MDM2:TP53 ≈ 0.5 → TP53 ≈ 2 → TIGAR UP.

**Pre-registered P21 (`share_el`: shared catalog + ε_lo 0.5 / w 0.2 / ε_hi
0.95, the clean config from specs/014, composition skipped) vs `share`:**
the 98 AKT1/AKT2-KO TP53 cases return to UP; held-out ≥ +20 (it was +29 on
`cat_os`); false change does not rise; concentration reported. If TP53 does
not return, the loop is landing in a third state and elasticity needs
re-examination on this catalog before anything else.

### P20 scored: `DS_OR_COMBINE=gate` fixes the masking and breaks the loops

| contrast | ALL | held-out | false change | largest nets |
|---|---|---|---|---|
| share_gate vs share | **0 / 24,100 changed** | 0 | 1,437 → 1,437 | — (specs/007 replicated on this catalog) |
| share+comp+gate vs share | −145 | **−79** (146/225) | 1,437 → **1,761 (+324)** | IFN α/β **+100**, DSB **−151**, TP53 −53, HDR −16; the 12 non-loop pathways that moved net **+70** |
| share+comp+gate vs share+comp | −143 | 0 (350/350) | 1,210 → 1,761 (+551) | IFN-γ **+33**, DAP12 **+30**, IFN α/β +26, MHC +17 · DSB −149, TP53 −141 |

The mechanism claim was exactly right: with the container multiplied by the
component's fold instead of `max`-ed against it, IFN-γ's 34 and DAP12's 30
masked DOWNs come back to the case, and Interferon α/β gains a further 26.
And the pre-registered failure clause fired at the same time: false change
rises 324 catalog-wide, almost all of it in **DSB Repair (−151) and TP53
(−141)** — the two loop-heavy pathways. A multiplicative composition input
inside a positive cycle is one more edge with gain 1, and the loops rail.

**Verdict on composition edges, second measurement.** The edge is right where
the graph is acyclic (+70 over 12 held-out pathways, with the oracle's targets
recovered by name) and destructive where it is not. Its held-out sign is set
by DSB alone. That is not a reason to adopt it and not a reason to delete it:
it is the third time today that the loop knife-edge has decided the sign of an
otherwise-correct structural fix (sharing's TP53 flip; composition's DSB
blow-up under both `max` and `gate`). The loop is upstream of everything.

**Pre-registered P22 (`shareplus_gate_el`: composition live, `gate`, plus the
specs/014 elasticity ε_lo 0.5 / w 0.2 / ε_hi 0.95) vs `share_el`:** if the
knife-edge is what converts composition's gains into DSB/TP53 losses, then
held-out > 0, DSB's false-change count returns to within 20 of `share_el`'s,
and Interferon α/β / IFN-γ / DAP12 keep their gains. If DSB still blows up
with the loops relaxed, the composition edge is broadcasting for a reason
elasticity does not touch, and it is refuted as an edge.

### P21 scored: elasticity replicates its gain, and does not return the TP53 cases

`share_el` vs `share`: **held-out +30** (34 fixed / 4 broke, p < 0.0001,
8 of 71 pathways, 33 readouts), tuning +5, false change 1,437 → 1,382 (−55);
83.20 → 83.35%, macro-F1 0.7959 → 0.7998. The +29 measured on `cat_os`
(specs/014) reproduces on a different catalog, distributed (PTK6 +13, ERBB2
+11, Chromatin +6, MET +5, EGFR +5). That is the second independent
measurement of the clean elasticity config and it is now adoptable on the
evidence.

But the 102 AKT1/AKT2-KO TP53 cases did **not** return to UP: 88 moved
DOWN → **NORM** (truth UP). The same trace with elasticity on:

| node | share | share_el |
|---|---|---|
| p-MDM2 [cytosol] | 0.5 | 0.5 |
| p-MDM2 [nucleoplasm] | 5.0 | **0.78** (translocation gives 0.52; three depletion inhibitors at 0.87–0.92 de-repress it back up) |
| MDM2 dimers | 1.9–5.0 | 0.79–0.85 |
| MDM2:TP53 | 2.65 | 0.89 |
| TP53 → TIGAR | 0.77 DOWN | **1.05 NORM** |
| converged | yes | **no** (500 iterations) |

The loop no longer rails or collapses — every node sits within 20% of
baseline — but a 2x source perturbation arrives as 1.05x. Three attenuators
stack: the ε = 0.5 exponent on every cycle-closing activator edge (the whole
MDM2–TP53 axis is inside the 1,250-node SCC, so the *main line* is damped, not
just the recycle); the depletion inhibitors (the dimers consume nuclear MDM2,
so fewer dimers means more free MDM2 — the model's own mass-action feedback
cancelling half the drop); and OR-mean with untouched sibling producers.
Elasticity is a global damper. It buys +30 by stopping rails and collapses and
pays for it in exactly the cases where the signal has to travel *through* a
loop.

The trace names a narrower lever. The cycle that mis-sets this pathway is a
**recycled catalyst**: the MDM2 dimer is both substrate and catalyst of its own
autoubiquitination, and USP7 regenerates it. A conserved moiety, not an
amplifier. `DS_SCC_BREAK_CATALYST` (2026-06, `steady_state.jl`) freezes a
recycling catalyst at its component-entry value — the level its *external*
supply sets — instead of letting it ride the cycle. It was judged on the
800-case small set (v1, catalyst at baseline, −9) and never measured on the
wide set with the held-out split, which is how two other decisions in this
project were reversed. Known reporting defect (specs/001 T028): with it on, the
final consistency check omits the freeze, so `converged` reads false
regardless — cosmetic for an A/B. One caveat read from the code before the
arm ran: `supply` is `copy(x)` taken when the component's iteration *starts*,
and a component-interior node has not been solved at that point, so for the
MDM2 dimer the "external supply" is in practice **baseline** — this is the v1
semantics for interior catalysts, with the v2 semantics only for catalysts
whose value was set upstream. The prediction below stands, but if it fails by
NORM rather than by a wrong sign, that is the first thing to check.

**Pre-registered P23 (`share_bc`: shared catalog, `DS_SCC_BREAK_CATALYST=1`,
composition skipped) vs `share`:** the AKT-KO TP53 cases return to UP
(dimer ≈ 0.5 with the recycle dissolved → MDM2:TP53 ≈ 0.5 → TP53 ≈ 2); held-out
≥ 0 with false change not rising; if TP53 returns but other loop-heavy
pathways lose (the 2026-06 pattern: ERBB2 gained, TP53's CDK12/Pol-II cycle
went flat), the freeze is right for recycling catalysts and wrong for
signal-carrying ones, and the two need to be told apart at the reaction level.

### P22 scored: the DSB blow-up is not the knife-edge — and the trace names the rule

`shareplus_gate_el` vs `share_el` (composition live, `gate`, elasticity on, vs
elasticity alone): **held-out −118** (136/254), tuning +3, false change 1,382 →
1,722 (+340). **DSB Repair −151 — the identical figure with loops relaxed**,
Interferon α/β +100 held, RUNX1 −38 new. P22's failure clause fired: the
composition edge broadcasts for a reason elasticity does not touch.

DSB's 178 discordant cases: **108 NORM→UP and 56 NORM→DOWN, truth NORM**, from
RAD52 / ERCC1 / ERCC4 / MUS81 / RTEL1 / XRCC5 over-expression and XRCC1 / LIG3
knockout, over 17 readouts. Under `gate` the container is multiplied by the
component's fold, so an 80x component makes an 80x container, and every
container above it in the hasComponent hierarchy in turn — the repair
complexes are deeply nested, and 365 composition edges sit in DSB. That is not
biology: **a complex cannot exceed its scarcest component**, so raising one
component leaves a container limited by the others, while removing one does
pull it down. RUNX1's MIR675-OE cases (25 of 46) are the same defect at the
gene-regulation layer.

The rule the trace asks for is a **limiter**, not a route or a multiplier:
composition inputs leave the AND cluster, each fold is capped at 1, and the
smallest multiplies the container's result. A knocked-out or reduced component
lowers the container (DAP12, IFN-γ, the severed Interferon α/β branch — all
DOWN cases), an over-expressed component changes nothing (DSB), and a
producing-route DOWN is never masked (`max` is not involved). Implemented as
`DS_COMPOSITION_MODE=limit` (`reaction_model.jl`; default `assembly` is
byte-identical), with 12 assertions in `test/test_propagator_invariants.jl`
(94 pass) and the mode in the guard-rail set (189 pass).

**Pre-registered P24 (`shareplus_limit`: composition live, `limit`, no gate, no
elasticity) vs `share`:** held-out > 0; Interferon α/β ≈ +100 (all its cases
are DOWN), IFN-γ ≈ +33 and DAP12 ≈ +30 restored; DSB within −56..+10 of `share`
(the 108 false UPs cannot occur; the 56 KO-driven false DOWNs — XRCC1, LIG3 —
may persist, because a limiter honours a knockout, and if the curator calls
those readouts NORMAL that is a hierarchy edge the curator did not walk);
false change not above `share`'s 1,437. If DSB still loses more than 56, the
limiter is pulling containers down through a hierarchy hop that is not
load-bearing and the edge needs a hop limit of 1.

### P23 scored: the catalyst freeze helps TP53 and hurts held-out, and the AKT cases still do not return

`share_bc` (`DS_SCC_BREAK_CATALYST=1`) vs `share`: ALL +22, **held-out −65**
(12/77, 3 pathways: DSB −48, Intrinsic Apoptosis −12, EGFR −5), tuning **+87**
(TP53 **+104**), false change 1,437 → 1,360 (−77). The 2026-06 pattern again,
mirrored: TP53's loops dissolve and its false changes fall, but the freeze reads
signal-carrying catalysts in DSB at baseline. Not adoptable. The 102 AKT-KO
TP53 cases: 84 DOWN → NORM, 4 to UP.

The node-level trace with the freeze on, where nothing rails and nothing
collapses, finally isolates the attenuator:

| node | fold |
|---|---|
| p-MDM2 [cytosol] → translocation | 0.50 |
| **p-MDM2 [nucleoplasm]** | **0.82** — one producer at 0.50, three depletion inhibitors: MDM2:TP53 0.90, MDM2 dimers 0.83, 0.83 → de-repression 1/(0.90·0.83·0.83) = **1.63x** |
| MDM2 dimers | 0.83–0.91 |
| MDM2:TP53 | 0.90 |
| TP53 → TIGAR | 1.03 NORM |

The depleters of nuclear MDM2 are its *own products*: the dimers are made
from it and MDM2:TP53 is made from the dimers. When its supply halves, they
fall with it, and the depletion edge then reads "fewer consumers → less
depletion → more free MDM2", cancelling most of the drop (0.50 → 0.82, and the
remaining 0.18 is halved again at each hop by the same edges). A mass balance
does not do this: the substrate's steady level is set by supply and rate
constants, and a product that is low *because* the substrate is low restores
nothing. The depletion edge is right when the complex moves because of its
**other** partner (TP53 KO → MDM2:TP53 falls → free MDM2 rises, the
de-repression it was built for; specs/011's EGFR:CBL / GRB2 case) and a
double count when the complex moves because of the substrate itself. The
propagator cannot tell the two apart at the node: both arrive as "the
depleter's fold fell".

Two rules were considered and one is testable: (a) read the depleter's fold
*relative* to the substrate's own supply fold (P/X ≈ the partner's fold) —
principled for heterodimers but wrong for a homodimer, whose fold is X², so it
would over-de-repress exactly the MDM2 case; (b) a structural flag: for a
depletion edge whose source is a direct product of a reaction consuming the
target (X → R → P, P ⊣ X), allow suppression (the abundant-complex case) but
not de-repression. (b) loses the legitimate partner-KO de-repression on those
edges; whether curators expect that de-repression more often than they expect
upstream signals to pass through is an empirical question, and an A/B. It is
the next arm after the limiter.

**Pre-registered P25 (`share_dep`: shared catalog, composition skipped,
`DS_DEPLETION_OWN_PRODUCT=suppress_only`) vs `share`.** Implemented as rule (b):
a per-depletion structural flag (`depletion_own_product`, set when the depleter
is reached from the target within two activator hops) and a factor cap of 1 on
those edges; default `full` is byte-identical; 13 assertions in
`test/test_propagator_invariants.jl` (107 pass), mode in the guard-rail set
(191 pass). Predictions: the AKT-KO TP53 cases return to UP (supply 0.50 now
reaches the dimers undiluted → MDM2:TP53 ≈ 0.5 → TP53 ≈ 2); held-out ≥ 0;
false change does not rise; the cost shows up as *missed* change on
partner-knockout cases (a KO of one complex partner no longer frees the
other), and if that cost exceeds the gain the two need separating by
stoichiometry rather than by structure.

### P24 scored: the limiter is the best of three composition semantics, and still not adoptable

`shareplus_limit` vs `share`: ALL −83, **held-out −30** (133/163, p 0.09, 13
pathways), tuning −53, false change 1,437 → 1,656 (+219). Non-loop pathways
**+56**. Interferon α/β **+100**; against `shareplus` (`max`), IFN-γ **+33** and
DAP12 **+30** restored as predicted; DSB **−89**, TP53 −54, RUNX1 −11, RUNX2 −10.

Across the three semantics, held-out for composition edges reads −79 (`max`),
−79 (`gate`), **−30 (`limit`)**, with the acyclic gains constant (+56..+70 over
12–14 pathways) and the loop-heavy losses shrinking. The remaining loss is not
what P24 allowed for (DSB ≤ −56 from honoured knockouts): DSB's 93 discordant
cases are 67 NORM→DOWN (truth NORM) and 19 UP→DOWN (truth UP), and their genes
are XRCC1/LIG3 **KO** (32) *and* TIMELESS / RAD52 / CHEK1 / KPNA2 / BRCA2
**OE** (53). A limiter cannot lower anything on an over-expression by
itself; the OE cases arrive through **depletion**: the over-expressed complex
drains its free subunits (the abundant-complex rule), and the limiter then
carries every drained subunit up into every container above it. In TP53 the
same edges did something P21 and P23 could not: the 98 AKT-KO cases went
**DOWN → UP, fixed (78)**, because a lowered MDM2 complex now lowers its
containers instead of being masked — while CDKN2A KO/OE produced 89 new false
UPs through the same route.

So the composition edge and the depletion edge are now coupled: the limiter
gives depletion's suppressions a way up the hierarchy. That is the same
double count P25 targets from the other side (a product depleting the
substrate it is built from), which is why the two are measured together
next.

**Pre-registered P26 (`shareplus_limit_dep`: composition live, `limit`,
`DS_DEPLETION_OWN_PRODUCT=suppress_only`) vs `share_dep` and vs `share`:**
if DSB's OE-driven NORM→DOWN cases are own-product depletion carried up by
the limiter, DSB moves to within −40 of `share_dep`, Interferon α/β / IFN-γ /
DAP12 keep their gains, and held-out vs `share` is > 0. If DSB still loses
more than 60 with own-product de-repression off, the drain is from foreign
complexes and the composition edge stays off pending a hop-1 variant.

### P25 scored: inert on its own — the depletion rule only matters once the loop is tame

`share_dep` vs `share`: ALL −18 (4/22), held-out −15 (RUNX2 −10, EGFR −5),
false change +12; **the 102 AKT-KO TP53 cases did not move** (88 still DOWN).
Fails as pre-registered, and for a reason the two earlier traces already
showed: without elasticity or the freeze the MDM2/USP7 loop collapses to zero
(shared) or rails (control) before any depletion arithmetic matters, so
capping own-product de-repression changes nothing there. The rule was
isolated on the wrong base. The measurement that tests it is on the tamed
loop: `share_el` + `suppress_only`.

**Pre-registered P27 (`share_el_dep`: elasticity clean config +
`DS_DEPLETION_OWN_PRODUCT=suppress_only`, composition skipped) vs `share_el`:**
nuclear MDM2 holds its 0.52 supply instead of recovering to 0.78, so the
AKT-KO TP53 cases move toward UP — at least 40 of the 88 DOWN/NORM cases
become UP (the ε = 0.5 damping on the main line still costs some); held-out
≥ 0 against `share_el`; false change does not rise. If fewer than 20 return,
elasticity's damping of the main line is the larger attenuator and the loop
work needs to distinguish recycle edges from main-line edges.

### P26 scored: DSB's loss is not depletion — it is an 80x inhibitor carried up the hierarchy

`shareplus_limit_dep` vs `share_dep`: held-out −5 (133/138, p 0.81), tuning
−60, false change +61; vs `share`: held-out −20 (p 0.26), all-other pathways
**+69 over 19**, **DSB −89 unchanged**, TP53 −63. P26's failure clause fired:
own-product de-repression is not what drains DSB. Structurally, **0 of DSB's
95 composition sources carry a depletion edge** (IFN α/β 0/8, TP53 0/40,
IFN-γ 0/6), so the depletion hypothesis for the limiter's DSB loss was wrong
before the arm ran.

The trace (RAD52 over-expressed to 80x on DSB, limiter on): 1,681 nodes off
baseline, **281 of 365 composition sources moved, and the moved ones sit at
fold 0.000**, chained through 35 nested containers. The origin is the
inhibition rule, not the hierarchy: `divide` gives an 80x inhibitor a factor
of (b+ε)/(x+ε) = 1/80 on every reaction it regulates, hill_sat AND multiplies
two such inputs to 1/6400, and the limiter — doing exactly what it should —
carries the near-zero component into every container built from it, through
every level of DSB's nesting. Under `max` those containers ignored the
component because their producing route sat at baseline (the masking that cost
IFN-γ and DAP12); under `limit` they honour it. TIMELESS at 80x, by contrast,
moves 2 nodes: its cases in the discordant list came in through the same
RAD52-style chain elsewhere, not through TIMELESS itself.

So the limiter's DSB cost is the price of no longer masking, applied to an
inhibitor over-expression the curators call inert for those readouts (truth
NORM on all 67). Whether an 80x negative regulator should zero its targets
is a question about the inhibition form's suppression bound (specs/006, /011
bounded the mirror cases), not about composition edges — and it is the same
bound question three times now: depletion was floored at 0 (fixed, +28),
de-repression is capped at h_max, and inhibition suppression by an
over-expressed regulator is still unbounded below.

### P27 scored: inert — the main-line damping is the larger attenuator

`share_el_dep` vs `share_el`: held-out −3 (0/3), tuning −12 (ERBB2 −12),
false change +11; **the AKT-KO TP53 cases: 94 NORM, unchanged.** Failure
clause: with the loop tamed by elasticity, removing own-product de-repression
does not let the 0.50 supply through either, so the ε = 0.5 exponent on the
main MDM2–TP53 axis (every edge of it is "in loop" inside the 1,250-node SCC)
is what flattens the signal. Elasticity keeps its +30 by stopping rails and
collapses; the cost is exactly the cases where a signal must *travel through*
a giant SCC. Distinguishing recycle edges (the USP7 cycle) from main-line
edges (AKT → MDM2 → TP53) inside one SCC is the open loop problem, restated
with a named case.

### Where the day landed — every arm, one table

All on the shared catalog (`LNG_SHARE_VARIANT_NODES=1`), 82 pathways, 24,100
cases, production defaults unless stated, held-out = the 70 non-paper
pathways. `ctrl` = current generator (`cat_os`), 23,908 cases.

| arm | vs | held-out | tuning | false change | note |
|---|---|---|---|---|---|
| share | ctrl | −8 (0/8, one gene) | −78 (TP53 −92 loop flip) | +7 | exp **+12**; 77/82 pathways identical; **structurally right** |
| share + composition (`max`) | share | **−79** | +77 | −227 | masks DOWNs (IFN-γ, DAP12) |
| share + composition + `gate` | share | −79 | −66 | +324 | OE multiplies up the hierarchy (DSB −151) |
| **share + elasticity** | share | **+30** (34/4, p<1e-4, 8 pw) | +5 | −55 | **replicates specs/014; adoptable** |
| share + comp + gate + el | share_el | −118 | +3 | +340 | DSB −151 unchanged: not the knife-edge |
| share + catalyst freeze | share | −65 | +87 (TP53 +104) | −77 | 2026-06 pattern |
| share + composition **`limit`** (new) | share | −30 (p 0.09) | −53 | +219 | best comp semantics; +56 in 13 acyclic pw; IFN α/β +100 |
| share + own-product dep cap (new) | share | −15 | −3 | +12 | inert until loop is tame |
| share + comp `limit` + dep cap | share | −20 (p 0.26) | −63 | +73 | DSB −89 = 80x inhibitor OE, unbounded |
| share + el + dep cap | share_el | −3 | −12 | +11 | inert; main-line damping dominates |

Adoptable on the evidence: **sharing** (generator default; reachability
identical in 92/92, one loop-flip loss in a tuning pathway, +12 experimental)
and **elasticity ε_lo 0.5 / w 0.2 / ε_hi 0.95** (+29 on `cat_os`, +30 here,
both p < 0.0001, distributed). Together: 83.20 → 83.35% (macro-F1 0.7959 →
0.7998) against `ctrl`'s 83.44 / 0.8011 on a smaller case set — i.e. net
flat on accuracy, with the structural defects fixed and the loop coin-flips
gone. Composition edges stay off: right where the graph is acyclic (+56..+70
across 12–19 pathways every time), and their sign is set by DSB Repair, where
they faithfully carry an unbounded 80x-inhibitor suppression into 99 nested
containers that the curators call inert.

Four mechanisms were traced to a node today and each is recorded above with
its case: the MDM2/USP7 recycling loop at gain 1; own-product depletion
cancelling upstream drops; `max(and, or)` masking DOWNs; and inhibitor
over-expression suppression with no lower bound. The last is the one not yet
measured: `divide` gives an 80x regulator a 1/80 factor per reaction and
hill_sat multiplies them, the mirror of the depletion asymmetry that specs/011
bounded for +28. That is the next arm, and it does not need composition edges
to be measured.

### Pre-registered P28: bound inhibitor suppression (the fourth traced mechanism)

The divide form gives an inhibitor at fold f a factor 1/f per regulated
reaction, with a de-repression ceiling of 10 per reaction and **no floor**:
an 80x over-expressed regulator suppresses to 1/80, two of them under hill_sat
to 1/6400. Depletion had the same asymmetry until specs/011 floored it
(+28 held-out). The per-edge floor already exists (`DS_INHIBITOR_FLOOR`,
default 0.0; scope `loops` default) and was judged once on the 9-pathway
experimental set with scope `loops` — never on the wide set, never with scope
`all`. Three arms on the shared catalog, production defaults otherwise:

- `share_fl01`: floor 0.1, scope all (mirror of `DS_DEPLETION_H_MIN`).
- `share_fl001`: floor 0.01, scope all (an 80x regulator may still suppress
  to 1/80 on one edge, but compounding across edges is bounded per edge).
- `shareplus_limit_fl01`: floor 0.1 + composition `limit`, to see whether DSB's
  −89 under the limiter was the unbounded RAD52 suppression.

Predictions: `share_fl01` vs `share` held-out ≥ 0 with false change *down*
(over-expressed regulators stop zeroing whole branches) and the losses, if
any, as *missed* DOWN on direct inhibitor-OE targets; 0.01 between 0.1 and
the default; `shareplus_limit_fl01` vs `shareplus_limit`: DSB recovers by
≥ 50 of its −89, IFN α/β / IFN-γ / DAP12 keep their gains, and held-out vs
`share` turns positive. If the floor at 0.1 loses held-out by missed change,
suppression by over-expressed regulators is load-bearing and the bound has
to be per-reaction (product) rather than per-edge to mirror the ceiling.

### P28 scored: the suppression bound is load-bearing where it binds, and DSB was never about it

| arm | vs | held-out | tuning | false change | note |
|---|---|---|---|---|---|
| floor 0.1, scope all | share | **−58** (10/68) | +12 (Cell Cycle Checkpoints) | +19 | RUNX1 −33, IFN-γ −17, RUNX2 −8: **missed change** on inhibitor-OE targets |
| floor 0.01, scope all | share | **0 / 24,100 changed** | 0 | 0 | the benchmark's OE is 80x → 1/80 = 0.0125 > 0.01; a per-edge floor below 0.0125 never binds |
| floor 0.1 + composition `limit` | shareplus_limit | −30 | +15 | +19 | **DSB unchanged** (not in the moved list) |

P28 fails on both clauses. Bounding per-edge suppression at 0.1 costs 58
held-out as missed DOWN — the curators *do* expect an over-expressed
repressor to shut its targets — so the 1/f factor is load-bearing, exactly
mirroring the depletion result where bounding *suppression* helped only
because the suppression there was compounding without limit across edges.
And DSB's −89 under the limiter did not move by a single case when the
suppression was floored, so the RAD52 trace's exact zeros are **not** the
1/80 factor. Exact 0.000 through a chain of 35 nested containers, unmoved by
any floor, is the signature of the third mechanism, not the fourth: a cyclic
component collapsing to the all-zero root (0 × anything = 0 is absorbing),
which the limiter then carries into every container above it — the same
thing that flipped TP53 under sharing, in the other loop-heavy pathway. The
limiter does not create the zero; it stops masking it. That is testable
directly: with the loops relaxed, the zeros should not form.

**Pre-registered P29 (`shareplus_limit_el`: composition `limit` + elasticity
clean config) vs `share_el`:** DSB within ±20 of `share_el` (the −89 vanishes
with the collapse), Interferon α/β ≈ +100, IFN-γ / DAP12 gains kept, held-out
> 0 — which would make composition edges adoptable *together with*
elasticity. Node-level check alongside: RAD52 OE on DSB under the same config
should show no composition source at fold 0.000. If DSB still loses ≥ 60,
the zeros come from something elasticity does not tame and the composition
edge stays off.

### The DSB zeros are a double count in the hierarchy, not a loop and not the floor

RAD52 OE on DSB under `limit` + elasticity: **identical** — 1,681 nodes off
baseline, 281 of 365 composition sources at fold 0.000. So the zeros are not
a loop collapse (P29's node-level clause fails before the arm reports) and
not the suppression floor (P28). They originate at "ATM phosphorylates NBN"
and "RNF4 ubiquitinates MDC1" (26 copies) and are carried through DNA
DSB:p-MRN:ATM:KAT5 → …:H2AFX-Nucleosome → …:MDC1 → …:RNF8 — the nested DSB
response complexes.

The count that explains it: **54% of the catalog's composition edges (3,026 of
5,625) connect a component to a container that a reaction already builds
from it** (S → R → T exists through activator edges; DSB 36%, TP53 51%,
IFN α/β 11 of 15). Under `limit` such a container reads the component's fold
twice — once through its producing reaction, once through the hierarchy edge
— so the fold is squared at every level, and DSB's nesting is deep enough to
reach hill_sat's zero. The 4 non-redundant Interferon α/β edges are the ones
that carried its +100: a hop with no reaction is the only kind that adds
information.

Implemented as `DS_COMPOSITION_MODE=limit_novel`: a per-activator structural
flag (`activator_comp_redundant`, two-hop check at index time) and the
redundant edges skipped; 7 assertions (114 pass), guard rails 192. Under it a
redundant edge is byte-identical to no edge, a novel one still limits.

**Pre-registered P30 (`shareplus_limitnovel` vs `share`; and `_el` vs
`share_el`):** DSB within ±20 of the base (the zeros cannot form), Interferon
α/β ≈ +100 kept (its edges are novel), IFN-γ / DAP12 within ±10 of `limit`
(theirs are mostly redundant, so part of that gain may go), held-out > 0
against both bases; false change not above base. If held-out is still
negative, the remaining loss is in the novel hops themselves and the
hierarchy edge is refuted as a carrier on this network.

### P29 scored: limiter + elasticity is the closest to positive yet, and DSB still carries the sign

`shareplus_limit_el` vs `share_el`: held-out −8 (121/129, p 0.66), tuning −11,
false change +204; non-loop pathways **+56**, loop-heavy −75; Interferon α/β
+100, **DSB −62** (was −89 without elasticity), TP53 −18, EGFR −15, RUNX1 −12.
vs `share`: held-out **+22** (141/119, p 0.19, 15 pathways), all-other **+95**;
macro-F1 **0.8011**, the highest of any arm on this catalog (ctrl 0.8011 on
the smaller case set). DSB moved 27 toward zero with the loops relaxed but
kept −62, so P29's DSB clause fails and the node-level clause already had:
the zeros are the hierarchy double count (54% redundant edges), which
`limit_novel` removes. P30 is the arm that decides composition edges.

### P30 scored: composition edges reduced to their irreducible parts

| arm | vs | held-out | tuning | false change | DSB | IFN α/β | TP53 |
|---|---|---|---|---|---|---|---|
| `limit_novel` | share | +12 (114/102, p 0.45) | −29 | +187 | **−62** | +100 | −25 |
| `limit_novel` | `limit` | **+42** (62/20, p<1e-4) | +24 | −32 | +27 | 0 | +29 |
| `limit_novel` + el | share_el | +10 (117/107) | **−98** | +181 | −62 | +100 | **−99** |
| `limit_novel` + el | share | **+40** (143/103, p 0.013, 15 pw) | −93 | +126 | −62 | +100 | −103 |

Skipping the redundant hierarchy edges did what it was built to do (+42
held-out over `limit`, DSB +27, TP53 +29), and DSB then stops at **−62 under
every configuration** — plain, with elasticity, or both. Its 62 cases are now
readable one by one: **XRCC1 KO 16 + LIG3 KO 16** (a knocked-out component
honestly lowers the complexes built from it; the curators call those readouts
NORMAL — a hierarchy hop the curator did not walk), and **RAD52 OE 14, KPNA2
OE 6, BRCA2 OE 5** (an over-expression reaching a container DOWN through a
*novel* hop, i.e. a suppressed component whose container has no producing
reaction in the network — 25 cases, the residual mechanism not traced
today). TP53 under `limit_novel` + elasticity: CDKN2A KO 49 / OE 48 all
false, in **both** directions from one gene — the MDM2:ARF loop's basin
flipping again, the same coin the AKT cases landed on this morning.

**Where composition edges stand.** Their entire held-out gain is Interferon
α/β: +100 on 25 readouts from 4 genes, every case NORM→DOWN fixed, exactly
the branch the oracle found severed at ISGF3 this morning. Against it: DSB
−62 (32 honoured knockouts + 25 OE-through-novel-hop + 5), and, when combined
with elasticity, a TP53 loop flip in the tuning set. By this repo's rule a
gain carried by one pathway is not a distributed improvement and is not
adopted on that basis; it is a *connectivity fix for one known-severed
pathway*, of the same kind as LNG#89's HSP90B1, and is honest to describe as
that. Configuration that carries it with the least damage: `LNG_COMPOSITION_EDGES=1`
+ `DS_COMPOSITION_MODE=limit_novel` (+12 held-out alone, +40 with elasticity).
Off by default in both repos; the decision is recorded here for Adam.

### Final table for 2026-09-19 (additions to the morning table)

| arm | vs | held-out | note |
|---|---|---|---|
| inhibitor floor 0.1 / all | share | −58 | suppression by over-expressed repressors is load-bearing |
| inhibitor floor 0.01 / all | share | 0 changed | never binds: the benchmark's OE is 80x = 1/80 |
| floor 0.1 + comp `limit` | shareplus_limit | −30 | DSB unmoved: not the floor |
| comp `limit` + el | share_el | −8 (p 0.66) | macro-F1 0.8011; DSB −62 |
| comp `limit_novel` | share | +12 | DSB −62, IFN +100 |
| comp `limit_novel` + el | share | **+40** (p 0.013) | one-pathway gain; TP53 −103 tuning |

Adoptable, unchanged from the morning: sharing + elasticity. New knobs, all
byte-identical by default: `DS_COMPOSITION_MODE` {assembly, limit,
limit_novel}, `DS_DEPLETION_OWN_PRODUCT` {full, suppress_only}. Structural
facts established: 54% of composition edges are redundant with a producing
reaction; 0 composition sources are depleted in DSB/IFN/TP53; the benchmark's
80x over-expression means any per-edge inhibitor floor below 1/80 is inert.

## 5b. Adversarial review, 2026-09-19 (three independent reviewers given code and data, not conclusions)

### Generator diff (reviewer B) — findings, verification, action

| # | finding | verified | action |
|---|---|---|---|
| 1 | Composition edges are emitted in parallel with a curated reaction that already joins the pair (`existing` checks direct pairs only). Interferon 10 of 15, Cell Cycle Checkpoints 19/32, MHC 10/22, Senescence 30/87 redundant. | Yes — independently counted from the solver side this afternoon: **54% catalog-wide** (§5, `limit_novel`). | Handled at solve time (`DS_COMPOSITION_MODE=limit_novel`). The right place is the generator (do not emit when `source → VR → target` exists); needs a regeneration and re-A/B, so recorded as the follow-up, not done today. |
| 2 | The edge targets **every positional copy** of the container (input copy, output copy, catalyst, regulator), and the logged fan-out counts stIds, not edges. | Yes — on the shared catalog 718 of 3,526 (component, container) stId pairs carry >1 edge, **up to 28** (DNA Repair), 24 (FOXO), 16 (RAF/MAPK). | Follow-up with #1: target the *produced* copy only (or one canonical copy). The "median 1, max 2" fan-out claim in §5 is stId-level and is now qualified. |
| 3 | `get_containing_complexes` Cypher leaves `entity` unlabelled → full-store scan, 660–800 ms per call vs 1–21 ms labelled; 30–90 min per catalog regeneration. | Reviewer measured on the live DB. | **Fixed**: `(entity:PhysicalEntity)`. |
| 4 | `except Exception: labels = []` swallows connection errors and silently drops a complex's edges, while the next lookup fails hard. | Read the code; correct. | **Fixed**: catches `(IndexError, KeyError)` only. |
| 5 | `rxn is not None` cannot catch a missing reaction id read back from a cached CSV as the string `"nan"`; every such variant would share `(eid, "nan", role)` across reactions. Unreachable today, wrong sentinel. | Correct. | **Fixed**: `rxn not in (None, "", "nan", "None")`, unmapped count returned and warned. |
| 6 | `(source, target)`-only dedupe suppresses a composition edge when a depletion or set_member edge joins the pair. | Correct, rare. | Left; pinned by an existing test. |
| 7 | Tests: `test_flag_default_is_off` was tautological; the boundary test could pass an implementation that ignored the exclusion; no test for a `"nan"` id or a partial map. | Correct. | **Fixed**: four tests added/strengthened; suite 1,015 pass (one pre-existing MHC failure). |

Checked and found correct by the reviewer: the sharing key cannot collapse two entities, roles or reactions; boundary membership is pathway-global so it cannot differ between variants; Phase 2 union-find, the final row dedupe (keyed with `edge_reaction_id`), `node_reaction_context`, `nodes.csv` and the uuid mapping are all safe under shared uuids; both functions are order-invariant in structure; the Cypher direction is right and cannot return the entity itself. One limitation named: containment through an EntitySet (`hasComponent → Set → hasMember → complex`) is not traversed.

### Methodology and numbers (reviewer C) — verified errata; these CORRECT the conclusions above

Every item below was re-derived by me from the dumps after the review; the
review's arithmetic held in every case checked.

**E1. The `ctrl` baseline was not `cat_os`.** `comp_ab.sh` line 6: `CAT=$S/cat_comp`
— the composition regeneration with composition edges skipped client-side, a
complete relabelling of `cat_os` (1 shared uuid of 2,308 in TP53). Against
the real `cat_os` dump (`ab_onesided.tsv`): `share` is **held-out −1 (0/1)**,
tuning **−64**, and differs in **3 pathways** (TP53 124, PIP3 14, RUNX2 12) —
"79 of 82 identical", not 77. The P17 verdict's "−8 held-out on STAT3 in MET"
was the relabelled baseline's own basin flip (`ctrl` vs `cat_os` differ in
22 predictions, +7 held-out, all STAT3/MET and PIK3CA/EGFR); `share`
reproduces `cat_os` on those six cases to seven figures. **Corrected P17:
sharing is prediction-neutral on held-out and costs 64 tuning cases, all the
TP53 loop flip.** The "Where the day landed" table's `ctrl = cat_os` is false.

**E2. The elasticity "replication" is the same 33 cases.** Held-out fixes on
`cat_os` (+29: 33/4) and on the shared catalog (+30: 34/4) overlap **33 of 33**
(the one new fix is a RUNX2 case); all 8 gaining pathways are
prediction-identical between the two catalogs. This is one measurement run
twice, not two. Its tuning cost — **−61 on `cat_os`** — reads as +5 on the
shared catalog only because `share` had already lost the 102 AKT-KO TP53
cases (ctrl 100/102 correct → share 2/102): the cost moved into the baseline.
DOK1 carries 13 of 34 fixes (38%, under the 50% warning). ERBB2 (+11) is a
tuning pathway and should not have been listed among the held-out gains.
**Corrected: elasticity is a single held-out measurement of +29 (7 pathways,
12 genes, p < 0.0001) with a −61 tuning cost, unmeasured on the experimental
axis until today** (arms launched after the review; results appended below).

**E3. The experimental "+12" for sharing is +1 under the repo's own
conditioning** (325 of 849 cases dropped for a changed perturbation set; 6/5,
p 1.0). No arm after `shareplus` had an experimental measurement, in breach of
P19's "experimental axis reported alongside as always".

**E4. `limitnovel_el` vs `share` held-out +40 (p 0.013) is −60 (43/103,
p < 0.0001) outside Interferon α/β**, and `limit` vs `share` −30 is **−130**
outside it. Composition edges are net **negative in every configuration once
IFN α/β is removed**, including the acyclic pathways: the sentence "right
where the graph is acyclic (+56..+70 across 12–19 pathways every time)" is
false; the acyclic net excluding IFN α/β is −44 (`limit`), −30 (`gate`), −44
(`limit_el` vs `share_el`), −31 (`limit_dep`), −5 (`limit_el` vs `share`).
P30's own clauses (held-out > 0 vs `share_el`: +10, p 0.55; false change not
above base: +126/+181) fail and were not scored as failing. IFN α/β scores 21%
under every base (268 of 448 predicted NORMAL against 338 true changes), so
any connectivity at all produces a large raw gain there. **Corrected:
composition edges are an Interferon α/β fix and a loss everywhere else.**

**E5. Headline comparisons used unequal denominators.** On the 23,908 shared
keys: `cat_os` 83.36%, relabelled ctrl 83.44%, `share` 83.09%, `share_el`
83.23%. Sharing + elasticity is **−0.13pp against the real baseline on
identical keys**, not "net flat". CD28 (192 cases, 97.9% correct) entered
only because `--max-edges` was raised to 40,000 for the share arms and adds
+0.11pp to every share-family headline.

**E6. The DSB traces did not perturb the benchmark's node set.** The benchmark
pins every node whose diagram entity or member leaves carry a gene's stable
id (RAD52 in DSB: 14 uuids, ten of them complexes); `trace_generic.py`
resolved by base stable id (4). The AKT1/TP53 trace used the same 2 uuids as
the benchmark and stands. **The DSB mechanism narrative (P26 "80x inhibitor
carried up the hierarchy", P28's "exact zeros are the loop", the double-count
trace) was derived from a perturbation the benchmark never ran** and is
unverified; the A/B numbers (DSB −89 / −62) and the 54% redundancy count are
unaffected. The `limit_novel` gain over `limit` (+42) is real; its stated
mechanism is not established.

**E7. `fl001` changes nothing for a reason knowable before the arm**:
h_k = (b + ε)/(x + ε) with x ≤ 1 is ≥ 0.01 for *any* inhibitor, so a per-edge
floor of 0.01 can never bind; "1/80 = 0.0125 > 0.01" was the wrong explanation.

**E8. Smaller corrections.** The AKT-KO TP53 cases number 51 + 51 = 102, of
which 96 are truth UP and 6 truth DOWN (49 was the discordant count).
Conditioned tuning for `share` vs ctrl is −92, not −78; the 14 PIP3 "mirror
fixes" are all in the dropped set and are not like-for-like. `TUNING_PATHWAYS`
has 11 entries under the label "the paper's ten"; "held-out = the 70" should
read 71 on the shared catalog. `limdep_share.sh` / `eldep_share.sh` did not echo
the defining knob in their `verified:` line (the values sit in BASE and the
arms behaved as their names say, but it was unguarded). Every "pre-registered"
block lives in this one uncommitted file, so the ordering of predictions and
results is not auditable from history; nothing was found that cites a later
number, but that is not the same as proof.

**What survived**: all 17 shared-catalog arms have identical key sets and
perturbation fingerprints; every spot-checked net, p, false-change count and
per-pathway figure reproduces; `share_gate` = `share` byte-identical;
`fl001` = `share` byte-identical; the AKT1 → TP53 basin-flip trace; the
per-arm tables above are numerically right — their *interpretation* is what
E1–E5 correct.

**Corrected standing recommendation.** Nothing measured today is adoptable on
accuracy grounds. Sharing is structurally right (reachability identical in
92/92, one node per entity per reaction) and prediction-neutral on held-out,
with a 64-case tuning cost that is a loop basin flip; it is a generator
correctness change, not an accuracy change, and should be judged as that.
Elasticity is one measurement: +29 held-out / −61 tuning; adopt only if the
experimental axis (running) is not negative and after a second, genuinely
independent test. Composition edges stay off.

### Solver diff (reviewer A) — findings, verification, action

Differential run by the reviewer: **default-config output byte-identical to
HEAD on 288 solves** (6 cyclic fixtures × 8 env configurations) plus the
sample network; all eight assertion test files pass. The claim broke in two
places that no knob gates:

| # | finding | verified | action |
|---|---|---|---|
| 1 | **Perf regression 3.7x time / 3.5x allocation** on the default path: two `Dict`s allocated and iterated on every `compute_reaction_output_vec` call regardless of mode (3,000-node cyclic solve: 0.05 s → 0.19 s). | Re-ran the reviewer's benchmark after the fix: **0.06 s / 0.09 GB**, same checksum — back to HEAD. | **Fixed**: both Dicts allocated only when their mode is on. |
| 2 | `composition` edges join the assembly-limiting cluster **by default** (`edge_type in ("assembly","composition")`), a behaviour change relative to HEAD for any network carrying the new edge type; "byte-identical" holds only for shipped catalogs, which have none. | Correct. | Documented in the code as a deliberate default change; not a knob. |
| 3 | `DS_DEDUP_ACTIVATORS` merged a plain `input` and a `composition` edge from one source with `d_grp = max(...)`, and `group > 0` is read as "is composition", so under `limit*` the merged edge became a ≤1 limiter (S=4 → X=1.0) and the input role vanished. | Reproduced by the reviewer with numbers. | **Fixed**: the group survives only if every duplicate is composition (a hierarchy edge parallel to a plain input is redundant, the plain role wins); test added. |
| 4 | Own-product detection walked **catalyst** edges too, so an enzyme's product counted as its "own product" and `suppress_only` also disabled enzyme–product de-repression (partner KO: 10.0 → 1.0). | Correct; the comment described a narrower relation than the code tested. | **Fixed**: chain adjacency excludes catalyst (and composition) edges; test added. **P25/P27 measured the broader definition**; both were inert, so the verdict stands, but the arm was not what its name said. |
| 5 | Redundancy two-hopped through **composition** edges too, so the "54%" included hierarchy-only paths. | Recounted with hops restricted to substrate/assembly chains: **33% (1,832 of 5,625)**; DSB 16%, TP53 44%, IFN α/β 10 of 15. | **Fixed** in code. **The `limit_novel` arm (P30) skipped the 54% set**, so its +42 over `limit` and DSB +27 belong to the old definition; the narrowed rule is unmeasured. |
| 6 | `DS_LOOP_ELASTICITY<1` with `DS_SCC_SOLVE=0` was a **silent no-op** (SCC ids not computed, so no edge is "in loop"). | Reproduced. | **Fixed**: components are computed whenever elasticity is on; test added. |
| 7 | `DS_SCC_OPTIMIZER`, `DS_GAMMA_ANNEAL*`, `DS_LM_ITERS` validated lazily (only when a cyclic component is reached), so a typo passes on every acyclic pathway. | Correct. | Left (minimise path is a recorded negative result); noted. |
| 8 | `lm_minimize!` returned "converged" when no downhill step could be found. | Correct; harmless because the caller recomputes the residual. | **Fixed**: returns false. |
| 9 | Elasticity and limiter folds divide by the **target's** baseline; correct only while every node's baseline is 0.01 (the parser hardcodes it). | Pre-existing convention shared by depletion. | Noted. |
| 10 | The determinism test's relabelling could collide silently. | No collision for current fixtures. | **Fixed**: relabel asserts node count preserved. |
| T | The "composition ≠ plain input" assertion was vacuous (A=0 makes both sides 0). | Correct. | **Fixed** (A=0.5, strict `<`). Several "default byte-identical to unset" tests compare two runs of the new code and cannot detect a change relative to HEAD — the differential run is what pins that, and it is not in the suite. |

After the fixes: propagator invariants **120**, guard rails 192, loop elasticity
150 (+1 broken), determinism 80, AND curves 71, cycle handling 34 (+2
broken), worked example 9, observation pinning 23 — all pass.

### Experimental axis for elasticity (missing until the review; now measured)

| | conditioned | note |
|---|---|---|
| `share_el` vs `share` | **0** (4/4, p 1.0) | one perturbation (BRCA1 in TP53) is half the discordant set |
| `cat_os` + el vs `cat_os` | **+12** (19/7, p 0.029) | PIP3 +9, TP53 +3; both tuning pathways (the experimental set is all tuning) |

Not negative on either catalog. Elasticity therefore stands as: one held-out
measurement of +29 (7 pathways, 12 genes, p < 0.0001, DOK1 38% of fixes), a
−61 tuning cost on `cat_os` that the shared baseline hid, and a neutral-to-
positive experimental axis. Adoptable *only* after a genuinely independent
held-out test — different cases, not a re-run — which the current benchmark
cannot supply for this config.

## 6. What this does not yet say

Reachability is necessary, not sufficient: a route that exists can still carry
the wrong sign or magnitude, and a route added can manufacture false change.
The oracle's next use is the sign check — for pairs where both agree on
reachability, does the solver's predicted direction match the oracle's signed
path? That is the part of the curator's reasoning that the current no_path
diff does not test.
